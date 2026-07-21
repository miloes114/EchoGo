.echogo_gost <- function(...) gprofiler2::gost(...)

.echogo_git_commit <- function() {
  value <- tryCatch(
    suppressWarnings(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)),
    error = function(e) character()
  )
  if (length(value) == 1L && nzchar(value)) value else NA_character_
}

#' Run g:Profiler enrichment across multiple species (customizable, cached GO depth)
#'
#' @param de_genes character vector of EggNOG-mapped DE gene IDs.
#' @param bg_genes character vector of EggNOG-mapped background gene IDs (used when do_no_bg = FALSE for with_bg runs).
#' @param species character vector of g:Profiler organism codes (e.g., "hsapiens"). Unlimited.
#'                 Named vector allowed for pretty labels (names = codes, values = labels). If unnamed, labels=codes.
#' @param outdir Output directory root.
#' @param do_no_bg logical; if TRUE, also run no-background (genome-wide) analyses.
#' @param sources g:Profiler sources to query.
#' @param user_threshold numeric p-value threshold.
#' @param correction_method FDR method ("fdr","gSCS","bonferroni").
#' @param evcodes logical; return evidence codes.
#' @param significant Logical; if TRUE, g:Profiler filters to significant results.
#' @param sleep_sec seconds to pause between species to avoid throttling.
#' @param verbose logical; print progress.
#' @param go_obo Optional path to a GO OBO file; if NULL, EchoGO caches one per session.
#' @param significance_rule DE significance rule saved in run metadata.
#' @param resolver_definition Shared canonical resolver saved in run metadata.
#' @return list with per-run data.frames and a `$paths` list of written files.
#' @export
run_gprofiler_cross_species <- function(
    de_genes,
    bg_genes = NULL,
    species = c("hsapiens","mmusculus","rnorvegicus","ggallus","drerio","dmelanogaster","celegans"),
    outdir = "cross_species_gprofiler",
    do_no_bg = TRUE,
    sources = c("GO:BP","GO:MF","GO:CC","KEGG"),
    user_threshold = 0.05,
    correction_method = "fdr",
    evcodes = TRUE,
    significant = FALSE,
    sleep_sec = 0.5,
    verbose = TRUE,
    go_obo = NULL,
    significance_rule = NULL,
    resolver_definition = NULL
) {
  stopifnot(is.character(de_genes), length(de_genes) > 0)
  if (is.null(bg_genes) || !is.character(bg_genes) || !length(bg_genes)) {
    stop("A non-empty custom background is required for the background-aware run.", call. = FALSE)
  }
  de_genes <- unique(trimws(de_genes[!is.na(de_genes) & nzchar(trimws(de_genes))]))
  bg_genes <- unique(trimws(bg_genes[!is.na(bg_genes) & nzchar(trimws(bg_genes))]))
  outside <- setdiff(de_genes, bg_genes)
  if (length(outside)) {
    stop(
      "The resolved foreground must be a subset of the custom background. Offending values: ",
      paste(utils::head(outside, 8L), collapse = ", "),
      call. = FALSE
    )
  }
  if (identical(de_genes, bg_genes)) {
    stop("The foreground and custom background are the same vector; the run is degenerate.", call. = FALSE)
  }

  # Short-circuit if outputs already exist and user didn't force a rerun
  if (!isTRUE(getOption("EchoGO.force_rerun_gprofiler", FALSE))) {
    have_cached <- length(list.files(
      outdir, pattern = "^gprofiler_.*_(with_bg|nobg).*\\.csv$", recursive = TRUE
    )) > 0
    if (have_cached) {
      if (!file.exists(file.path(outdir, "run_manifest.json"))) {
        stop(
          "Cached g:Profiler CSVs were found without a v0.1.3 run manifest. ",
          "Historical response metadata cannot be reconstructed. Use a new output directory ",
          "or force a new run explicitly.",
          call. = FALSE
        )
      }
      if (isTRUE(verbose)) message("   - using cached g:Profiler CSVs in: ", outdir)
      return(invisible(list(paths = list(
        written = file.path(outdir, "run_manifest.json"), outdir = outdir
      ))))
    }
  }

  # Normalize & prepare output dirs (canonical substructure)
  outdir   <- normalizePath(outdir, winslash = "/", mustWork = FALSE)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  nobg_dir <- file.path(outdir, "no_background_genome_wide")
  bg_dir   <- file.path(outdir, "with_custom_background")
  dir.create(nobg_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(bg_dir,   showWarnings = FALSE, recursive = TRUE)

  # --- species labels: support named vector (names=codes, values=labels) or plain vector ---
  if (is.null(names(species))) {
    sp_codes  <- species
    sp_labels <- species
  } else {
    sp_codes  <- names(species)
    sp_labels <- as.character(unname(species))
  }

  # --- Prepare GO depth function (prefer GO.db if installed, else cached OBO) ---
  depth_fun <- NULL
  depth_backend <- NA_character_
  if (requireNamespace("GO.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
    if (verbose) message("Using GO.db for depth...")
    depth_backend <- "GO.db"
    depth_fun <- function(term_ids) {
      term_ids <- as.character(term_ids)
      out <- rep(NA_integer_, length(term_ids))
      valid <- !is.na(term_ids) & startsWith(term_ids, "GO:")
      keys <- unique(term_ids[valid])
      if (!length(keys)) return(out)

      safe_mget <- function(map) {
        tryCatch(
          AnnotationDbi::mget(keys, map, ifnotfound = NA),
          error = function(e) setNames(rep(list(NA), length(keys)), keys)
        )
      }
      anc_bp <- safe_mget(GO.db::GOBPANCESTOR)
      anc_mf <- safe_mget(GO.db::GOMFANCESTOR)
      anc_cc <- safe_mget(GO.db::GOCCANCESTOR)

      depth_by_key <- vapply(keys, function(key) {
        ancestors <- unique(c(
          unlist(anc_bp[[key]], use.names = FALSE),
          unlist(anc_mf[[key]], use.names = FALSE),
          unlist(anc_cc[[key]], use.names = FALSE)
        ))
        ancestors <- ancestors[!is.na(ancestors) & ancestors != "all"]
        as.integer(length(ancestors))
      }, integer(1))
      out[valid] <- unname(depth_by_key[term_ids[valid]])
      out
    }
  } else {
    if (is.null(go_obo)) {
      go_obo <- system.file("extdata","go-basic.obo", package = "EchoGO")
      if (go_obo == "") {
        go_obo <- file.path(tempdir(), "go-basic.obo")
        if (!file.exists(go_obo)) {
          if (verbose) message("Downloading GO OBO for depth (cached): ", go_obo)
          utils::download.file("http://purl.obolibrary.org/obo/go.obo",
                               destfile = go_obo, mode = "wb", quiet = !verbose)
        }
      }
    }
    if (verbose) message("Using ontologyIndex with cached OBO for depth...")
    depth_backend <- "ontologyIndex"
    go_ont <- ontologyIndex::get_ontology(go_obo, extract_tags = "minimal")
    depth_fun <- function(term_ids) {
      vapply(as.character(term_ids), function(term_id) {
        if (is.na(term_id) || !startsWith(term_id, "GO:")) return(NA_integer_)
        if (!term_id %in% go_ont$id) return(NA_integer_)
        as.integer(length(ontologyIndex::get_ancestors(go_ont, term_id)))
      }, integer(1))
    }
  }

  results_list <- list()
  paths <- list(written = character(0),
                with_bg_dir = bg_dir,
                nobg_dir    = nobg_dir,
                outdir      = outdir,
                depth_backend = depth_backend)
  summary_log <- data.frame(
    species_code = character(), species_label = character(),
    mode = character(), n_sig = integer(),
    stringsAsFactors = FALSE
  )
  manifest_runs <- list()
  run_timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  relative_to_outdir <- function(path) {
    if (is.null(path)) return(NULL)
    normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
    prefix <- paste0(outdir, "/")
    if (startsWith(normalized, prefix)) substring(normalized, nchar(prefix) + 1L) else normalized
  }
  scalar_or_values <- function(x) {
    x <- unique(stats::na.omit(x))
    if (!length(x)) return(NULL)
    if (length(x) == 1L) unname(x[[1]]) else unname(x)
  }
  record_run <- function(response, organism, label, mode, directory, background) {
    suffix <- if (identical(mode, "custom_background")) "with_bg" else "nobg"
    stem <- file.path(directory, paste0("gprofiler_", label, "_", suffix))
    query_file <- paste0(stem, "_query.txt")
    background_file <- if (identical(mode, "custom_background")) paste0(stem, "_background.txt") else NULL
    metadata_file <- paste0(stem, "_metadata.json")
    result_file <- paste0(stem, ".csv")
    if (!file.exists(result_file)) result_file <- NULL
    writeLines(de_genes, query_file)
    if (!is.null(background_file)) writeLines(background, background_file)
    result_table <- if (!is.null(response) && is.data.frame(response$result)) response$result else data.frame()
    effective_query <- if (nrow(result_table) && "query_size" %in% names(result_table)) {
      scalar_or_values(result_table$query_size)
    } else if (!is.null(response$meta$effective_query_size)) {
      response$meta$effective_query_size
    } else NULL
    effective_domain <- if (nrow(result_table) && "effective_domain_size" %in% names(result_table)) {
      scalar_or_values(result_table$effective_domain_size)
    } else if (!is.null(response$meta$effective_domain_size)) {
      response$meta$effective_domain_size
    } else NULL
    query_hash <- unname(tools::md5sum(query_file))
    background_hash <- if (!is.null(background_file)) unname(tools::md5sum(background_file)) else NULL
    metadata <- list(
      echogo_version = tryCatch(as.character(utils::packageVersion("EchoGO")), error = function(e) NA_character_),
      echogo_git_commit = .echogo_git_commit(),
      run_timestamp = run_timestamp,
      timezone = Sys.timezone(),
      r_version = R.version.string,
      gprofiler2_version = as.character(utils::packageVersion("gprofiler2")),
      execution = "live",
      organism_code = organism,
      organism_label = label,
      requested_sources = sources,
      correction_method = correction_method,
      user_threshold = user_threshold,
      significant = significant,
      ordered_query = FALSE,
      multi_query = FALSE,
      background_mode = mode,
      submitted_foreground_count = length(de_genes),
      submitted_background_count = if (is.null(background)) NULL else length(background),
      effective_query_size = effective_query,
      effective_domain_size = effective_domain,
      de_significance_rule = significance_rule,
      resolver_definition = resolver_definition,
      query_hash = query_hash,
      background_hash = background_hash,
      response_metadata = if (!is.null(response)) response$meta else NULL,
      status = if (is.null(response)) "request_failed" else if (nrow(result_table)) "completed" else "no_results"
    )
    jsonlite::write_json(metadata, metadata_file, auto_unbox = TRUE, pretty = TRUE, na = "null")
    entry <- c(
      metadata[c(
        "organism_code", "organism_label", "background_mode",
        "submitted_foreground_count", "submitted_background_count",
        "effective_query_size", "effective_domain_size", "query_hash",
        "background_hash", "status"
      )],
      list(
        result_file = relative_to_outdir(result_file),
        metadata_file = relative_to_outdir(metadata_file),
        query_file = relative_to_outdir(query_file),
        background_file = relative_to_outdir(background_file)
      )
    )
    list(entry = entry, files = c(query_file, background_file, metadata_file))
  }

  # helper to run a single gost call safely
  .run_gost <- function(query, organism, custom_bg, mode_label) {
    if (verbose) message("  - ", organism, " [", mode_label, "]")
    request_failed <- FALSE
    res <- tryCatch(
      .echogo_gost(
        query = query,
        organism = organism,
        custom_bg = custom_bg,
        sources = sources,
        correction_method = correction_method,
        user_threshold = user_threshold,
        evcodes = evcodes,
        significant = significant,
        ordered_query = FALSE,
        multi_query = FALSE
      ),
      error = function(e) {
        request_failed <<- TRUE
        warning("gost failed for ", organism, " (", mode_label, "): ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(res) && !request_failed) {
      # gprofiler2 returns NULL for a successful request with no reportable
      # terms. Preserve that outcome as an explicit empty response so metadata
      # does not conflate it with an exception or transport failure.
      res <- list(
        result = data.frame(),
        meta = list(service_outcome = "no_results_returned")
      )
    }
    res
  }

  for (i in seq_along(sp_codes)) {
    sp <- sp_codes[i]
    lab <- sp_labels[i]
    if (verbose) message(">> g:Profiler for ", sp, " (", lab, ")")

    # --- WITH custom background ---
    res_bg <- .run_gost(de_genes, sp, custom_bg = bg_genes, mode_label = "with_bg")
    out_subdir <- bg_dir
    if (!is.null(res_bg) && nrow(res_bg$result) > 0) {
      tbl <- res_bg$result
      tbl$fold_enrichment <- (tbl$intersection_size / tbl$query_size) / (tbl$term_size / tbl$effective_domain_size)
      tbl$depth <- depth_fun(tbl$term_id)
      tbl$species_code  <- sp
      tbl$species_label <- lab
      csv  <- file.path(out_subdir, paste0("gprofiler_", lab, "_with_bg.csv"))
      xlsx <- file.path(out_subdir, paste0("gprofiler_", lab, "_with_bg.xlsx"))
      readr::write_csv(tbl, csv)
      openxlsx::write.xlsx(tbl, xlsx, asTable = TRUE, overwrite = TRUE)
      results_list[[paste0(lab, "_with_bg")]] <- tbl
      paths$written <- c(paths$written, csv, xlsx)
      summary_log <- rbind(summary_log, data.frame(
        species_code = sp, species_label = lab, mode = "with_bg",
        n_sig = sum(!is.na(tbl$p_value) & tbl$p_value <= user_threshold)
      ))
    } else {
      file.create(file.path(out_subdir, paste0("gprofiler_", lab, "_with_bg_NO_RESULTS.txt")))
      results_list[[paste0(lab, "_with_bg")]] <- NULL
      summary_log <- rbind(summary_log, data.frame(
        species_code = sp, species_label = lab, mode = "with_bg", n_sig = 0
      ))
    }
    recorded <- record_run(
      res_bg, sp, lab, "custom_background", bg_dir, bg_genes
    )
    manifest_runs[[length(manifest_runs) + 1L]] <- recorded$entry
    paths$written <- c(paths$written, recorded$files)

    # --- NO background (genome-wide) ---
    if (isTRUE(do_no_bg)) {
      res_nb <- .run_gost(de_genes, sp, custom_bg = NULL, mode_label = "nobg")
      out_subdir <- nobg_dir
      if (!is.null(res_nb) && nrow(res_nb$result) > 0) {
        tbl <- res_nb$result
        tbl$fold_enrichment <- (tbl$intersection_size / tbl$query_size) / (tbl$term_size / tbl$effective_domain_size)
        tbl$depth <- depth_fun(tbl$term_id)
        tbl$species_code  <- sp
        tbl$species_label <- lab
        csv  <- file.path(out_subdir, paste0("gprofiler_", lab, "_nobg.csv"))
        xlsx <- file.path(out_subdir, paste0("gprofiler_", lab, "_nobg.xlsx"))
        readr::write_csv(tbl, csv)
        openxlsx::write.xlsx(tbl, xlsx, asTable = TRUE, overwrite = TRUE)
        results_list[[paste0(lab, "_nobg")]] <- tbl
        paths$written <- c(paths$written, csv, xlsx)
        summary_log <- rbind(summary_log, data.frame(
          species_code = sp, species_label = lab, mode = "nobg",
          n_sig = sum(!is.na(tbl$p_value) & tbl$p_value <= user_threshold)
        ))
      } else {
        file.create(file.path(out_subdir, paste0("gprofiler_", lab, "_nobg_NO_RESULTS.txt")))
        results_list[[paste0(lab, "_nobg")]] <- NULL
        summary_log <- rbind(summary_log, data.frame(
          species_code = sp, species_label = lab, mode = "nobg", n_sig = 0
        ))
      }
      recorded <- record_run(
        res_nb, sp, lab, "no_background_genome_wide", nobg_dir, NULL
      )
      manifest_runs[[length(manifest_runs) + 1L]] <- recorded$entry
      paths$written <- c(paths$written, recorded$files)
    }

    if (sleep_sec > 0) Sys.sleep(sleep_sec)
  }

  summary_file <- file.path(outdir, "gprofiler_enrichment_summary.csv")
  readr::write_csv(summary_log, summary_file)
  manifest_file <- file.path(outdir, "run_manifest.json")
  run_manifest <- list(
    schema_version = "1.0",
    echogo_version = tryCatch(as.character(utils::packageVersion("EchoGO")), error = function(e) NA_character_),
    echogo_git_commit = .echogo_git_commit(),
    generated = run_timestamp,
    timezone = Sys.timezone(),
    execution = "live",
    vector_contract = "shared_portable_canonical_organism_context",
    explicit_species_specific_ortholog_mapping = FALSE,
    submitted_foreground_count = length(de_genes),
    submitted_background_count = length(bg_genes),
    significance_rule = significance_rule,
    resolver_definition = resolver_definition,
    historical_cached_output_limitation = paste(
      "Pre-v0.1.3 cached outputs did not preserve response metadata;",
      "their exact historical g:Profiler data release cannot be recovered."
    ),
    runs = manifest_runs
  )
  jsonlite::write_json(run_manifest, manifest_file, auto_unbox = TRUE, pretty = TRUE, na = "null")
  paths$written <- unique(c(paths$written, summary_file, manifest_file))
  results_list$summary <- summary_log
  results_list$paths <- paths
  return(results_list)
}

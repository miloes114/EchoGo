#' Run g:Profiler enrichment across multiple species (customizable, cached GO depth)
#'
#' @param de_genes character vector of EggNOG-mapped DE gene IDs.
#' @param bg_genes character vector of EggNOG-mapped background gene IDs (used when do_no_bg = FALSE for with_bg runs).
#' @param species character vector of g:Profiler organism codes (e.g., "hsapiens"). Unlimited.
#'                 Named vector allowed for pretty labels (names = codes, values = labels). If unnamed, labels=codes.
#' @param outdir output directory root (keeps your structure).
#' @param do_no_bg logical; if TRUE, also run no-background (genome-wide) analyses.
#' @param sources g:Profiler sources to query.
#' @param user_threshold numeric p-value threshold.
#' @param correction_method FDR method ("fdr","gSCS","bonferroni").
#' @param evcodes logical; return evidence codes.
#' @param significant logical; if TRUE g:Profiler filters to significant; we keep FALSE then filter ourselves as needed.
#' @param sleep_sec seconds to pause between species to avoid throttling.
#' @param verbose logical; print progress.
#' @param go_obo optional path to a GO OBO file; if NULL we auto-cache one per session.
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
    go_obo = NULL
) {
  stopifnot(is.character(de_genes), length(de_genes) > 0)
  if (!is.null(bg_genes)) stopifnot(is.character(bg_genes))

  # Short-circuit if outputs already exist and user didn't force a rerun
  if (!isTRUE(getOption("EchoGO.force_rerun_gprofiler", FALSE))) {
    have_cached <- length(list.files(
      outdir, pattern = "^gprofiler_.*_(with_bg|nobg).*\\.csv$", recursive = TRUE
    )) > 0
    if (have_cached) {
      if (isTRUE(verbose)) message("   · using cached g:Profiler CSVs in: ", outdir)
      return(invisible(NULL))  # results are on disk; upstream will load them
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
    if (verbose) message("Using GO.db for depth…")
    depth_backend <- "GO.db"
    depth_fun <- function(term_id) {
      if (!startsWith(term_id, "GO:")) return(NA_integer_)
      get_anc <- function(map) {
        x <- tryCatch(AnnotationDbi::mget(term_id, map, ifnotfound = NA)[[1]],
                      error = function(e) NA)
        x <- unlist(x, use.names = FALSE)
        x[!is.na(x) & x != "all"]
      }
      anc <- unique(c(
        get_anc(GO.db::GOBPANCESTOR),
        get_anc(GO.db::GOMFANCESTOR),
        get_anc(GO.db::GOCCANCESTOR)
      ))
      if (!length(anc)) 0L else length(anc)
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
    if (verbose) message("Using ontologyIndex with cached OBO for depth…")
    depth_backend <- "ontologyIndex"
    go_ont <- ontologyIndex::get_ontology(go_obo, extract_tags = "minimal")
    depth_fun <- function(term_id) {
      if (!startsWith(term_id, "GO:")) return(NA_integer_)
      if (!term_id %in% go_ont$id) return(NA_integer_)
      length(ontologyIndex::get_ancestors(go_ont, term_id))
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

  # helper to run a single gost call safely
  .run_gost <- function(query, organism, custom_bg, mode_label) {
    if (verbose) message("  - ", organism, " [", mode_label, "]")
    res <- tryCatch(
      gprofiler2::gost(
        query = query,
        organism = organism,
        custom_bg = custom_bg,
        sources = sources,
        correction_method = correction_method,
        user_threshold = user_threshold,
        evcodes = evcodes,
        significant = significant
      ),
      error = function(e) {
        warning("gost failed for ", organism, " (", mode_label, "): ", conditionMessage(e))
        NULL
      }
    )
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
      tbl$depth <- vapply(tbl$term_id, depth_fun, integer(1))
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

    # --- NO background (genome-wide) ---
    if (isTRUE(do_no_bg)) {
      res_nb <- .run_gost(de_genes, sp, custom_bg = NULL, mode_label = "nobg")
      out_subdir <- nobg_dir
      if (!is.null(res_nb) && nrow(res_nb$result) > 0) {
        tbl <- res_nb$result
        tbl$fold_enrichment <- (tbl$intersection_size / tbl$query_size) / (tbl$term_size / tbl$effective_domain_size)
        tbl$depth <- vapply(tbl$term_id, depth_fun, integer(1))
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
    }

    if (sleep_sec > 0) Sys.sleep(sleep_sec)
  }

  readr::write_csv(summary_log, file.path(outdir, "gprofiler_enrichment_summary.csv"))
  results_list$summary <- summary_log
  results_list$paths <- paths
  return(results_list)
}

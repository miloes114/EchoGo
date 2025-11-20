#' Load g:Profiler enrichment results
#'
#' Loads per-species enrichment tables from a directory, filtering for FDR < 0.05,
#' and adds species-specific logical columns (`in_<species>` or `in_<species>_nobg`).
#'
#' @param path Directory containing enrichment CSVs.
#' @param mode Either "with_bg" or "nobg" to define file naming logic.
#' @param species_map Named character vector mapping g:Profiler species codes -> labels.
#'   The filenames are assumed to use the *labels* (values of species_map).
#' @return A data frame of significant enrichment hits with species flags.
#' @keywords internal
load_gprofiler_results <- function(path, mode, species_map) {
  result_list <- list()
  sp_labels <- unname(species_map)

  for (sp_label in sp_labels) {
    # Accept both filename conventions:
    # - gprofiler_<label>_<mode>.csv                (current writer)
    # - gprofiler_<label>_<mode>_enrichment.csv     (legacy)
    candidates <- c(
      file.path(path, paste0("gprofiler_", sp_label, "_", mode, ".csv")),
      file.path(path, paste0("gprofiler_", sp_label, "_", mode, "_enrichment.csv"))
    )
    file_path <- candidates[file.exists(candidates)][1]
    if (!length(file_path) || is.na(file_path)) next

    df <- readr::read_csv(file_path, show_col_types = FALSE)
    if (!"p_value" %in% names(df)) stop("g:Profiler CSV missing 'p_value': ", basename(file_path))
    df <- dplyr::filter(df, .data$p_value <= 0.05)

    if (nrow(df) > 0) {
      flag_col <- paste0("in_", sp_label, if (mode == "nobg") "_nobg" else "")
      df[[flag_col]] <- TRUE
      result_list[[sp_label]] <- df
    }
  }

  if (!length(result_list)) return(tibble::tibble())
  dplyr::bind_rows(result_list)
}


#' @name summarize_gprofiler_by_term
#' @title Summarize g:Profiler enrichment results by term
#' @description Aggregates species-specific enrichment hits into a single row per GO/KEGG term.
#' @param gprof_df g:Profiler result data frame (one or more species).
#' @param species_labels Character vector of species *labels* (e.g., c("hsapiens","mmusculus") or display names).
#' @param mode Either "bg" or "nobg", for with- or no-background enrichment.
#' @return A summary data frame with fold enrichment, p-value, species flags, and merged gene sets.
#' @keywords internal
#' @importFrom rlang sym
summarize_gprofiler_by_term <- function(gprof_df, species_labels, mode = "bg") {
  if (nrow(gprof_df) == 0) return(NULL)

  # standardize column names from writer
  if (!all(c("term_id","term_name","source","fold_enrichment","p_value","intersection") %in% names(gprof_df))) {
    stop("g:Profiler frame missing required columns (need term_id, term_name, source, fold_enrichment, p_value, intersection).")
  }

  # target column names expected by downstream scoring
  avg_fold_col <- if (mode == "bg") "avg_fold_gprof_bg" else "avg_fold_gprof_nobg"
  min_pval_col <- if (mode == "bg") "min_pval_gprof_bg" else "min_pval_gprof_nobg"
  gene_col_out <- if (mode == "bg") "genes_gprofiler_bg" else "genes_gprofiler_nobg"

  # keep original numeric columns but also create deterministic aliases for summarise()
  gprof_df <- gprof_df %>%
    dplyr::rename(
      ontology_gprof = source
    )

  # ensure per-species flags exist for all requested labels
  for (sp in species_labels) {
    col <- paste0("in_", sp, if (mode == "nobg") "_nobg" else "")
    if (!col %in% names(gprof_df)) gprof_df[[col]] <- FALSE
  }

  gprof_df %>%
    dplyr::group_by(term_id, term_name, ontology_gprof) %>%
    dplyr::summarise(
      !!avg_fold_col := mean(.data$fold_enrichment, na.rm = TRUE),
      !!min_pval_col := min(.data$p_value, na.rm = TRUE),
      dplyr::across(dplyr::starts_with("in_"), ~ any(.x, na.rm = TRUE)),
      !!gene_col_out := {
        # merge all intersections across species; keep unique, non-empty
        vals <- paste(.data$intersection, collapse = ",")
        paste(unique(stats::na.omit(trimws(strsplit(vals, ",")[[1]]))), collapse = ", ")
      },
      .groups = "drop"
    )
}


#' Merge enrichment results from all sources
#'
#' Combines GOseq, g:Profiler with background, and g:Profiler without background
#' into a unified consensus table before scoring.
#'
#' @param goseq_df GOseq data frame.
#' @param gprof_bg_df Summarized g:Profiler (with background).
#' @param gprof_nobg_df Summarized g:Profiler (no background).
#' @return A merged data frame containing all enrichment terms.
#' @keywords internal
merge_enrichment_sources <- function(goseq_df, gprof_bg_df, gprof_nobg_df) {
  dplyr::full_join(gprof_bg_df, gprof_nobg_df, by = "term_id", suffix = c(".bg", ".nobg")) %>%
    dplyr::full_join(goseq_df, by = "term_id") %>%
    dplyr::mutate(
      term_name   = dplyr::coalesce(.data$term_name, .data$term_name.bg, .data$term_name.nobg),
      all_genes   = dplyr::coalesce(.data$gene_names, .data$genes_gprofiler_bg, .data$genes_gprofiler_nobg),
      ontology_gprof = dplyr::coalesce(.data$ontology_gprof.bg, .data$ontology_gprof.nobg),
      ontology    = dplyr::coalesce(.data$ontology, .data$ontology_gprof)
    ) %>%
    dplyr::filter(!is.na(.data$term_name))
}


#' Score and classify consensus enrichment terms
#'
#' Applies metrics such as species support, consensus scores, robustness, origin tags, etc.
#'
#' @param df Data frame output from `merge_enrichment_sources()`.
#' @return Scored consensus enrichment table.
#' @keywords internal
score_consensus_terms <- function(df,
                                  P_CAP = 6,
                                  FOLD_CAP = 8,
                                  W_STRICT = list(
                                    w_goseq   = 1.0,
                                    w_bgprev  = 0.8,
                                    w_foldbg  = 0.40,
                                    w_foldgs  = 0.40,
                                    w_p       = 0.40
                                  ),
                                  W_ALL = list(
                                    w_goseq   = 1.0,
                                    w_bgprev  = 0.7,
                                    w_nbprev  = 0.3,
                                    w_foldbg  = 0.35,
                                    w_foldnb  = 0.25,
                                    w_foldgs  = 0.35,
                                    w_p       = 0.35
                                  )) {

  # Identify species flags robustly (exclude the GOseq indicator)
  all_in_cols <- grep("^in_", names(df), value = TRUE)
  nobg_cols   <- grep("^in_.*_nobg$", all_in_cols, value = TRUE)
  bg_cols     <- setdiff(setdiff(all_in_cols, nobg_cols), "in_goseq")

  n_bg_cols   <- length(bg_cols)
  n_nobg_cols <- length(nobg_cols)

  # helpers
  cap_logp <- function(p, P_CAP) {
    p <- dplyr::coalesce(p, 1)
    p <- pmax(p, .Machine$double.xmin)
    pmin(-log10(p), P_CAP) / P_CAP
  }
  log2p1_cap <- function(x, cap) log2(1 + pmin(dplyr::coalesce(x, 0), cap))

  # NEW: shallow terms get less weight; NA depth -> 1 (no penalty); KEGG stays NA so weight=1
  depth_w <- function(d, maxd = 12L) {
    d <- suppressWarnings(as.integer(d))
    ifelse(is.na(d), 1, pmin(d / maxd, 1))
  }

  df %>%
    dplyr::mutate(
      # species counts
      num_species_gprof_bg   = if (n_bg_cols)   rowSums(dplyr::across(dplyr::all_of(bg_cols)),   na.rm = TRUE) else 0,
      num_species_gprof_nobg = if (n_nobg_cols) rowSums(dplyr::across(dplyr::all_of(nobg_cols)), na.rm = TRUE) else 0,

      # normalized prevalence (0..1); treat 0 cols as 0
      bg_prev   = if (n_bg_cols)   num_species_gprof_bg   / n_bg_cols else 0,
      nobg_prev = if (n_nobg_cols) num_species_gprof_nobg / n_nobg_cols else 0,

      # robustness flags
      is_robust_nobg = !is.na(min_pval_gprof_nobg) & min_pval_gprof_nobg <= 0.05 & num_species_gprof_nobg >= 2,

      # source count (use isTRUE() to ensure only TRUE counts, not NA/FALSE)
      sources_count = as.integer(isTRUE(.data$in_goseq)) +
        dplyr::coalesce(rowSums(dplyr::across(dplyr::all_of(bg_cols)),   na.rm = TRUE), 0L) +
        dplyr::coalesce(rowSums(dplyr::across(dplyr::all_of(nobg_cols)), na.rm = TRUE), 0L),

      # capped/normalized components
      comp_p_strict = cap_logp(min_pval_gprof_bg, P_CAP),
      comp_p_all    = pmax(
        cap_logp(min_pval_gprof_bg,   P_CAP),
        cap_logp(min_pval_gprof_nobg, P_CAP),
        cap_logp(min_pval_goseq,      P_CAP)
      ),

      # CHANGED: apply depth weight to fold components
      comp_fold_bg   = depth_w(depth) * log2p1_cap(avg_fold_gprof_bg,   FOLD_CAP),
      comp_fold_nb   = depth_w(depth) * log2p1_cap(avg_fold_gprof_nobg, FOLD_CAP),
      comp_fold_gs   = depth_w(depth) * log2p1_cap(fold_enrichment_goseq, FOLD_CAP),

      # base indicator (count only explicit TRUE)
      comp_goseq = as.integer(isTRUE(in_goseq)),

      # final scores
      consensus_score =
        W_STRICT$w_goseq  * comp_goseq +
        W_STRICT$w_bgprev * bg_prev +
        W_STRICT$w_foldbg * comp_fold_bg +
        W_STRICT$w_foldgs * comp_fold_gs +
        W_STRICT$w_p      * comp_p_strict,

      consensus_score_all =
        W_ALL$w_goseq   * comp_goseq +
        W_ALL$w_bgprev  * bg_prev +
        W_ALL$w_nbprev  * nobg_prev +
        W_ALL$w_foldbg  * comp_fold_bg +
        W_ALL$w_foldnb  * comp_fold_nb +
        W_ALL$w_foldgs  * comp_fold_gs +
        W_ALL$w_p       * comp_p_all,

      origin = dplyr::case_when(
        !is.na(in_goseq) & !is.na(min_pval_gprof_bg)   & min_pval_gprof_bg   <= 0.05 ~ "GO terms - Consensus (with BG)",
        !is.na(in_goseq) & !is.na(min_pval_gprof_nobg) & min_pval_gprof_nobg <= 0.05 ~ "GO terms - Consensus (no BG)",
        !is.na(in_goseq) ~ "GO terms - GOseq only",
        is.na(in_goseq)  & !is.na(min_pval_gprof_bg)   & min_pval_gprof_bg   <= 0.05 ~ "GO terms - g:Profiler only (with BG)",
        is.na(in_goseq)  & !is.na(min_pval_gprof_nobg) & min_pval_gprof_nobg <= 0.05 ~ "GO terms - g:Profiler only (no BG)",
        TRUE ~ "Other"
      ),

      significant_in_any = dplyr::if_else(
        (!is.na(min_pval_goseq)       & min_pval_goseq       <= 0.05) |
          (!is.na(min_pval_gprof_bg)    & min_pval_gprof_bg    <= 0.05) |
          (!is.na(min_pval_gprof_nobg)  & min_pval_gprof_nobg  <= 0.05),
        TRUE, FALSE
      ),

      source_origin = purrr::pmap_chr(
        list(in_goseq, min_pval_gprof_bg, min_pval_gprof_nobg),
        function(in_goseq, bg, nobg) {
          sources <- character(0)
          if (!is.na(in_goseq)) sources <- c(sources, "GOseq")
          if (!is.na(bg)   && bg   <= 0.05) sources <- c(sources, "g:Profiler_BG")
          if (!is.na(nobg) && nobg <= 0.05) sources <- c(sources, "g:Profiler_noBG")
          paste(sources, collapse = "+")
        }
      )
    ) %>%
    dplyr::select(-dplyr::starts_with("."))
}

audit_consensus_flags <- function(df, bg_cols = NULL, nobg_cols = NULL) {
  if (is.null(bg_cols) || is.null(nobg_cols)) {
    all_in_cols <- grep("^in_", names(df), value = TRUE)
    nobg_cols   <- grep("^in_.*_nobg$", all_in_cols, value = TRUE)
    bg_cols     <- setdiff(setdiff(all_in_cols, nobg_cols), "in_goseq")
  }

  # Compute sums without dplyr::across() (safe even if 0 columns)
  bg_sum <- if (length(bg_cols))   rowSums(as.data.frame(df[, bg_cols,   drop = FALSE]), na.rm = TRUE) else rep(0L, nrow(df))
  nb_sum <- if (length(nobg_cols)) rowSums(as.data.frame(df[, nobg_cols, drop = FALSE]), na.rm = TRUE) else rep(0L, nrow(df))
  src_calc <- as.integer(isTRUE(df$in_goseq)) + bg_sum + nb_sum

  dplyr::bind_cols(
    df,
    tibble::tibble(.bg_sum = bg_sum, .nobg_sum = nb_sum, .src_calc = src_calc)
  ) %>%
    dplyr::filter(
      (origin == "GO terms - GOseq only" & (.bg_sum > 0 | .nobg_sum > 0)) |
        (sources_count != .src_calc)
    ) %>%
    dplyr::select(term_id, term_name, origin, sources_count, .src_calc, .bg_sum, .nobg_sum)
}

# ---- Utility: compute GO depth for a data.frame with term_id ----
.echogo_add_go_depth <- function(df, id_col = "term_id") {
  if (!id_col %in% names(df)) return(df)

  n <- nrow(df)

  # Robust per-row string extraction (handles list-columns, factors, etc.)
  extract_first_string <- function(x) {
    vapply(seq_len(n), function(i) {
      xi <- df[[id_col]][[i]]
      if (length(xi) == 0) "" else as.character(xi)[1]
    }, character(1), USE.NAMES = FALSE)
  }

  # Keep GO-only (others like KEGG/REAC remain NA)
  normalize_goid <- function(vec_chr) {
    m <- regmatches(vec_chr, regexpr("GO:\\d{7}", vec_chr))
    out <- ifelse(nzchar(m), m, NA_character_)
    out
  }

  term_chr <- extract_first_string(df[[id_col]])
  goids    <- normalize_goid(term_chr)

  depth_vec <- rep(NA_integer_, n)
  used <- "none"

  if (requireNamespace("GO.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {

    valid <- !is.na(goids)
    keys  <- unique(goids[valid])

    safe_mget <- function(keys, map) {
      if (!length(keys)) return(setNames(vector("list", 0), character(0)))
      AnnotationDbi::mget(keys, map, ifnotfound = NA)
    }

    anc_bp <- safe_mget(keys, GO.db::GOBPANCESTOR)
    anc_mf <- safe_mget(keys, GO.db::GOMFANCESTOR)
    anc_cc <- safe_mget(keys, GO.db::GOCCANCESTOR)

    depth_vec[valid] <- vapply(goids[valid], function(gi) {
      ai_raw <- c(
        unlist(anc_bp[[gi]], use.names = FALSE),
        unlist(anc_mf[[gi]], use.names = FALSE),
        unlist(anc_cc[[gi]], use.names = FALSE)
      )
      ai <- ai_raw[!is.na(ai_raw) & ai_raw != "all"]
      if (!length(ai)) 0L else length(unique(ai))
    }, integer(1))

    used <- "GO.db"

  } else {
    # Fallback if you ship go-basic.obo with the package (inst/extdata/go-basic.obo)
    go_obo <- system.file("extdata","go-basic.obo", package = "EchoGO")
    if (requireNamespace("ontologyIndex", quietly = TRUE) && file.exists(go_obo)) {
      go_ont <- ontologyIndex::get_ontology(go_obo, extract_tags = "minimal")
      depth_vec <- vapply(goids, function(term_id) {
        if (!is.na(term_id) && term_id %in% go_ont$id) {
          length(ontologyIndex::get_ancestors(go_ont, term_id))
        } else NA_integer_
      }, integer(1))
      used <- "ontologyIndex"
    }
  }

  # --- Safe assignment (handles existing depth + tibble rules) ---
  depth_vec <- as.integer(depth_vec)

  df <- dplyr::ungroup(df)

  if ("depth" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(.depth_calc = depth_vec) %>%
      dplyr::mutate(depth = dplyr::coalesce(depth, .depth_calc)) %>%
      dplyr::select(-.depth_calc)
  } else {
    df <- df %>% dplyr::mutate(depth = depth_vec)
  }

  message(switch(used,
                 "GO.db"         = "Using GO.db for GO depth.",
                 "ontologyIndex" = "Using bundled go-basic.obo for GO depth (ontologyIndex).",
                 "none"          = "No GO depth backend available; depth set to NA."
  ))

  return(df)
}


#' @name build_consensus_table
#' @title Build Consensus Enrichment Table
#' @description Merges GOseq and g:Profiler enrichment results (with and without background)
#' into a scored consensus table with logical annotations, supporting downstream visualizations.
#' @param goseq_file Path to GOseq annotated enrichment CSV file.
#' @param gprofiler_dir Directory with g:Profiler enrichment subfolders (default: "cross_species_gprofiler").
#' @param species_map Named character vector mapping g:Profiler species codes -> labels.
#' @param output_dir Output folder to write final consensus Excel tables.
#' @return A scored consensus data frame.
#' @export
build_consensus_table <- function(
    goseq_file,
    gprofiler_dir = "cross_species_gprofiler",
    species_map = c("hsapiens" = "human","mmusculus" = "mouse","drerio" = "zebrafish"),
    output_dir = "consensus_enrichment"
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # GOseq input (already annotated)
  goseq_df <- readr::read_csv(goseq_file, show_col_types = FALSE) %>%
    dplyr::mutate(term_id = .data$clean_go_term, in_goseq = TRUE) %>%
    dplyr::select(term_id, term, ontology, foldEnrichment, over_represented_FDR, gene_names, in_goseq, depth) %>%
    dplyr::rename(
      fold_enrichment_goseq = foldEnrichment,
      min_pval_goseq        = over_represented_FDR,
      term_name             = term
    )

  # g:Profiler inputs (respect custom species set end-to-end)
  gprof_all_bg   <- load_gprofiler_results(file.path(gprofiler_dir, "with_custom_background"), "with_bg", species_map)
  gprof_all_nobg <- load_gprofiler_results(file.path(gprofiler_dir, "no_background_genome_wide"), "nobg", species_map)

  message("✅ g:Profiler with_bg terms: ", if (nrow(gprof_all_bg)) nrow(gprof_all_bg) else 0L)
  message("✅ g:Profiler no_bg terms: ",   if (nrow(gprof_all_nobg)) nrow(gprof_all_nobg) else 0L)

  # summaries (pass *labels* only)
  sp_labels <- unname(species_map)
  gprof_summary_bg   <- summarize_gprofiler_by_term(gprof_all_bg,   sp_labels, mode = "bg")
  gprof_summary_nobg <- summarize_gprofiler_by_term(gprof_all_nobg, sp_labels, mode = "nobg")

  # Fallback: ensure minimal structure exists with all required 'in_' columns
  if (is.null(gprof_summary_bg)) {
    gprof_summary_bg <- tibble::tibble(
      term_id           = character(),
      term_name         = character(),
      ontology_gprof    = character(),
      avg_fold_gprof_bg = numeric(),
      min_pval_gprof_bg = numeric(),
      genes_gprofiler_bg= character()
    )
    for (sp in sp_labels) gprof_summary_bg[[paste0("in_", sp)]] <- logical()
  }
  if (is.null(gprof_summary_nobg)) {
    gprof_summary_nobg <- tibble::tibble(
      term_id               = character(),
      term_name             = character(),
      ontology_gprof        = character(),
      avg_fold_gprof_nobg   = numeric(),
      min_pval_gprof_nobg   = numeric(),
      genes_gprofiler_nobg  = character()
    )
    for (sp in sp_labels) gprof_summary_nobg[[paste0("in_", sp, "_nobg")]] <- logical()
  }

  raw_consensus <- merge_enrichment_sources(goseq_df, gprof_summary_bg, gprof_summary_nobg)

  # NEW: add GO depth for all terms (GO only; others remain NA)
  raw_consensus <- .echogo_add_go_depth(raw_consensus, id_col = "term_id")

  # Ensure required columns exist with correct types (avoid mutate errors when empty)
  num_cols <- c("avg_fold_gprof_bg","min_pval_gprof_bg",
                "avg_fold_gprof_nobg","min_pval_gprof_nobg",
                "fold_enrichment_goseq","min_pval_goseq")
  chr_cols <- c("genes_gprofiler_bg","genes_gprofiler_nobg","gene_names")

  for (nm in num_cols) {
    if (!nm %in% names(raw_consensus)) raw_consensus[[nm]] <- NA_real_
    raw_consensus[[nm]] <- suppressWarnings(as.numeric(raw_consensus[[nm]]))
  }
  for (nm in chr_cols) {
    if (!nm %in% names(raw_consensus)) raw_consensus[[nm]] <- NA_character_
    raw_consensus[[nm]] <- as.character(raw_consensus[[nm]])
  }
  if (!"in_goseq" %in% names(raw_consensus)) raw_consensus[["in_goseq"]] <- NA
  raw_consensus[["in_goseq"]] <- as.logical(raw_consensus[["in_goseq"]])


  # Ensure columns for scoring exist (prevents mutate crashes)
  required_columns <- c(
    "avg_fold_gprof_bg", "min_pval_gprof_bg", "genes_gprofiler_bg",
    "avg_fold_gprof_nobg", "min_pval_gprof_nobg", "genes_gprofiler_nobg",
    "fold_enrichment_goseq", "min_pval_goseq", "in_goseq", "gene_names"
  )
  for (col in required_columns) if (!col %in% names(raw_consensus)) raw_consensus[[col]] <- NA

  consensus_df <- score_consensus_terms(raw_consensus) %>%
    dplyr::mutate(
      ontology = dplyr::case_when(
        .data$ontology %in% c("GO:BP", "BP") ~ "BP",
        .data$ontology %in% c("GO:MF", "MF") ~ "MF",
        .data$ontology %in% c("GO:CC", "CC") ~ "CC",
        TRUE ~ .data$ontology
      )
    ) %>%
    dplyr::select(-dplyr::any_of(c("term", "ontology_gprof.bg", "ontology_gprof.nobg", "ontology_gprof")))

  # Optional sanity check
  bad_rows <- audit_consensus_flags(consensus_df)
  if (nrow(bad_rows)) {
    message("⚠️ Inconsistencies detected in sources/species flags:")
    print(bad_rows, n = 20)
  }

  # Save full table
  openxlsx::write.xlsx(
    consensus_df,
    file = file.path(output_dir, "consensus_enrichment_results_with_and_without_bg.xlsx"),
    asTable = TRUE
  )

  # Save filtered exploration table
  exploration_cols <- c(
    "term_id", "term_name", "ontology", "depth",
    "fold_enrichment_goseq", "min_pval_goseq",
    "avg_fold_gprof_bg", "min_pval_gprof_bg",
    "avg_fold_gprof_nobg", "min_pval_gprof_nobg",
    grep("^in_", names(consensus_df), value = TRUE),
    "all_genes", "consensus_score", "consensus_score_all",
    "sources_count", "num_species_gprof_bg", "num_species_gprof_nobg",
    "is_robust_nobg", "significant_in_any", "origin"
  )
  openxlsx::write.xlsx(
    consensus_df[, exploration_cols, drop = FALSE] %>%
      dplyr::mutate(dplyr::across(dplyr::starts_with("in_"), as.logical)),
    file = file.path(output_dir, "consensus_enrichment_exploration_clean.xlsx"),
    asTable = TRUE
  )

  message("✅ Consensus table saved to: ", output_dir)
  consensus_df
}

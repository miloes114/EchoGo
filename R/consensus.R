# Scoreless term-evidence assembly -------------------------------------------

.echogo_validate_support_columns <- function(df, columns) {
  columns <- intersect(columns, names(df))
  invalid <- columns[!vapply(df[columns], is.logical, logical(1))]
  if (length(invalid)) {
    stop(
      "All logical support columns must contain only TRUE, FALSE, or NA. Invalid column(s): ",
      paste(invalid, collapse = ", "), call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Deprecated score adapter
#'
#' `consensus_score` and `consensus_score_all` were retired in EchoGO v0.1.4.
#' This adapter provides evidence profiles for legacy wide tables but deliberately
#' does not recreate either rejected composite score.
#'
#' @param df Legacy wide consensus table.
#' @param ... Ignored retired score parameters.
#' @return A legacy-compatible scoreless evidence view without score columns.
#' @export
score_consensus_terms <- function(df, ...) {
  warning(
    "score_consensus_terms() is deprecated in v0.1.4: canonical composite scores ",
    "are retired. Use build_consensus_table() or assemble_echogo_evidence().",
    call. = FALSE
  )
  all_in <- grep("^in_", names(df), value = TRUE)
  nobg_in <- grep("^in_.*_nobg$", all_in, value = TRUE)
  bg_in <- setdiff(setdiff(all_in, nobg_in), "in_goseq")
  .echogo_validate_support_columns(df, c("in_goseq", bg_in, nobg_in))
  if (!"term_id" %in% names(df)) stop("Legacy table requires term_id.", call. = FALSE)
  if (!"term_name" %in% names(df)) df$term_name <- df$term_id
  if (!"ontology" %in% names(df)) df$ontology <- NA_character_
  goseq_p <- if ("min_pval_goseq" %in% names(df)) suppressWarnings(as.numeric(df$min_pval_goseq)) else rep(NA_real_, nrow(df))
  target_supported <- (df$in_goseq %in% TRUE) & !is.na(goseq_p) & goseq_p <= 0.05
  bg_n <- if (length(bg_in)) rowSums(as.data.frame(df[, bg_in, drop = FALSE]), na.rm = TRUE) else rep(0L, nrow(df))
  bg_total <- length(bg_in)
  profile <- ifelse(
    target_supported & bg_n == 0L, "TARGET_ONLY",
    ifelse(target_supported & bg_n > 0L, "TARGET_PLUS_CONTEXT",
      ifelse(!target_supported & bg_n > 0L, "ALTERNATIVE_CONTEXT", "NO_PRIMARY_SUPPORT")
    )
  )
  out <- tibble::as_tibble(df) %>%
    dplyr::transmute(
      term_id = .data$term_id,
      term_name = .data$term_name,
      ontology = .data$ontology,
      target_goseq_supported = target_supported,
      target_goseq_adjusted_p = goseq_p,
      gprof_custom_context_support_n = as.integer(bg_n),
      gprof_custom_context_support_fraction = if (bg_total) bg_n / bg_total else NA_real_,
      alternative_context_support_n = as.integer(bg_n),
      alternative_context_support_fraction = if (bg_total) bg_n / bg_total else NA_real_,
      queried_context_n = bg_total,
      alternative_queried_context_n = bg_total,
      evidence_profile = profile,
      legacy_wide_table = TRUE
    )
  .echogo_add_display_order(out)
}

#' Add GO depth as ontology metadata
#'
#' Depth is retained for display or ontology inspection only. It is not a
#' component of v0.1.4 evidence classification, recurrence, or ordering.
#' @keywords internal
.echogo_add_go_depth <- function(df, id_col = "term_id") {
  if (!id_col %in% names(df)) return(df)
  raw_id <- vapply(df[[id_col]], function(x) if (length(x)) as.character(x[[1]]) else NA_character_, character(1))
  # str_extract preserves one output position per input row. Unlike a
  # vectorized regmatches/regexpr pair, it cannot recycle a neighbouring GO ID
  # into non-GO or missing rows.
  go_id <- stringr::str_extract(raw_id, "GO:\\d{7}")
  depth <- rep(NA_integer_, nrow(df))
  if (requireNamespace("GO.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
    keys <- unique(go_id[!is.na(go_id)])
    if (length(keys)) {
      maps <- lapply(
        list(GO.db::GOBPANCESTOR, GO.db::GOMFANCESTOR, GO.db::GOCCANCESTOR),
        function(map) AnnotationDbi::mget(keys, map, ifnotfound = NA)
      )
      depth[!is.na(go_id)] <- vapply(go_id[!is.na(go_id)], function(id) {
        ancestors <- unique(unlist(lapply(maps, function(map) map[[id]]), use.names = FALSE))
        as.integer(sum(!is.na(ancestors) & ancestors != "all"))
      }, integer(1))
    }
  }
  if ("depth" %in% names(df)) {
    df$depth <- dplyr::coalesce(suppressWarnings(as.integer(df$depth)), depth)
  } else {
    df$depth <- depth
  }
  df
}

.echogo_normalize_ontology <- function(x) {
  dplyr::case_when(
    x %in% c("GO:BP", "BP") ~ "BP",
    x %in% c("GO:MF", "MF") ~ "MF",
    x %in% c("GO:CC", "CC") ~ "CC",
    TRUE ~ as.character(x)
  )
}

#' Build the canonical scoreless exact-term evidence table
#'
#' The returned table is a complete, provenance-linked evidence landscape. It
#' never combines p-values across contexts. `target_goseq_supported` is derived
#' from the configured source-specific adjusted-p threshold, not from GOseq row
#' presence.
#'
#' @param goseq_file Annotated GOseq result CSV.
#' @param gprofiler_dir g:Profiler output root.
#' @param species_map Named g:Profiler organism code to label mapping.
#' @param output_dir Output directory.
#' @param target_context Optional researcher-designated target g:Profiler code or label.
#' @param context_metadata Optional researcher-supplied context rationale provenance.
#' @param adjusted_p_threshold Source-specific adjusted-p qualification threshold.
#' @param annotation_provenance Optional resolver mapping table.
#' @return Exact-term evidence data frame, with long provenance attached as an attribute.
#' @export
build_consensus_table <- function(
    goseq_file,
    gprofiler_dir = "cross_species_gprofiler",
    species_map = c("hsapiens" = "human", "mmusculus" = "mouse", "drerio" = "zebrafish"),
    output_dir = "consensus_enrichment",
    target_context = NULL,
    context_metadata = NULL,
    adjusted_p_threshold = 0.05,
    annotation_provenance = NULL
) {
  if (!file.exists(goseq_file)) stop("GOseq file not found: ", goseq_file, call. = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  contexts <- .echogo_context_configuration(species_map, target_context, context_metadata)
  goseq_long <- .echogo_read_goseq_long(goseq_file)
  custom_long <- .echogo_read_gprofiler_long(
    gprofiler_dir, species_map,
    mode = "custom_experimental_background", target_context = target_context,
    context_metadata = context_metadata
  )
  exploratory_long <- .echogo_read_gprofiler_long(
    gprofiler_dir, species_map,
    mode = "default_domain_exploratory", target_context = target_context,
    context_metadata = context_metadata
  )
  assembled <- assemble_echogo_evidence(
    goseq_long = goseq_long,
    gprofiler_custom_long = custom_long,
    gprofiler_exploratory_long = exploratory_long,
    context_configuration = contexts,
    adjusted_p_threshold = adjusted_p_threshold
  )
  evidence <- assembled$exact_terms %>%
    dplyr::mutate(
      ontology = .echogo_normalize_ontology(.data$ontology),
      primary_evidence = .data$evidence_profile != "NO_PRIMARY_SUPPORT",
      source_adjusted_p_threshold = adjusted_p_threshold
    )
  evidence <- .echogo_add_go_depth(evidence, id_col = "term_id")

  readr::write_csv(evidence, file.path(output_dir, "term_evidence_exact.csv"))
  readr::write_csv(assembled$source_provenance, file.path(output_dir, "term_source_provenance_long.csv"))
  readr::write_csv(contexts, file.path(output_dir, "context_configuration.csv"))
  # The historical filename remains as a compatibility alias, but the content
  # is scoreless v0.1.4 exact-term evidence.
  openxlsx::write.xlsx(
    evidence,
    file = file.path(output_dir, "consensus_enrichment_results_with_and_without_bg.xlsx"),
    asTable = TRUE
  )
  openxlsx::write.xlsx(
    evidence,
    file = file.path(output_dir, "term_evidence_exact.xlsx"),
    asTable = TRUE
  )
  if (!is.null(annotation_provenance)) {
    .echogo_write_annotation_provenance(annotation_provenance, evidence, output_dir)
  }
  attr(evidence, "source_provenance") <- assembled$source_provenance
  attr(evidence, "context_configuration") <- contexts
  evidence
}

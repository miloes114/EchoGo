# File-aware g:Profiler scientific-input resolution ---------------------------
#
# Plot/output directories may exist even when the scientific tables remain in a
# legacy compatibility directory. Readers therefore resolve expected files,
# never directory existence alone.

.echogo_gprofiler_result_candidates <- function(gprofiler_dir, mode, context_label) {
  custom <- identical(mode, "custom_experimental_background")
  suffix <- if (custom) "with_bg" else "nobg"
  roots <- .echogo_gprofiler_mode_paths(gprofiler_dir, mode)

  unlist(lapply(roots, function(root) {
    stem <- file.path(root, paste0("gprofiler_", context_label, "_", suffix))
    c(
      paste0(stem, ".csv"),
      paste0(stem, "_enrichment.csv")
    )
  }), use.names = FALSE)
}

.echogo_find_gprofiler_result_file <- function(gprofiler_dir, mode, context_label) {
  candidates <- .echogo_gprofiler_result_candidates(gprofiler_dir, mode, context_label)
  hit <- candidates[file.exists(candidates)][1]
  if (length(hit)) normalizePath(hit, winslash = "/", mustWork = FALSE) else NA_character_
}

.echogo_gprofiler_sidecar_from_csv <- function(csv, suffix) {
  if (is.null(csv) || is.na(csv) || !nzchar(csv)) return(NA_character_)
  stem <- sub("_enrichment\\.csv$", "", csv, ignore.case = TRUE)
  stem <- sub("\\.csv$", "", stem, ignore.case = TRUE)
  paste0(stem, suffix)
}

# Override the earlier directory-first implementation with a file-aware reader.
.echogo_read_gprofiler_long <- function(gprofiler_dir, species_map,
                                        mode = c("custom_experimental_background", "default_domain_exploratory"),
                                        target_context = NULL,
                                        context_metadata = NULL) {
  mode <- match.arg(mode)
  contexts <- .echogo_context_configuration(species_map, target_context, context_metadata)
  rows <- vector("list", nrow(contexts))

  for (i in seq_len(nrow(contexts))) {
    label <- contexts$context_label[[i]]
    csv <- .echogo_find_gprofiler_result_file(gprofiler_dir, mode, label)
    if (is.na(csv)) next

    tab <- readr::read_csv(csv, show_col_types = FALSE)
    if (!all(c("term_id", "p_value") %in% names(tab))) {
      stop("g:Profiler CSV lacks term_id or p_value: ", csv, call. = FALSE)
    }

    metadata <- .echogo_gprofiler_sidecar_from_csv(csv, "_metadata.json")
    scalar <- function(field, fallback = NA) {
      value <- .echogo_read_json_scalar(metadata, field)
      if (length(value) == 0L || all(is.na(value))) fallback else value[[1]]
    }

    rows[[i]] <- tibble::tibble(
      term_id = as.character(tab$term_id),
      term_name = if ("term_name" %in% names(tab)) as.character(tab$term_name) else NA_character_,
      ontology = if ("source" %in% names(tab)) as.character(tab$source) else NA_character_,
      source_type = "gprofiler",
      context_code = contexts$context_code[[i]],
      context_label = label,
      context_role = contexts$context_role[[i]],
      context_rationale = contexts$context_rationale[[i]],
      selection_basis = contexts$selection_basis[[i]],
      background_mode = mode,
      domain_scope = if (identical(mode, "custom_experimental_background")) "custom" else "annotated",
      source_adjusted_p = suppressWarnings(as.numeric(tab$p_value)),
      fold_enrichment = if ("fold_enrichment" %in% names(tab)) suppressWarnings(as.numeric(tab$fold_enrichment)) else NA_real_,
      contributing_genes = if ("intersection" %in% names(tab)) as.character(tab$intersection) else NA_character_,
      intersection_size = if ("intersection_size" %in% names(tab)) suppressWarnings(as.numeric(tab$intersection_size)) else NA_real_,
      query_size = if ("query_size" %in% names(tab)) suppressWarnings(as.numeric(tab$query_size)) else NA_real_,
      effective_domain_size = if ("effective_domain_size" %in% names(tab)) suppressWarnings(as.numeric(tab$effective_domain_size)) else NA_real_,
      submitted_foreground_n = suppressWarnings(as.numeric(scalar("submitted_foreground_count"))),
      submitted_background_n = suppressWarnings(as.numeric(scalar("submitted_background_count"))),
      effective_query_n = suppressWarnings(as.numeric(scalar("effective_query_size"))),
      effective_background_or_domain_n = suppressWarnings(as.numeric(scalar("effective_domain_size"))),
      provenance_file = csv,
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(rows)
}

# Locate cached exact vectors without allowing an empty plot directory to mask
# a populated compatibility directory.
.echogo_find_cached_gprofiler_vectors <- function(gprofiler_dir, mode = "custom_experimental_background") {
  roots <- .echogo_gprofiler_mode_paths(gprofiler_dir, mode)
  for (root in roots) {
    if (!dir.exists(root)) next
    queries <- list.files(root, pattern = "_query\\.txt$", full.names = TRUE)
    backgrounds <- list.files(root, pattern = "_background\\.txt$", full.names = TRUE)
    if (length(queries) && length(backgrounds)) {
      return(list(directory = root, queries = queries, backgrounds = backgrounds))
    }
  }
  list(directory = NA_character_, queries = character(), backgrounds = character())
}

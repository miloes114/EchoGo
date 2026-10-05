# Cache-directory ordering for pipeline compatibility -------------------------
#
# The evidence reader resolves exact result files in canonical-first order. The
# older cached-vector validation in run_echogo_pipeline(), however, still asks
# for a directory before looking for *_query.txt and *_background.txt files.
# Reorder only that directory list so a plot-only canonical directory can never
# mask the directory that actually contains the frozen scientific vectors.

.echogo_gprofiler_mode_paths <- function(gprofiler_dir, mode) {
  paths <- if (identical(mode, "custom_experimental_background")) {
    c(
      file.path(gprofiler_dir, "custom_experimental_background"),
      file.path(gprofiler_dir, "bg"),
      file.path(gprofiler_dir, "with_custom_background")
    )
  } else {
    c(
      file.path(gprofiler_dir, "default_domain_exploratory"),
      file.path(gprofiler_dir, "nobg"),
      file.path(gprofiler_dir, "no_background_genome_wide")
    )
  }

  vector_score <- vapply(paths, function(root) {
    if (!dir.exists(root)) return(0L)
    has_query <- length(list.files(root, pattern = "_query\\.txt$", full.names = FALSE)) > 0L
    has_background <- length(list.files(root, pattern = "_background\\.txt$", full.names = FALSE)) > 0L
    if (has_query && has_background) return(2L)
    1L
  }, integer(1))

  paths[order(-vector_score, seq_along(paths))]
}

# Exact scientific result-table resolution remains canonical-first by FILE,
# independent of the compatibility ordering above.
.echogo_gprofiler_result_candidates <- function(gprofiler_dir, mode, context_label) {
  custom <- identical(mode, "custom_experimental_background")
  suffix <- if (custom) "with_bg" else "nobg"
  roots <- if (custom) {
    c(
      file.path(gprofiler_dir, "custom_experimental_background"),
      file.path(gprofiler_dir, "bg"),
      file.path(gprofiler_dir, "with_custom_background")
    )
  } else {
    c(
      file.path(gprofiler_dir, "default_domain_exploratory"),
      file.path(gprofiler_dir, "nobg"),
      file.path(gprofiler_dir, "no_background_genome_wide")
    )
  }

  unlist(lapply(roots, function(root) {
    stem <- file.path(root, paste0("gprofiler_", context_label, "_", suffix))
    c(
      paste0(stem, ".csv"),
      paste0(stem, "_enrichment.csv")
    )
  }), use.names = FALSE)
}

# Source-aware g:Profiler report status -------------------------------------
# Keeps genuine biological empty states distinct from missing plot assets.

.echogo_gprofiler_report_status <- function(sources, context_label, ontology,
                                            exploratory = FALSE) {
  ontology_code <- toupper(as.character(ontology[[1]]))
  if (!ontology_code %in% c("BP", "MF", "CC")) stop("ontology must be BP, MF, or CC", call. = FALSE)
  mode <- if (isTRUE(exploratory)) "default_domain_exploratory" else "custom_experimental_background"

  if (is.null(sources) || !is.data.frame(sources) || !nrow(sources)) {
    return(tibble::tibble(
      status = "NO_SOURCE_ROWS",
      source_row_n = 0L,
      qualifying_term_n = 0L,
      exploratory_term_n = 0L,
      message = "No source-level g:Profiler result rows were recorded for this context and ontology."
    ))
  }

  ont <- if ("ontology" %in% names(sources)) sub("^GO:", "", toupper(as.character(sources$ontology))) else rep(NA_character_, nrow(sources))
  bg <- if ("background_mode" %in% names(sources)) as.character(sources$background_mode) else rep(NA_character_, nrow(sources))
  ctx <- if ("context_label" %in% names(sources)) as.character(sources$context_label) else rep(NA_character_, nrow(sources))
  typ <- if ("source_type" %in% names(sources)) as.character(sources$source_type) else rep(NA_character_, nrow(sources))

  # Compatibility aliases are accepted only for locating the same run mode.
  bg_ok <- if (isTRUE(exploratory)) {
    bg %in% c("default_domain_exploratory", "no_background_genome_wide", "no_background")
  } else {
    bg %in% c("custom_experimental_background", "with_custom_background", "custom_background")
  }

  d <- sources[typ == "gprofiler" & ctx == context_label & ont == ontology_code & bg_ok, , drop = FALSE]
  source_n <- nrow(d)
  qual_n <- if (source_n && "source_qualifies" %in% names(d)) {
    length(unique(as.character(d$term_id[d$source_qualifies %in% TRUE])))
  } else 0L
  exploratory_n <- NA_integer_

  if (!source_n) {
    status <- "NO_SOURCE_ROWS"
    message <- "No source-level g:Profiler result rows were recorded for this context and ontology."
  } else if (isTRUE(exploratory)) {
    exploratory_n <- if ("term_id" %in% names(d)) length(unique(stats::na.omit(as.character(d$term_id)))) else 0L
    if (exploratory_n) {
      status <- "EXPLORATORY_PLOT_EXPECTED"
      message <- paste0(
        exploratory_n,
        " default-domain exploratory GO-term result(s) are available. Source statistics remain context-local and are not primary EchoGO support."
      )
    } else {
      status <- "NO_EXPLORATORY_RESULTS"
      message <- "No default-domain exploratory GO-term results were recorded for this context and ontology."
    }
  } else if (!qual_n) {
    status <- "NO_QUALIFYING_TERMS"
    message <- "No source-qualified matched-background GO terms were available for this context and ontology."
  } else {
    status <- "FIGURE_EXPECTED"
    message <- paste0(qual_n, " source-qualified matched-background GO term", if (qual_n == 1L) "" else "s", " are present; a source-level figure is expected.")
  }

  tibble::tibble(
    status = status,
    source_row_n = as.integer(source_n),
    qualifying_term_n = if (isTRUE(exploratory)) NA_integer_ else as.integer(qual_n),
    exploratory_term_n = if (isTRUE(exploratory)) as.integer(exploratory_n) else NA_integer_,
    message = message
  )
}

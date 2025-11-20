#' Plotting utilities for EchoGO
#'
#' Generic plots used across the EchoGO workflow:
#' - Lollipop plots for top enriched terms
#' - Density plots of EQI / enrichment by method and ontology
#' - Term-origin stacked bar plots
#'
#' @name echogo_plotting_utils
#' @keywords internal
#' @noRd
NULL

#' Generic lollipop plot for enrichment results
#'
#' Creates a lollipop plot of top ranked GO terms or KEGG pathways using fold enrichment or
#' consensus scores, colored and sized by optional variables (e.g., -log10(p), origin).
#'
#' @param df A data frame of enrichment results.
#' @param score_col Name of the column used for ranking (e.g. "fold_enrichment", "consensus_score").
#' @param label Plot title.
#' @param output_path File path to save the plot.
#' @param top_n Number of top terms to show.
#' @param color_by Column name to use for coloring (e.g. "-log10(p)", "origin").
#' @param size_by Column name to use for point size (e.g. "intersection_size").
#' @export
#' @importFrom stringr str_trunc
plot_lollipop <- function(df, score_col, label, output_path, top_n = 50,
                          color_by = NULL, size_by = NULL) {

  # Basic checks
  if (!is.data.frame(df) || !nrow(df)) {
    stop("plot_lollipop(): `df` is empty or not a data.frame.")
  }
  if (!score_col %in% names(df)) {
    stop("plot_lollipop(): `score_col` '", score_col, "' not found in df.")
  }
  if (!is.null(color_by) && !color_by %in% names(df)) {
    warning("plot_lollipop(): `color_by` column '", color_by, "' not found; ignoring.")
    color_by <- NULL
  }
  if (!is.null(size_by) && !size_by %in% names(df)) {
    warning("plot_lollipop(): `size_by` column '", size_by, "' not found; ignoring.")
    size_by <- NULL
  }

  # Pick the first label column that exists (consensus tables don't have `term`)
  pick_label <- function(d) {
    for (nm in c("term_name", "label", "term_id")) {
      if (nm %in% names(d)) return(as.character(d[[nm]]))
    }
    rep(NA_character_, nrow(d))
  }

  df <- df %>%
    dplyr::filter(!is.na(.data[[score_col]])) %>%
    dplyr::arrange(dplyr::desc(.data[[score_col]])) %>%
    dplyr::slice_head(n = top_n)

  df$term_display <- pick_label(df)
  if (all(is.na(df$term_display))) {
    stop("plot_lollipop(): couldn't find any of term_name/label/term_id to use as labels.")
  }

  df <- df %>%
    dplyr::mutate(
      term_display = stringr::str_trunc(term_display, 40),
      term_display = factor(term_display, levels = rev(unique(term_display)))
    )

  # Build aesthetics programmatically (avoid conditionals inside aes())
  p <- ggplot2::ggplot(df, ggplot2::aes_string(x = "term_display", y = score_col))

  seg_map <- if (is.null(color_by)) {
    ggplot2::aes_string(xend = "term_display", y = 0, yend = score_col)
  } else {
    ggplot2::aes_string(xend = "term_display", y = 0, yend = score_col, colour = color_by)
  }

  pt_map <- if (is.null(color_by) && is.null(size_by)) {
    ggplot2::aes()
  } else if (is.null(color_by)) {
    ggplot2::aes_string(size = size_by)
  } else if (is.null(size_by)) {
    ggplot2::aes_string(colour = color_by)
  } else {
    ggplot2::aes_string(colour = color_by, size = size_by)
  }

  p <- p +
    ggplot2::geom_segment(mapping = seg_map, linewidth = 1.2) +
    ggplot2::geom_point(mapping = pt_map) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title  = label,
      x      = "Term",
      y      = score_col,
      colour = if (is.null(color_by)) NULL else color_by,
      size   = if (is.null(size_by))   NULL else size_by
    ) +
    ggplot2::theme_minimal(base_size = 13)

  ggplot2::ggsave(output_path, p, width = 10, height = 7)
  invisible(p)
}

#' Density plot of EQI or enrichment values by method and ontology
#'
#' @param df A long-format data frame with columns:
#'   - \code{metric} column (\code{"EQI"} or \code{"enrichment"}) present by name
#'   - \code{method} (e.g., "GOseq", "Consensus")
#'   - \code{ontology} (BP/MF/CC)
#' @param metric One of \code{"EQI"} or \code{"enrichment"} (default: "EQI").
#' @param output_path File path to save the PDF.
#'
#' @return (Invisibly) the ggplot object; writes the PDF.
#' @export
#' @importFrom ggplot2 ggplot aes geom_density facet_wrap labs scale_fill_manual theme_minimal ggsave
#' @importFrom rlang .data
plot_density_by_ontology <- function(df, metric = "EQI", output_path) {
  if (!metric %in% c("EQI", "enrichment")) {
    stop("plot_density_by_ontology(): metric must be 'EQI' or 'enrichment'.")
  }
  req <- c("method", "ontology", metric)
  miss <- setdiff(req, names(df))
  if (length(miss)) {
    stop("plot_density_by_ontology(): missing columns: ", paste(miss, collapse = ", "))
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[metric]], fill = .data$method)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::facet_wrap(~ontology, scales = "free") +
    ggplot2::labs(
      title = paste(metric, "Distribution by Method and Ontology"),
      x = metric, y = "Density"
    ) +
    ggplot2::scale_fill_manual(values = c("GOseq" = "steelblue", "Consensus" = "darkgreen")) +
    ggplot2::theme_minimal(base_size = 13)

  ggplot2::ggsave(output_path, plot = p, width = 10, height = 6)
  invisible(p)
}

#' Bar plot of term origin frequencies
#'
#' @param df Consensus enrichment data frame with columns \code{origin} and \code{ontology}.
#' @param output_path File path to save the plot.
#'
#' @return (Invisibly) the ggplot object; writes the figure.
#' @export
#' @importFrom ggplot2 ggplot aes geom_col labs theme_minimal ggsave
#' @importFrom dplyr count
#' @importFrom rlang .data
plot_term_origin_distribution <- function(df, output_path) {
  if (!all(c("origin", "ontology") %in% names(df))) {
    stop("plot_term_origin_distribution(): df must have 'origin' and 'ontology'.")
  }

  p <- df %>%
    dplyr::count(.data$origin, .data$ontology, name = "n") %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$ontology, y = .data$n, fill = .data$origin)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::labs(
      title = "GO Term Origin by Ontology",
      x = "Ontology", y = "Count of Terms", fill = "Origin"
    ) +
    ggplot2::theme_minimal(base_size = 13)

  ggplot2::ggsave(output_path, p, width = 8, height = 6)
  invisible(p)
}

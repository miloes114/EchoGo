# EchoGO key-findings report layer ------------------------------------------
#
# This file provides a compact, scoreless summary view of the exact-term
# evidence table.  It intentionally does not alter evidence membership or
# source-local statistics.  The only cross-reference quantity displayed is the
# already-defined alternative-context recurrence n/N.

.echogo_key_findings_gene_count <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  unname(vapply(x, function(one) {
    genes <- trimws(unlist(strsplit(one, "[,;/|]")))
    genes <- unique(genes[nzchar(genes)])
    length(genes)
  }, integer(1)))
}

.echogo_key_findings_group_labels <- c(
  TARGET_PLUS_CONTEXT = "Target + context",
  TARGET_ONLY = "Target only",
  ALTERNATIVE_CONTEXT = "Context-derived hypotheses"
)

.echogo_key_findings_cap_for_profile <- function(profile,
                                                  max_target_recovered,
                                                  max_target_only,
                                                  max_hypotheses) {
  switch(
    as.character(profile),
    TARGET_PLUS_CONTEXT = as.integer(max_target_recovered),
    TARGET_ONLY = as.integer(max_target_only),
    ALTERNATIVE_CONTEXT = as.integer(max_hypotheses),
    0L
  )
}

#' Prepare a compact, scoreless EchoGO key-findings display table
#'
#' This is a presentation helper.  It never creates a composite statistic.
#' Target-supported streams use target GOseq adjusted p-values only as a
#' within-stream display tie-breaker.  Alternative-context hypotheses are
#' ordered by exact-term recurrence and then by the existing deterministic
#' display order/GO ID; g:Profiler p-values are never compared across contexts.
#'
#' @param evidence Exact-term EchoGO evidence table.
#' @param ontology Optional ontology (`BP`, `MF`, or `CC`).
#' @param max_target_recovered Maximum displayed TARGET_PLUS_CONTEXT terms.
#' @param max_target_only Maximum displayed TARGET_ONLY terms.
#' @param max_hypotheses Maximum displayed ALTERNATIVE_CONTEXT terms.
#' @return A tibble containing only terms selected for the compact display.
#' @keywords internal
.echogo_key_findings_data <- function(evidence,
                                      ontology = NULL,
                                      max_target_recovered = 8L,
                                      max_target_only = 4L,
                                      max_hypotheses = 8L) {
  if (is.null(evidence) || !is.data.frame(evidence) || !nrow(evidence)) {
    return(tibble::tibble())
  }

  required <- c(
    "term_id", "term_name", "ontology", "evidence_profile",
    "target_goseq_adjusted_p", "alternative_context_support_n",
    "alternative_queried_context_n", "alternative_context_support_fraction",
    "contributing_genes"
  )
  missing <- setdiff(required, names(evidence))
  if (length(missing)) {
    stop(
      "Key-findings display requires columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!"display_order" %in% names(evidence)) {
    evidence <- .echogo_add_display_order(evidence)
  }

  allowed_profiles <- c("TARGET_PLUS_CONTEXT", "TARGET_ONLY", "ALTERNATIVE_CONTEXT")
  profile_order <- c("TARGET_PLUS_CONTEXT", "TARGET_ONLY", "ALTERNATIVE_CONTEXT")
  profile_levels <- unname(.echogo_key_findings_group_labels[profile_order])

  d <- tibble::as_tibble(evidence) |>
    dplyr::filter(.data$evidence_profile %in% allowed_profiles)

  if (!is.null(ontology)) {
    ontology <- toupper(as.character(ontology[[1]]))
    d <- dplyr::filter(d, .data$ontology %in% c(ontology, paste0("GO:", ontology)))
  }
  if (!nrow(d)) return(tibble::tibble())

  d <- d |>
    dplyr::mutate(
      ontology = sub("^GO:", "", toupper(as.character(.data$ontology))),
      alternative_context_support_n = dplyr::coalesce(
        suppressWarnings(as.integer(.data$alternative_context_support_n)), 0L
      ),
      alternative_queried_context_n = dplyr::coalesce(
        suppressWarnings(as.integer(.data$alternative_queried_context_n)), 0L
      ),
      alternative_context_support_fraction = dplyr::if_else(
        .data$alternative_queried_context_n > 0L,
        .data$alternative_context_support_n / .data$alternative_queried_context_n,
        0
      ),
      target_goseq_adjusted_p = suppressWarnings(as.numeric(.data$target_goseq_adjusted_p)),
      contributing_gene_n = .echogo_key_findings_gene_count(.data$contributing_genes),
      finding_group = unname(.echogo_key_findings_group_labels[.data$evidence_profile]),
      profile_order = match(.data$evidence_profile, profile_order),
      support_label = dplyr::if_else(
        .data$alternative_queried_context_n > 0L,
        paste0(.data$alternative_context_support_n, "/", .data$alternative_queried_context_n),
        "no alternative contexts"
      )
    )

  # The ordering is intentionally transparent and stream-specific.
  # TARGET_PLUS_CONTEXT: broader reference recovery first, then the target's
  # own GOseq FDR as a source-local tie-breaker.
  # TARGET_ONLY: the target's own GOseq FDR, because recurrence is necessarily 0.
  # ALTERNATIVE_CONTEXT: recurrence, then pre-existing deterministic order/GO ID.
  plus <- d |>
    dplyr::filter(.data$evidence_profile == "TARGET_PLUS_CONTEXT") |>
    dplyr::arrange(
      dplyr::desc(.data$alternative_context_support_n),
      .data$target_goseq_adjusted_p,
      .data$display_order,
      .data$term_id
    ) |>
    dplyr::slice_head(n = as.integer(max_target_recovered))

  target_only <- d |>
    dplyr::filter(.data$evidence_profile == "TARGET_ONLY") |>
    dplyr::arrange(
      .data$target_goseq_adjusted_p,
      .data$display_order,
      .data$term_id
    ) |>
    dplyr::slice_head(n = as.integer(max_target_only))

  hypothesis <- d |>
    dplyr::filter(.data$evidence_profile == "ALTERNATIVE_CONTEXT") |>
    dplyr::arrange(
      dplyr::desc(.data$alternative_context_support_n),
      .data$display_order,
      .data$term_id
    ) |>
    dplyr::slice_head(n = as.integer(max_hypotheses))

  out <- dplyr::bind_rows(plus, target_only, hypothesis) |>
    dplyr::mutate(
      finding_group = factor(
        .data$finding_group,
        levels = profile_levels
      )
    ) |>
    dplyr::arrange(.data$profile_order, dplyr::desc(.data$alternative_context_support_n),
                   .data$target_goseq_adjusted_p, .data$display_order, .data$term_id) |>
    dplyr::mutate(key_findings_display_order = dplyr::row_number())

  out
}

#' Plot the compact EchoGO Key Findings Map
#'
#' The x-axis is the already-defined fraction of selected alternative reference
#' species that recover each exact GO term under matched-background g:Profiler.
#' It is not an EchoGO significance score.  Rows are split into target-supported
#' findings and new hypotheses so that recurrence is never interpreted as a
#' replacement for target evidence.
#'
#' @param evidence Exact-term EchoGO evidence table.
#' @param ontology Optional ontology (`BP`, `MF`, or `CC`).
#' @param max_target_recovered Maximum displayed TARGET_PLUS_CONTEXT terms.
#' @param max_target_only Maximum displayed TARGET_ONLY terms.
#' @param max_hypotheses Maximum displayed ALTERNATIVE_CONTEXT terms.
#' @param wrap_width Width used to wrap GO term labels.
#' @return A ggplot or NULL when there is nothing to display.
#' @keywords internal
plot_echogo_key_findings_map <- function(evidence,
                                         ontology = NULL,
                                         max_target_recovered = 8L,
                                         max_target_only = 4L,
                                         max_hypotheses = 8L,
                                         wrap_width = 44L) {
  d <- .echogo_key_findings_data(
    evidence = evidence,
    ontology = ontology,
    max_target_recovered = max_target_recovered,
    max_target_only = max_target_only,
    max_hypotheses = max_hypotheses
  )
  if (!nrow(d)) return(NULL)

  d <- d |>
    dplyr::mutate(
      term_label = paste0(
        stringr::str_wrap(.data$term_name, width = wrap_width),
        "\n", .data$term_id
      ),
      point_fraction = pmax(0, pmin(1, .data$alternative_context_support_fraction))
    )

  # Keep visual rows in the exact prepared display order.  Reversing factor
  # levels makes the first selected term appear at the top of each facet.
  d$term_label <- factor(d$term_label, levels = rev(unique(d$term_label)))

  title_ontology <- if (is.null(ontology)) {
    "EchoGO Key Findings Map"
  } else {
    paste0("EchoGO Key Findings Map - ", toupper(as.character(ontology[[1]])))
  }

  has_other_refs <- any(d$alternative_queried_context_n > 0L, na.rm = TRUE)
  subtitle_text <- if (has_other_refs) {
    paste(
      "Position shows how many selected alternative annotation contexts recover the exact GO term.",
      "Rows are selected for display only; no EchoGO significance score is calculated."
    )
  } else {
    paste(
      "No alternative annotation contexts were declared for this run.",
      "The map therefore shows target findings without a context-recovery axis."
    )
  }

  p <- ggplot2::ggplot(d, ggplot2::aes(y = .data$term_label)) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = .data$point_fraction,
        yend = .data$term_label,
        colour = .data$evidence_profile
      ),
      linewidth = 1.05,
      alpha = .65
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data$point_fraction,
        size = pmax(1L, .data$contributing_gene_n),
        fill = .data$evidence_profile
      ),
      shape = 21,
      stroke = .75,
      colour = "#314854",
      alpha = .95
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = .data$point_fraction,
        label = .data$support_label
      ),
      nudge_x = .035,
      hjust = 0,
      size = 3.15,
      colour = "#435665"
    ) +
    ggplot2::facet_grid(
      finding_group ~ .,
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    ggplot2::scale_colour_manual(values = .echogo_profile_colours, guide = "none") +
    ggplot2::scale_fill_manual(
      values = .echogo_profile_colours,
      breaks = c("TARGET_PLUS_CONTEXT", "TARGET_ONLY", "ALTERNATIVE_CONTEXT"),
      labels = c(
        "Target + context",
        "Target only",
        "Context-derived hypotheses"
      ),
      name = "Where the finding comes from"
    ) +
    ggplot2::scale_size_continuous(
      range = c(3.0, 8.5),
      name = "Unique contributing genes\nreported across sources"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1.16),
      breaks = c(0, .25, .5, .75, 1),
      labels = c("0", "25%", "50%", "75%", "100%")
    ) +
    ggplot2::labs(
      x = "Annotation contexts with qualifying evidence",
      y = NULL,
      title = title_ontology,
      subtitle = subtitle_text,
      caption = paste(
        "For target-supported terms, the target GOseq result remains the experimental anchor.",
        "Context recovery is descriptive n/N support, not evidence of conservation or independent replication."
      )
    ) +
    ggplot2::theme_minimal(base_size = 11.5) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(
        angle = 0,
        face = "bold",
        size = 9.4,
        colour = "#193247"
      ),
      axis.text.y = ggplot2::element_text(size = 9.4, colour = "#263d4b"),
      plot.title = ggplot2::element_text(face = "bold", colour = "#0d4258"),
      plot.subtitle = ggplot2::element_text(colour = "#526574"),
      plot.caption = ggplot2::element_text(hjust = 0, colour = "#60727e"),
      legend.position = "bottom"
    )

  p
}

#' Prepare a compact RRvGO theme-level recovery summary
#'
#' This helper summarizes an already-computed RRvGO reduced-term table.  It does
#' not rerun semantic similarity and does not score clusters.  Each cluster is
#' represented by the RRvGO representative term with the highest existing
#' non-inferential representative order (`score`).
#'
#' @param reduced_terms One `rrvgo_*_clusters.csv` table.
#' @param evidence Exact-term EchoGO evidence table.
#' @param sources Long-form EchoGO source provenance table.
#' @param max_clusters Maximum clusters to display.
#' @return A compact theme-level tibble.
#' @keywords internal
.echogo_theme_summary_data <- function(reduced_terms, evidence, sources,
                                       max_clusters = 10L) {
  if (is.null(reduced_terms) || !is.data.frame(reduced_terms) || !nrow(reduced_terms)) {
    return(tibble::tibble())
  }
  required_rr <- c("go", "term", "cluster", "score")
  missing_rr <- setdiff(required_rr, names(reduced_terms))
  if (length(missing_rr)) {
    stop("RRvGO theme summary requires columns: ", paste(missing_rr, collapse = ", "), call. = FALSE)
  }
  if (is.null(evidence) || !is.data.frame(evidence) || !nrow(evidence)) return(tibble::tibble())
  if (is.null(sources) || !is.data.frame(sources)) sources <- tibble::tibble()

  rr <- tibble::as_tibble(reduced_terms) |>
    dplyr::mutate(go = trimws(as.character(.data$go)))

  ev <- tibble::as_tibble(evidence) |>
    dplyr::select(
      .data$term_id,
      .data$alternative_queried_context_n,
      .data$contributing_genes,
      .data$evidence_profile
    )

  reps <- rr |>
    dplyr::group_by(.data$cluster) |>
    dplyr::arrange(dplyr::desc(.data$score), .data$go, .by_group = TRUE) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      cluster = .data$cluster,
      representative_go = .data$go,
      representative_term = .data$term,
      representative_order = .data$score
    )

  membership <- rr |>
    dplyr::select(.data$cluster, term_id = .data$go) |>
    dplyr::distinct() |>
    dplyr::left_join(ev, by = "term_id")

  alt_support <- if (nrow(sources)) {
    sources |>
      dplyr::filter(
        .data$source_type == "gprofiler",
        .data$background_mode == "custom_experimental_background",
        .data$context_role == "ALTERNATIVE",
        .data$source_qualifies %in% TRUE
      ) |>
      dplyr::select(.data$term_id, .data$context_label) |>
      dplyr::distinct() |>
      dplyr::inner_join(rr |> dplyr::select(.data$cluster, term_id = .data$go) |> dplyr::distinct(),
                        by = "term_id") |>
      dplyr::group_by(.data$cluster) |>
      dplyr::summarise(reference_support_n = dplyr::n_distinct(.data$context_label), .groups = "drop")
  } else {
    tibble::tibble(cluster = unique(rr$cluster), reference_support_n = 0L)
  }

  summary <- membership |>
    dplyr::group_by(.data$cluster) |>
    dplyr::summarise(
      exact_term_n = dplyr::n_distinct(.data$term_id),
      contributing_genes = paste(stats::na.omit(.data$contributing_genes), collapse = ";"),
      alternative_queried_context_n = suppressWarnings(max(.data$alternative_queried_context_n, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      alternative_queried_context_n = dplyr::if_else(
        is.infinite(.data$alternative_queried_context_n), 0, .data$alternative_queried_context_n
      ),
      contributing_gene_n = .echogo_key_findings_gene_count(.data$contributing_genes)
    ) |>
    dplyr::left_join(alt_support, by = "cluster") |>
    dplyr::mutate(
      reference_support_n = dplyr::coalesce(.data$reference_support_n, 0L),
      reference_support_fraction = dplyr::if_else(
        .data$alternative_queried_context_n > 0,
        .data$reference_support_n / .data$alternative_queried_context_n,
        0
      ),
      support_label = dplyr::if_else(
        .data$alternative_queried_context_n > 0,
        paste0(.data$reference_support_n, "/", .data$alternative_queried_context_n),
        "no alternative contexts"
      )
    ) |>
    dplyr::left_join(reps, by = "cluster") |>
    dplyr::arrange(dplyr::desc(.data$exact_term_n), dplyr::desc(.data$reference_support_n),
                   dplyr::desc(.data$representative_order), .data$representative_go) |>
    dplyr::slice_head(n = as.integer(max_clusters))

  summary
}

#' Plot a compact biological-theme summary from an existing RRvGO product
#'
#' @param reduced_terms One RRvGO reduced-term table.
#' @param evidence Exact-term EchoGO evidence table.
#' @param sources Long-form source provenance.
#' @param product_label Human-facing product label.
#' @param max_clusters Maximum clusters displayed.
#' @param wrap_width Label wrapping width.
#' @return A ggplot or NULL.
#' @keywords internal
plot_echogo_theme_summary <- function(reduced_terms, evidence, sources,
                                      product_label = "Functional themes",
                                      max_clusters = 10L,
                                      wrap_width = 42L) {
  d <- .echogo_theme_summary_data(
    reduced_terms = reduced_terms,
    evidence = evidence,
    sources = sources,
    max_clusters = max_clusters
  )
  if (!nrow(d)) return(NULL)

  d <- d |>
    dplyr::mutate(
      theme_label = paste0(
        stringr::str_wrap(.data$representative_term, width = wrap_width),
        "\n", .data$representative_go
      ),
      point_fraction = pmax(0, pmin(1, .data$reference_support_fraction))
    )
  d$theme_label <- factor(d$theme_label, levels = rev(d$theme_label))

  ggplot2::ggplot(d, ggplot2::aes(y = .data$theme_label)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = .data$point_fraction, yend = .data$theme_label),
      linewidth = 1.1, colour = "#89aeb5"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$point_fraction, size = .data$exact_term_n),
      shape = 21, fill = "#2a8f88", colour = "#314854", stroke = .75
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = .data$point_fraction, label = .data$support_label),
      nudge_x = .035, hjust = 0, size = 3.15, colour = "#435665"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1.16), breaks = c(0, .25, .5, .75, 1),
      labels = c("0", "25%", "50%", "75%", "100%")
    ) +
    ggplot2::scale_size_continuous(range = c(3.4, 9), name = "Exact GO terms\nin semantic theme") +
    ggplot2::labs(
      x = "Annotation contexts contributing at least one exact term to the theme",
      y = NULL,
      title = product_label,
      subtitle = "RRvGO themes are summarized after exact-term evidence is fixed; this is not a cluster score.",
      caption = "Theme labels use the existing RRvGO representative term. Context recovery remains descriptive."
    ) +
    ggplot2::theme_minimal(base_size = 11.5) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 9.5),
      plot.title = ggplot2::element_text(face = "bold", colour = "#0d4258"),
      legend.position = "bottom"
    )
}

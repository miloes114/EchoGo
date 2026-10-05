# Universal report-facing visualizations -------------------------------------
#
# These helpers are presentation views of already assembled EchoGO evidence.
# They must not alter evidence membership, source statistics, recurrence, or
# semantic products.

.echogo_profile_colours <- c(
  TARGET_ONLY = "#176b87",
  TARGET_PLUS_CONTEXT = "#2a8f88",
  ALTERNATIVE_CONTEXT = "#d99a26",
  NO_PRIMARY_SUPPORT = "#b5c0c6"
)

.echogo_profile_labels <- c(
  TARGET_ONLY = "Target only",
  TARGET_PLUS_CONTEXT = "Target + context",
  ALTERNATIVE_CONTEXT = "Context-derived hypotheses",
  NO_PRIMARY_SUPPORT = "No primary support"
)

.echogo_safe_count <- function(x) sum(x %in% TRUE, na.rm = TRUE)

.echogo_mapping_summary <- function(mapping) {
  if (is.null(mapping) || !is.data.frame(mapping) || !nrow(mapping)) {
    return(tibble::tibble())
  }

  portable_col <- intersect(c("resolved_name", "portable_name"), names(mapping))[1]
  portable <- if (!is.na(portable_col)) {
    trimws(as.character(mapping[[portable_col]]))
  } else {
    rep(NA_character_, nrow(mapping))
  }
  missing_portable <- is.na(portable) | !nzchar(portable)

  duplicate_rows <- if ("exclusion_reason" %in% names(mapping)) {
    grepl("duplicate", as.character(mapping$exclusion_reason), ignore.case = TRUE)
  } else {
    rep(FALSE, nrow(mapping))
  }

  tibble::tibble(
    tested_entities = .echogo_safe_count(mapping$tested),
    submitted_background_names = .echogo_safe_count(mapping$included_background),
    tested_without_portable_name = sum(mapping$tested %in% TRUE & missing_portable, na.rm = TRUE),
    tested_duplicate_rows_collapsed = sum(mapping$tested %in% TRUE & duplicate_rows, na.rm = TRUE),
    significant_entities = .echogo_safe_count(mapping$significant),
    submitted_foreground_names = .echogo_safe_count(mapping$included_foreground)
  )
}

#' Plot representability from the target experiment to comparative submission
#' @keywords internal
plot_echogo_representability <- function(mapping) {
  s <- .echogo_mapping_summary(mapping)
  if (!nrow(s)) return(NULL)

  p <- ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = .2, xmax = 2.8, ymin = 2.25, ymax = 3.25,
                      fill = "#eef7f8", colour = "#c8d7df", linewidth = .5) +
    ggplot2::annotate("rect", xmin = 5.2, xmax = 7.8, ymin = 2.25, ymax = 3.25,
                      fill = "#eef7f8", colour = "#c8d7df", linewidth = .5) +
    ggplot2::annotate("rect", xmin = .2, xmax = 2.8, ymin = .25, ymax = 1.25,
                      fill = "#f2f8fb", colour = "#c8d7df", linewidth = .5) +
    ggplot2::annotate("rect", xmin = 5.2, xmax = 7.8, ymin = .25, ymax = 1.25,
                      fill = "#f2f8fb", colour = "#c8d7df", linewidth = .5) +
    ggplot2::annotate(
      "segment", x = 2.9, xend = 5.05, y = 2.75, yend = 2.75,
      arrow = grid::arrow(length = grid::unit(.16, "inches")), colour = "#587181"
    ) +
    ggplot2::annotate(
      "segment", x = 2.9, xend = 5.05, y = .75, yend = .75,
      arrow = grid::arrow(length = grid::unit(.16, "inches")), colour = "#587181"
    ) +
    ggplot2::annotate("text", x = 1.5, y = 2.92,
                      label = format(s$tested_entities, big.mark = ","),
                      size = 7, fontface = "bold", colour = "#0d4258") +
    ggplot2::annotate("text", x = 1.5, y = 2.55,
                      label = "tested target entities", size = 4.2, colour = "#435665") +
    ggplot2::annotate("text", x = 6.5, y = 2.92,
                      label = format(s$submitted_background_names, big.mark = ","),
                      size = 7, fontface = "bold", colour = "#0d4258") +
    ggplot2::annotate("text", x = 6.5, y = 2.55,
                      label = "unique names submitted for comparison", size = 4.0, colour = "#435665") +
    ggplot2::annotate(
      "text", x = 4.0, y = 2.02,
      label = paste0(
        format(s$tested_without_portable_name, big.mark = ","), " without a portable name\n",
        format(s$tested_duplicate_rows_collapsed, big.mark = ","), " duplicate-name rows collapsed"
      ),
      size = 3.6, colour = "#6a7680"
    ) +
    ggplot2::annotate("text", x = 1.5, y = .92,
                      label = format(s$significant_entities, big.mark = ","),
                      size = 7, fontface = "bold", colour = "#0d4258") +
    ggplot2::annotate("text", x = 1.5, y = .55,
                      label = "significant target entities", size = 4.2, colour = "#435665") +
    ggplot2::annotate("text", x = 6.5, y = .92,
                      label = format(s$submitted_foreground_names, big.mark = ","),
                      size = 7, fontface = "bold", colour = "#0d4258") +
    ggplot2::annotate("text", x = 6.5, y = .55,
                      label = "unique foreground names submitted", size = 4.0, colour = "#435665") +
    ggplot2::coord_cartesian(xlim = c(0, 8), ylim = c(0, 3.5), clip = "off") +
    ggplot2::theme_void() +
    ggplot2::labs(title = "How much of your experiment is represented across selected annotation contexts?")

  p
}

#' Plot primary EchoGO evidence-profile counts
#' @keywords internal
plot_echogo_profile_counts <- function(evidence) {
  keep <- c("TARGET_ONLY", "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT")
  if (is.null(evidence) || !is.data.frame(evidence) || !nrow(evidence)) return(NULL)

  d <- evidence |>
    dplyr::filter(.data$evidence_profile %in% keep) |>
    dplyr::count(.data$evidence_profile, name = "n") |>
    tidyr::complete(evidence_profile = keep, fill = list(n = 0L)) |>
    dplyr::mutate(
      evidence_profile = factor(.data$evidence_profile, levels = keep),
      label = unname(.echogo_profile_labels[as.character(.data$evidence_profile)])
    )

  ggplot2::ggplot(d, ggplot2::aes(x = .data$label, y = .data$n, fill = .data$evidence_profile)) +
    ggplot2::geom_col(width = .68) +
    ggplot2::geom_text(ggplot2::aes(label = .data$n), vjust = -.35, fontface = "bold") +
    ggplot2::scale_fill_manual(values = .echogo_profile_colours, guide = "none") +
    ggplot2::labs(x = NULL, y = "Exact GO terms", title = "What is supported by your target, and what is new?") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold")
    )
}

.echogo_context_summary_data <- function(sources, evidence, contexts) {
  if (is.null(sources) || !is.data.frame(sources) || !nrow(sources) ||
      is.null(contexts) || !is.data.frame(contexts) || !nrow(contexts)) {
    return(tibble::tibble())
  }

  custom <- sources |>
    dplyr::filter(
      .data$source_type == "gprofiler",
      .data$background_mode == "custom_experimental_background"
    )
  if (!nrow(custom)) return(tibble::tibble())

  safe_max <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x)) max(x) else NA_real_
  }

  meta <- custom |>
    dplyr::group_by(.data$context_label, .data$context_role) |>
    dplyr::summarise(
      submitted_foreground_n = safe_max(.data$submitted_foreground_n),
      effective_query_n = safe_max(.data$effective_query_n),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      recognition_fraction = dplyr::if_else(
        !is.na(.data$submitted_foreground_n) & .data$submitted_foreground_n > 0,
        .data$effective_query_n / .data$submitted_foreground_n,
        NA_real_
      )
    )

  target_ids <- evidence$term_id[evidence$target_goseq_supported %in% TRUE]
  hypothesis_ids <- evidence$term_id[evidence$evidence_profile == "ALTERNATIVE_CONTEXT"]

  support <- custom |>
    dplyr::filter(.data$source_qualifies %in% TRUE) |>
    dplyr::group_by(.data$context_label, .data$context_role) |>
    dplyr::summarise(
      target_terms_recovered = dplyr::n_distinct(.data$term_id[.data$term_id %in% target_ids]),
      hypothesis_terms_recovered = dplyr::n_distinct(.data$term_id[.data$term_id %in% hypothesis_ids]),
      .groups = "drop"
    )

  order_tbl <- contexts |>
    dplyr::mutate(context_order = dplyr::row_number()) |>
    dplyr::select(.data$context_label, .data$context_role, .data$context_order)

  order_tbl |>
    dplyr::left_join(meta, by = c("context_label", "context_role")) |>
    dplyr::left_join(support, by = c("context_label", "context_role")) |>
    dplyr::mutate(
      target_terms_recovered = dplyr::coalesce(.data$target_terms_recovered, 0L),
      hypothesis_terms_recovered = dplyr::coalesce(.data$hypothesis_terms_recovered, 0L),
      display_label = dplyr::if_else(
        .data$context_role == "TARGET",
        paste0(.data$context_label, "  [target context]"),
        .data$context_label
      )
    ) |>
    dplyr::arrange(.data$context_order)
}

#' Plot how selected annotation contexts represent the same experimental query
#' @keywords internal
plot_echogo_context_summary <- function(sources, evidence, contexts) {
  d <- .echogo_context_summary_data(sources, evidence, contexts)
  if (!nrow(d)) return(NULL)

  d$display_label <- factor(d$display_label, levels = rev(d$display_label))
  fmt_pct <- function(x) ifelse(is.na(x), "NA", paste0(round(100 * x), "%"))
  max_target <- max(c(1, d$target_terms_recovered), na.rm = TRUE)
  max_hyp <- max(c(1, d$hypothesis_terms_recovered), na.rm = TRUE)

  base_theme <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 11.5)
    )

  p1 <- ggplot2::ggplot(d, ggplot2::aes(.data$recognition_fraction, .data$display_label)) +
    ggplot2::geom_col(width = .62, fill = "#176b87", na.rm = TRUE) +
    ggplot2::geom_text(
      ggplot2::aes(label = dplyr::if_else(
        is.na(.data$effective_query_n) | is.na(.data$submitted_foreground_n),
        "not recorded",
        paste0(.data$effective_query_n, "/", .data$submitted_foreground_n)
      )),
      hjust = -.08, size = 3.1, na.rm = TRUE
    ) +
    ggplot2::scale_x_continuous(
      labels = fmt_pct,
      limits = c(0, 1.12),
      breaks = c(0, .25, .5, .75, 1)
    ) +
    ggplot2::labs(x = "Foreground recognized", y = NULL, title = "A  Query recognition") +
    base_theme

  p2 <- ggplot2::ggplot(d, ggplot2::aes(.data$target_terms_recovered, .data$display_label)) +
    ggplot2::geom_col(width = .62, fill = "#2a8f88") +
    ggplot2::geom_text(ggplot2::aes(label = .data$target_terms_recovered), hjust = -.15, size = 3.1) +
    ggplot2::scale_x_continuous(limits = c(0, max_target * 1.18 + 1)) +
    ggplot2::labs(x = "Target findings recovered", y = NULL, title = "B  What the context also recovers") +
    base_theme +
    ggplot2::theme(axis.text.y = ggplot2::element_blank())

  p3 <- ggplot2::ggplot(d, ggplot2::aes(.data$hypothesis_terms_recovered, .data$display_label)) +
    ggplot2::geom_col(width = .62, fill = "#d99a26") +
    ggplot2::geom_text(ggplot2::aes(label = .data$hypothesis_terms_recovered), hjust = -.15, size = 3.1) +
    ggplot2::scale_x_continuous(limits = c(0, max_hyp * 1.18 + 1)) +
    ggplot2::labs(x = "New hypotheses", y = NULL, title = "C  What the context adds") +
    base_theme +
    ggplot2::theme(axis.text.y = ggplot2::element_blank())

  p1 + p2 + p3 +
    patchwork::plot_layout(widths = c(1.45, 1, 1)) +
    patchwork::plot_annotation(
      title = "What do selected annotation contexts recover?",
      subtitle = "The same experimental query is shown in the declared order; contexts are not ranked."
    )
}

.echogo_landscape_terms <- function(evidence, ontology, max_terms_per_profile = 12L) {
  profile_order <- c("TARGET_ONLY", "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT")
  evidence |>
    dplyr::filter(
      .data$ontology %in% c(ontology, paste0("GO:", ontology)),
      .data$evidence_profile %in% profile_order
    ) |>
    dplyr::mutate(profile_order = match(.data$evidence_profile, profile_order)) |>
    dplyr::arrange(.data$profile_order, .data$display_order, .data$term_id) |>
    dplyr::group_by(.data$evidence_profile) |>
    dplyr::slice_head(n = max_terms_per_profile) |>
    dplyr::ungroup()
}

.echogo_landscape_long <- function(evidence, sources, contexts, ontology,
                                   max_terms_per_profile = 12L) {
  terms <- .echogo_landscape_terms(evidence, ontology, max_terms_per_profile)
  if (!nrow(terms)) {
    return(list(terms = terms, matrix = tibble::tibble(), columns = tibble::tibble()))
  }

  ctx <- contexts |>
    dplyr::mutate(context_order = dplyr::row_number()) |>
    dplyr::arrange(.data$context_order)

  cols <- tibble::tibble(
    source_column = c("Target GOseq", ctx$context_label),
    source_order = seq_len(nrow(ctx) + 1L),
    source_role = c("TARGET_GOSEQ", ctx$context_role)
  )

  grid <- tidyr::crossing(term_id = terms$term_id, source_column = cols$source_column) |>
    dplyr::left_join(cols, by = "source_column")

  gs <- terms |>
    dplyr::transmute(
      term_id = .data$term_id,
      source_column = "Target GOseq",
      qualifies = .data$target_goseq_supported
    )

  gp <- sources |>
    dplyr::filter(
      .data$source_type == "gprofiler",
      .data$background_mode == "custom_experimental_background",
      .data$context_label %in% ctx$context_label,
      .data$term_id %in% terms$term_id
    ) |>
    dplyr::group_by(.data$term_id, source_column = .data$context_label) |>
    dplyr::summarise(qualifies = any(.data$source_qualifies %in% TRUE), .groups = "drop")

  matrix <- grid |>
    dplyr::left_join(dplyr::bind_rows(gs, gp), by = c("term_id", "source_column")) |>
    dplyr::mutate(qualifies = dplyr::coalesce(.data$qualifies, FALSE)) |>
    dplyr::left_join(
      terms |>
        dplyr::select(
          .data$term_id, .data$term_name, .data$evidence_profile,
          .data$alternative_context_support_n, .data$alternative_queried_context_n,
          .data$alternative_context_support_fraction, .data$display_order
        ),
      by = "term_id"
    )

  list(terms = terms, matrix = matrix, columns = cols)
}

#' Plot the scoreless EchoGO Evidence Landscape for one ontology
#' @keywords internal
plot_echogo_evidence_landscape <- function(evidence, sources, contexts, ontology,
                                           max_terms_per_profile = 12L,
                                           wrap_width = 44L) {
  dat <- .echogo_landscape_long(
    evidence, sources, contexts, ontology,
    max_terms_per_profile = max_terms_per_profile
  )
  d <- dat$matrix
  if (!nrow(d)) return(NULL)

  profile_order <- c("TARGET_ONLY", "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT")
  profile_labels <- unname(.echogo_profile_labels[profile_order])

  d <- d |>
    dplyr::mutate(
      profile_label = factor(
        unname(.echogo_profile_labels[.data$evidence_profile]),
        levels = profile_labels
      ),
      source_column = factor(.data$source_column, levels = dat$columns$source_column),
      source_role_plot = dplyr::case_when(
        .data$source_role == "TARGET_GOSEQ" ~ "Target GOseq",
        .data$source_role == "TARGET" ~ "Target context",
        TRUE ~ "Alternative context"
      ),
      term_label = paste0(stringr::str_wrap(.data$term_name, width = wrap_width), "\n", .data$term_id)
    )

  term_order <- d |>
    dplyr::distinct(.data$term_id, .data$term_label, .data$evidence_profile, .data$display_order) |>
    dplyr::mutate(profile_order = match(.data$evidence_profile, profile_order)) |>
    dplyr::arrange(.data$profile_order, .data$display_order, .data$term_id) |>
    dplyr::pull(.data$term_label)

  d$term_label <- factor(d$term_label, levels = rev(unique(term_order)))

  support_colours <- c(
    "Target GOseq" = "#176b87",
    "Target context" = "#6f3c9f",
    "Alternative context" = "#d99a26"
  )

  p_matrix <- ggplot2::ggplot(d, ggplot2::aes(.data$source_column, .data$term_label)) +
    ggplot2::geom_point(
      shape = 21, size = 3.0, stroke = .55,
      fill = "white", colour = "#d3dde2"
    ) +
    ggplot2::geom_point(
      data = d[d$qualifies %in% TRUE, , drop = FALSE],
      ggplot2::aes(fill = .data$source_role_plot),
      shape = 21, size = 3.3, stroke = .65, colour = "#314854"
    ) +
    ggplot2::scale_fill_manual(values = support_colours, name = "Qualifying source evidence") +
    ggplot2::facet_grid(profile_label ~ ., scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = paste0("EchoGO Evidence Landscape - ", ontology),
      subtitle = "Filled cells show qualifying source-local evidence. No cross-context significance score is calculated."
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", size = 9),
      legend.position = "bottom"
    )

  recurrence <- d |>
    dplyr::distinct(
      .data$term_id, .data$term_label, .data$evidence_profile,
      .data$alternative_context_support_n, .data$alternative_queried_context_n,
      .data$alternative_context_support_fraction
    ) |>
    dplyr::mutate(
      profile_label = factor(
        unname(.echogo_profile_labels[.data$evidence_profile]),
        levels = profile_labels
      ),
      term_label = factor(.data$term_label, levels = levels(d$term_label)),
      recurrence_label = paste0(
        .data$alternative_context_support_n, "/", .data$alternative_queried_context_n
      )
    )

  p_rec <- ggplot2::ggplot(
    recurrence,
    ggplot2::aes(
      .data$alternative_context_support_fraction, .data$term_label,
      fill = .data$evidence_profile
    )
  ) +
    ggplot2::geom_col(width = .55) +
    ggplot2::geom_text(ggplot2::aes(label = .data$recurrence_label), hjust = -.12, size = 3.0) +
    ggplot2::scale_fill_manual(values = .echogo_profile_colours, guide = "none") +
    ggplot2::scale_x_continuous(
      limits = c(0, 1.20), breaks = c(0, .5, 1), labels = c("0", "50%", "100%")
    ) +
    ggplot2::facet_grid(profile_label ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(x = "Annotation contexts with qualifying evidence", y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
    )

  p_matrix + p_rec + patchwork::plot_layout(widths = c(4.2, 1.25))
}

#' Plot a universal provenance trace for one hypothesis
#' @keywords internal
plot_echogo_hypothesis_trace <- function(drilldown, max_genes = 10L) {
  if (is.null(drilldown) || is.null(drilldown$term) || !nrow(drilldown$term)) return(NULL)

  term <- drilldown$term
  src <- drilldown$sources
  is_hypothesis <- identical(as.character(term$evidence_profile[[1]]), "ALTERNATIVE_CONTEXT")

  ctx <- unique(src$context_label[
    src$source_type == "gprofiler" &
      src$background_mode == "custom_experimental_background" &
      src$context_role == "ALTERNATIVE" &
      src$source_qualifies %in% TRUE
  ])
  ctx <- ctx[!is.na(ctx) & nzchar(ctx)]

  genes <- unique(trimws(unlist(strsplit(
    paste(stats::na.omit(src$contributing_genes), collapse = ";"), "[,;]"
  ))))
  genes <- genes[nzchar(genes)]
  if (length(genes) > max_genes) {
    genes <- c(genes[seq_len(max_genes)], paste0("+", length(genes) - max_genes, " more"))
  }

  ctx_text <- if (length(ctx)) paste(ctx, collapse = " | ") else "No qualifying alternative annotation context"
  gene_text <- if (length(genes)) paste(genes, collapse = " | ") else "No contributing-gene labels available"
  target_text <- if (isTRUE(term$target_goseq_supported[[1]])) {
    "Qualifying target GOseq support"
  } else {
    "No qualifying target GOseq support"
  }
  header_text <- if (is_hypothesis) {
    "Why did this new hypothesis appear?"
  } else {
    "Why did this GO term appear in the report?"
  }

  ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin=.3, xmax=7.7, ymin=5.2, ymax=6.35,
                      fill="#f2f8fb", colour="#cbdbe3") +
    ggplot2::annotate("text", x=.6, y=5.98, hjust=0,
                      label=term$term_name[[1]], size=5.4, fontface="bold", colour="#0d4258") +
    ggplot2::annotate("text", x=.6, y=5.53, hjust=0,
                      label=paste(term$term_id[[1]], term$ontology[[1]], sep="  |  "),
                      size=3.8, colour="#526574") +
    ggplot2::annotate("segment", x=4, xend=4, y=5.15, yend=4.65,
                      arrow=grid::arrow(length=grid::unit(.15,"inches")), colour="#7d909c") +
    ggplot2::annotate("rect", xmin=.6, xmax=7.4, ymin=3.55, ymax=4.6,
                      fill="#fff8eb", colour="#e4c98f") +
    ggplot2::annotate("text", x=.9, y=4.26, hjust=0,
                      label="Your target analysis", size=3.9, fontface="bold", colour="#0d4258") +
    ggplot2::annotate("text", x=.9, y=3.86, hjust=0,
                      label=target_text, size=3.7, colour="#526574") +
    ggplot2::annotate("segment", x=4, xend=4, y=3.5, yend=3.0,
                      arrow=grid::arrow(length=grid::unit(.15,"inches")), colour="#7d909c") +
    ggplot2::annotate("rect", xmin=.6, xmax=7.4, ymin=1.85, ymax=2.95,
                      fill="#eef7f8", colour="#bcd9d6") +
    ggplot2::annotate("text", x=.9, y=2.64, hjust=0,
                      label=paste0(
                        "Recovered in ", term$alternative_context_support_n[[1]], "/",
                        term$alternative_queried_context_n[[1]], " selected alternative annotation contexts"
                      ), size=3.9, fontface="bold", colour="#0d4258") +
    ggplot2::annotate("text", x=.9, y=2.18, hjust=0,
                      label=stringr::str_wrap(ctx_text, 90), size=3.45, colour="#526574") +
    ggplot2::annotate("segment", x=4, xend=4, y=1.8, yend=1.3,
                      arrow=grid::arrow(length=grid::unit(.15,"inches")), colour="#7d909c") +
    ggplot2::annotate("rect", xmin=.6, xmax=7.4, ymin=.2, ymax=1.25,
                      fill="#f7f3fb", colour="#d8c8e5") +
    ggplot2::annotate("text", x=.9, y=.95, hjust=0,
                      label="Contributing genes reported by qualifying sources",
                      size=3.8, fontface="bold", colour="#0d4258") +
    ggplot2::annotate("text", x=.9, y=.52, hjust=0,
                      label=stringr::str_wrap(gene_text, 90), size=3.4, colour="#526574") +
    ggplot2::coord_cartesian(xlim=c(0,8), ylim=c(0,6.5), clip="off") +
    ggplot2::theme_void() +
    ggplot2::labs(title=header_text)
}

#' Save report-facing plot variants
#' @keywords internal
.echogo_save_plot_variants <- function(plot, directory, stem,
                                       width = 12, height = 8, dpi = 320) {
  if (is.null(plot)) return(list(png = NA_character_, pdf = NA_character_, svg = NA_character_))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  png <- file.path(directory, paste0(stem, ".png"))
  pdf <- file.path(directory, paste0(stem, ".pdf"))
  svg <- file.path(directory, paste0(stem, ".svg"))

  ggplot2::ggsave(png, plot, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(pdf, plot, width = width, height = height, bg = "white", limitsize = FALSE)

  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(
      svg, plot, width = width, height = height,
      device = svglite::svglite, bg = "white", limitsize = FALSE
    )
  } else {
    svg <- NA_character_
  }

  list(png = png, pdf = pdf, svg = svg)
}

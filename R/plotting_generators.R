#' Plot generation orchestrators (internal)
#'
#' Makers that write the standard EchoGO lollipop figures to disk. These are not exported;
#' they are called by the pipeline to populate the canonical layout:
#' - goseq/                                    (GOseq)
#' - gprofiler/{with_custom_background,no_background_genome_wide}  (g:Profiler)
#' - consensus/{plots_strict,plots_exploratory} (Consensus lollipops)
#'
#' Functions in this file:
#' - goseq_make_lollipops()
#' - gprofiler_make_lollipops()
#' - consensus_make_lollipops()
#'
#' @name echogo_plot_generators
#' @keywords internal
#' @noRd
#' @importFrom readr read_csv
#' @importFrom dplyr mutate filter arrange slice_head across case_when coalesce select distinct
#' @importFrom ggplot2 ggplot aes geom_segment geom_point coord_flip labs theme_minimal ggsave scale_color_gradient expand_limits
NULL

`%||%` <- function(a, b) if (!is.null(a)) a else b

.make_dir <- function(...) {
  p <- file.path(...); dir.create(p, recursive = TRUE, showWarnings = FALSE); p
}

# ------------------------------------------------------------------------------
# Shared helpers (internal)
# ------------------------------------------------------------------------------

# Ontology palette (matches Rmd)
.echogo_onto_col <- function(ont) {
  c(BP = "firebrick", MF = "forestgreen", CC = "dodgerblue", KEGG = "darkorange")[ont] %||% "purple"
}
.echogo_low_col <- function() "yellow"

# Choose a human-readable term column
.pick_label <- function(d) {
  for (nm in c("term_name","term","label","term_id")) if (nm %in% names(d)) return(as.character(d[[nm]]))
  rep(NA_character_, nrow(d))
}

# Internal: draw + save a lollipop (no aes_string; no deprecation warnings)
.plot_lollipop_core <- function(df,
                                y_col,
                                title, outfile,
                                x_lab = "Term", y_lab = NULL,
                                colour_col = NULL, size_col = NULL,
                                high_colour = NULL,
                                label_width = 35,
                                width = 12, height = 9,
                                # --- consensus-only extras ---
                                subtitle_text = NULL,
                                legend_color_title = NULL,
                                legend_size_title  = NULL,
                                low_colour_override = NULL) {

  stopifnot(is.data.frame(df), nrow(df) > 0, y_col %in% names(df))
  y_lab <- y_lab %||% y_col

  df$.__label__ <- .pick_label(df)
  if (all(is.na(df$.__label__))) stop(".plot_lollipop_core(): no term labels found.")

  trunc_fun <- function(x, w) ifelse(nchar(x) > w, paste0(substr(x, 1, max(1, w - 3)), "..."), x)
  df$.__short__ <- trunc_fun(df$.__label__, label_width)
  df$.__short__ <- make.unique(as.character(df$.__short__))
  df$.__short__ <- factor(df$.__short__, levels = rev(df$.__short__))

  # Create fixed aux columns so we can use standard aes()
  df$.__y__ <- df[[y_col]]
  has_col   <- !is.null(colour_col) && (colour_col %in% names(df))
  has_size  <- !is.null(size_col)   && (size_col   %in% names(df))
  if (has_col)  df$.__col__  <- df[[colour_col]]
  if (has_size) df$.__size__ <- df[[size_col]]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$.__short__, y = .data$.__y__)) +
    ggplot2::geom_segment(
      mapping = if (has_col)
        ggplot2::aes(xend = .data$.__short__, y = 0, yend = .data$.__y__, colour = .data$.__col__)
      else
        ggplot2::aes(xend = .data$.__short__, y = 0, yend = .data$.__y__),
      linewidth = 1.2, lineend = "round"
    ) +
    ggplot2::geom_point(
      mapping = if (has_col && has_size)
        ggplot2::aes(colour = .data$.__col__, size = .data$.__size__)
      else if (has_col)
        ggplot2::aes(colour = .data$.__col__)
      else if (has_size)
        ggplot2::aes(size = .data$.__size__)
      else
        ggplot2::aes(),
      shape = 16, stroke = 0
    ) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle_text,
      x = x_lab, y = y_lab,
      colour = if (has_col) (legend_color_title %||% colour_col) else NULL,
      size   = if (has_size) (legend_size_title  %||% size_col)   else NULL
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::expand_limits(y = 0)

  if (has_col && !is.null(high_colour)) {
    p <- p + ggplot2::scale_color_gradient(
      low  = low_colour_override %||% .echogo_low_col(),
      high = high_colour
    )
  }

  ggplot2::ggsave(outfile, p, width = width, height = height)
  invisible(p)
}

# ---- helpers (safe access) ----------------------------------------------------

# Safe coalesce over maybe-missing columns (vectorized)
.coalesce_cols <- function(d, cols, fallback) {
  present <- intersect(cols, names(d))
  if (!length(present)) return(rep(fallback, nrow(d)))
  out <- d[[present[1]]]
  if (length(present) > 1) {
    for (nm in present[-1]) out <- dplyr::coalesce(out, d[[nm]])
  }
  dplyr::coalesce(out, fallback)
}

# Flexible TRUE-consensus detector (unused here, retained for compatibility)
.detect_true_consensus <- function(d) {
  n <- nrow(d)
  if (!n) return(logical(0))

  has_bool <- all(c("in_goseq","num_species_gprof_bg") %in% names(d))
  if (has_bool) {
    return( (as.logical(d$in_goseq) %in% TRUE) &
              (suppressWarnings(as.numeric(d$num_species_gprof_bg)) > 0) )
  }

  if ("origin" %in% names(d)) {
    return(grepl("Consensus.*with.*BG|with.*background", d$origin, ignore.case = TRUE))
  }

  has_any_score <- ("consensus_score" %in% names(d)) & !is.na(d$consensus_score)
  has_gs  <- "fold_enrichment_goseq" %in% names(d)     & !is.na(d$fold_enrichment_goseq)
  has_bg  <- "avg_fold_gprof_bg"     %in% names(d)     & !is.na(d$avg_fold_gprof_bg)
  return( (has_any_score & (has_gs | has_bg)) )
}

# ------------------------------------------------------------------------------
#  A) GOseq lollipops (BP/CC/MF)
# ------------------------------------------------------------------------------

#' @rdname echogo_plot_generators
#' @keywords internal
#' @noRd
#' @param goseq_csv Path to the annotated GOseq CSV (e.g., GOseq_enrichment_full_annotated.csv).
#' @param outdir_goseq Output directory (usually <base>/goseq).
goseq_make_lollipops <- function(goseq_csv, outdir_goseq) {
  if (!file.exists(goseq_csv)) {
    warning("goseq_make_lollipops(): missing ", goseq_csv)
    return(invisible(NULL))
  }
  df <- suppressMessages(readr::read_csv(goseq_csv, show_col_types = FALSE))
  req_any <- c("foldEnrichment","term","ontology","over_represented_FDR")
  if (!all(req_any %in% names(df))) {
    warning("goseq_make_lollipops(): CSV missing columns: ",
            paste(setdiff(req_any, names(df)), collapse = ", "))
    return(invisible(NULL))
  }

  outdir <- .make_dir(outdir_goseq)  # typically <base>/goseq

  df <- dplyr::mutate(
    df,
    term_name   = .data$term,
    negLogFDR   = -log10(pmax(.data$over_represented_FDR, .Machine$double.xmin)),
    ontology    = dplyr::case_when(
      .data$ontology %in% c("BP","GO:BP") ~ "BP",
      .data$ontology %in% c("MF","GO:MF") ~ "MF",
      .data$ontology %in% c("CC","GO:CC") ~ "CC",
      TRUE ~ as.character(.data$ontology)
    )
  )

  ontos <- c("BP","CC","MF")

  # ---- Per-depth top 50 (only if depth present)
  if ("depth" %in% names(df)) {
    depths <- sort(unique(stats::na.omit(df$depth)))
    for (depth_val in depths) {
      for (ont in ontos) {
        d <- df |>
          dplyr::filter(.data$depth == depth_val, .data$ontology == ont, !is.na(.data$foldEnrichment)) |>
          dplyr::arrange(dplyr::desc(.data$foldEnrichment)) |>
          dplyr::slice_head(n = 50)
        if (!nrow(d)) next
        .plot_lollipop_core(
          df = d,
          y_col = "foldEnrichment",
          title = sprintf("Top Enriched GO Terms — Depth %s — %s", depth_val, ont),
          outfile = file.path(outdir, sprintf("GO_lollipop_depth%s_%s.pdf", depth_val, ont)),
          x_lab = "GO Term",
          y_lab = "Fold Enrichment",
          colour_col = "negLogFDR",
          size_col   = if ("numDEInCat" %in% names(d)) "numDEInCat" else NULL,
          high_colour = .echogo_onto_col(ont),
          label_width = 35
        )
      }
    }
  }

  # ---- All depths top 50 per ontology
  for (ont in ontos) {
    d <- df |>
      dplyr::filter(.data$ontology == ont, !is.na(.data$foldEnrichment)) |>
      dplyr::arrange(dplyr::desc(.data$foldEnrichment)) |>
      dplyr::slice_head(n = 50)
    if (!nrow(d)) next
    .plot_lollipop_core(
      df = d,
      y_col = "foldEnrichment",
      title = sprintf("Top GO Terms (All Depths) — %s", ont),
      outfile = file.path(outdir, sprintf("GO_lollipop_allDepths_top50_%s.pdf", ont)),
      x_lab = "GO Term",
      y_lab = "Fold Enrichment",
      colour_col = "negLogFDR",
      size_col   = if ("numDEInCat" %in% names(d)) "numDEInCat" else NULL,
      high_colour = .echogo_onto_col(ont),
      label_width = 35
    )
  }
  invisible(NULL)
}

# ------------------------------------------------------------------------------
#  B) g:Profiler lollipops (per species; with background & no background)
# ------------------------------------------------------------------------------

# helper (internal)
.gp_find_species_csv <- function(base_dir, mode, sp_label) {
  # Either written as gprofiler_<label>_<mode>.csv or ..._enrichment.csv
  subdir <- if (mode == "with_bg") "with_custom_background" else "no_background_genome_wide"
  root   <- file.path(base_dir, subdir)
  cands  <- c(
    file.path(root, sprintf("gprofiler_%s_%s.csv",           sp_label, if (mode=="with_bg") "with_bg" else "nobg")),
    file.path(root, sprintf("gprofiler_%s_%s_enrichment.csv", sp_label, if (mode=="with_bg") "with_bg" else "nobg"))
  )
  hit <- cands[file.exists(cands)][1]
  if (length(hit)) hit else NA_character_
}

#' @rdname echogo_plot_generators
#' @keywords internal
#' @noRd
#' @param gprof_base_dir Base directory containing g:Profiler outputs
#'   (expects subfolders: with_custom_background/ and no_background_genome_wide/).
#' @param species_map Named character vector mapping species codes to labels,
#'   e.g. c(hsapiens="human", mmusculus="mouse", drerio="zebrafish").
gprofiler_make_lollipops <- function(gprof_base_dir, species_map) {
  ont_tags <- c("GO:BP","GO:MF","GO:CC")

  for (mode in c("with_bg","nobg")) {
    subdir <- if (mode == "with_bg") "with_custom_background" else "no_background_genome_wide"
    outdir <- .make_dir(gprof_base_dir, subdir)

    for (sp_code in names(species_map)) {
      sp_label <- species_map[[sp_code]]
      csv <- .gp_find_species_csv(gprof_base_dir, mode, sp_label)
      if (is.na(csv)) next

      df <- suppressMessages(readr::read_csv(csv, show_col_types = FALSE))
      req <- c("term_id","term_name","source","fold_enrichment","p_value","intersection_size")
      if (!all(req %in% names(df))) next

      # mirror Rmd: plot only significant terms if present; otherwise create a note file
      sig_tbl <- df |> dplyr::filter(!is.na(.data$p_value) & .data$p_value <= 0.05)
      if (!nrow(sig_tbl)) {
        file.create(file.path(outdir, paste0("gprofiler_", sp_label, "_NO_SIGNIFICANT_RESULTS.txt")))
        next
      }

      dff <- dplyr::mutate(sig_tbl, negLogP = -log10(pmax(.data$p_value, .Machine$double.xmin)))

      for (ont in ont_tags) {
        ont_short <- sub("^GO:", "", ont)
        col_hi    <- .echogo_onto_col(ont_short)

        dd <- dff |>
          dplyr::filter(.data$source %in% c(ont, ont_short)) |>
          dplyr::arrange(dplyr::desc(.data$fold_enrichment)) |>
          dplyr::slice_head(n = 50)
        if (!nrow(dd)) next

        .plot_lollipop_core(
          df = dd,
          y_col = "fold_enrichment",
          title = sprintf("g:Profiler (%s) — %s Top 50 Enriched Terms", tools::toTitleCase(sp_label), ont),
          outfile = file.path(outdir, sprintf("gprofiler_%s_%s_%s_lollipop.pdf",
                                              sp_label, if (mode=="with_bg") "with_bg" else "nobg",
                                              gsub(":", "_", ont))),
          x_lab = "Term",
          y_lab = "Fold Enrichment",
          colour_col = "negLogP",
          size_col   = "intersection_size",
          high_colour = col_hi,
          label_width = 35
        )
      }
    }
  }
  invisible(NULL)
}

# ------------------------------------------------------------------------------
#  C) Consensus lollipops — True Consensus vs Exploratory
# ------------------------------------------------------------------------------
#' Plot EchoGO consensus lollipops (True Consensus & Exploratory) to mirror Rmd
#' @rdname echogo_plot_generators
#' @keywords internal
#' @noRd
#' @param consensus_df data.frame produced by build_consensus_table()
#' @param base_outdir canonical base dir for consensus outputs (expects subfolders we create)
#' @param top_n integer, default 50
#' @param legacy_root optional; if provided, mirrors plots into legacy folders
#'        (`consensus_plots_strict_true_consensus/` and `consensus_plots_all_exploratory/`)
consensus_make_lollipops <- function(consensus_df, base_outdir, top_n = 50, legacy_root = NULL) {
  if (!is.data.frame(consensus_df) || !nrow(consensus_df)) {
    warning("consensus_make_lollipops(): empty consensus_df")
    return(invisible(NULL))
  }

  # -- Precompute a robust label (no cur_data() inside mutate)
  cdf <- consensus_df
  cdf$term_name <- .coalesce_cols(
    cdf,
    c("term_name", "term_name.bg", "term_name.nobg", "term_id"),
    NA_character_
  )

  # -- Now do the rest safely in mutate
  cdf <- dplyr::mutate(
    cdf,
    ontology = dplyr::case_when(
      .data$ontology %in% c("GO:BP","BP") ~ "BP",
      .data$ontology %in% c("GO:CC","CC") ~ "CC",
      .data$ontology %in% c("GO:MF","MF") ~ "MF",
      grepl("^KEGG", .data$ontology %||% "", ignore.case = TRUE) ~ "KEGG",
      TRUE ~ as.character(.data$ontology)
    ),
    # color by best available fold enrichment
    fold_color = dplyr::coalesce(.data$fold_enrichment_goseq,
                                 .data$avg_fold_gprof_bg,
                                 .data$avg_fold_gprof_nobg),
    # size by strongest (minimum) p across methods
    .pmin_any = pmin(
      dplyr::coalesce(.data$min_pval_goseq,      Inf),
      dplyr::coalesce(.data$min_pval_gprof_bg,   Inf),
      dplyr::coalesce(.data$min_pval_gprof_nobg, Inf),
      na.rm = TRUE
    ),
    negLogP = -log10(pmax(.data$.pmin_any, .Machine$double.xmin))
  )

  if (".pmin_any" %in% names(cdf)) cdf <- dplyr::select(cdf, -".pmin_any")

  if ("significant_in_any" %in% names(cdf)) {
    cdf <- dplyr::filter(cdf, .data$significant_in_any == TRUE)
  }

  # Canonical output roots (changed!)
  out_true <- .make_dir(base_outdir, "plots_strict")
  out_expl <- .make_dir(base_outdir, "plots_exploratory")

  # ---------- TRUE CONSENSUS ----------
  true_df <- dplyr::filter(
    cdf,
    .data$ontology %in% c("BP","MF","CC","KEGG"),
    if ("origin" %in% names(cdf)) .data$origin %in% "GO terms - Consensus (with BG)" else TRUE
  )

  ontos_no_kegg <- c("BP","MF","CC")
  depths_true <- if ("depth" %in% names(true_df)) sort(unique(stats::na.omit(true_df$depth))) else numeric(0)

  for (ont in ontos_no_kegg) {
    hi <- .echogo_onto_col(ont)

    if (length(depths_true)) {
      for (depth_val in depths_true) {
        dd <- true_df |>
          dplyr::filter(.data$ontology == ont, .data$depth == depth_val) |>
          dplyr::arrange(dplyr::desc(.data$consensus_score)) |>
          dplyr::slice_head(n = top_n)
        if (!nrow(dd)) next
        .plot_lollipop_core(
          df = dd,
          y_col = "consensus_score",
          title = sprintf("Top %d Consensus Terms - %s (Depth %s)", top_n, ont, depth_val),
          outfile = file.path(out_true, sprintf("lollipop_top%d_%s_depth%s.pdf", top_n, ont, depth_val)),
          x_lab = "Term",
          y_lab = "Consensus Score",
          colour_col = "fold_color",
          size_col   = "negLogP",
          high_colour = hi,
          label_width = 35,
          # Rmd-aligned legends/subtitle (CONSENSUS ONLY)
          subtitle_text       = "Sorted by consensus score; color = fold enrichment",
          legend_color_title  = "Fold Enrichment",
          legend_size_title   = expression(-log[10](p-value)),
          low_colour_override = "gray"
        )
      }
    }

    dd_all <- true_df |>
      dplyr::filter(.data$ontology == ont) |>
      dplyr::arrange(dplyr::desc(.data$consensus_score)) |>
      dplyr::slice_head(n = top_n)
    if (nrow(dd_all)) {
      .plot_lollipop_core(
        df = dd_all,
        y_col = "consensus_score",
        title = sprintf("Top %d Consensus Terms (All Depths) - %s", top_n, ont),
        outfile = file.path(out_true, sprintf("lollipop_top%d_%s_alldepths.pdf", top_n, ont)),
        x_lab = "Term",
        y_lab = "Consensus Score",
        colour_col = "fold_color",
        size_col   = "negLogP",
        high_colour = hi,
        label_width = 35,
        subtitle_text       = "Sorted by consensus score; color = fold enrichment",
        legend_color_title  = "Fold Enrichment",
        legend_size_title   = expression(-log[10](p-value)),
        low_colour_override = "gray"
      )
    }
  }

  # ---------- EXPLORATORY ----------
  ex_df <- cdf
  if ("consensus_score_all" %in% names(ex_df)) ex_df$consensus_score <- ex_df$consensus_score_all
  depths_ex <- if ("depth" %in% names(ex_df)) sort(unique(stats::na.omit(ex_df$depth))) else numeric(0)

  for (ont in ontos_no_kegg) {
    hi <- .echogo_onto_col(ont)

    if (length(depths_ex)) {
      for (depth_val in depths_ex) {
        dd <- ex_df |>
          dplyr::filter(.data$ontology == ont, .data$depth == depth_val) |>
          dplyr::arrange(dplyr::desc(.data$consensus_score)) |>
          dplyr::slice_head(n = top_n)
        if (!nrow(dd)) next
        .plot_lollipop_core(
          df = dd,
          y_col = "consensus_score",
          title = sprintf("Top %d Exploratory Terms - %s (Depth %s)", top_n, ont, depth_val),
          outfile = file.path(out_expl, sprintf("lollipop_top%d_%s_depth%s.pdf", top_n, ont, depth_val)),
          x_lab = "Term",
          y_lab = "Consensus Score",
          colour_col = "fold_color",
          size_col   = "negLogP",
          high_colour = hi,
          label_width = 35,
          subtitle_text       = "Sorted by consensus score; color = fold enrichment",
          legend_color_title  = "Fold Enrichment",
          legend_size_title   = expression(-log[10](p-value)),
          low_colour_override = "gray"
        )
      }
    }

    dd_all <- ex_df |>
      dplyr::filter(.data$ontology == ont) |>
      dplyr::arrange(dplyr::desc(.data$consensus_score)) |>
      dplyr::slice_head(n = top_n)
    if (nrow(dd_all)) {
      .plot_lollipop_core(
        df = dd_all,
        y_col = "consensus_score",
        title = sprintf("Top %d Exploratory Terms (All Depths) - %s", top_n, ont),
        outfile = file.path(out_expl, sprintf("lollipop_top%d_%s_alldepths.pdf", top_n, ont)),
        x_lab = "Term",
        y_lab = "Consensus Score",
        colour_col = "fold_color",
        size_col   = "negLogP",
        high_colour = hi,
        label_width = 35,
        subtitle_text       = "Sorted by consensus score; color = fold enrichment",
        legend_color_title  = "Fold Enrichment",
        legend_size_title   = expression(-log[10](p-value)),
        low_colour_override = "gray"
      )
    }
  }

  # ---------- KEGG: single "Significant Only" plot ----------
  kegg_dd <- cdf |>
    dplyr::filter(.data$ontology == "KEGG") |>
    dplyr::mutate(
      fold_color = dplyr::coalesce(.data$avg_fold_gprof_bg, .data$avg_fold_gprof_nobg),
      .p_kegg = dplyr::coalesce(.data$min_pval_gprof_bg, .data$min_pval_gprof_nobg, Inf),
      negLogP = -log10(pmax(.data$.p_kegg, .Machine$double.xmin))
    ) |>
    dplyr::arrange(dplyr::desc(.data$consensus_score)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::select(-.p_kegg)

  if (nrow(kegg_dd)) {
    .plot_lollipop_core(
      df = kegg_dd,
      y_col = "consensus_score",
      title = sprintf("Consensus KEGG Terms - Top %d (Significant Only)", top_n),
      outfile = file.path(out_expl, sprintf("lollipop_top%d_KEGG_alldepths.pdf", top_n)),
      x_lab = "KEGG Term",
      y_lab = "Consensus Score",
      colour_col = "fold_color",
      size_col   = "negLogP",
      high_colour = .echogo_onto_col("KEGG"),
      label_width = 35,
      subtitle_text       = "Sorted by consensus score; color = fold enrichment",
      legend_color_title  = "Fold Enrichment",
      legend_size_title   = expression(-log[10](p-value)),
      low_colour_override = "gray"
    )
  }

  invisible(NULL)
}


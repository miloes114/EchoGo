# RRvGO report-output wrapper -------------------------------------------------
# Preserve the existing semantic calculation exactly, then add robust report
# sidecars from the already-reduced cluster tables. No term membership or
# semantic clustering is recalculated here.

.echogo_rrvgo_consensus_analysis_base <- run_rrvgo_consensus_analysis

.echogo_rrvgo_close_new_devices <- function(before_ids) {
  after <- grDevices::dev.list()
  if (is.null(after)) return(invisible(NULL))
  after_ids <- unname(after)
  new_ids <- setdiff(after_ids, before_ids)
  for (id in rev(new_ids)) {
    try(grDevices::dev.off(which = id), silent = TRUE)
  }
  invisible(NULL)
}

.echogo_rrvgo_write_report_sidecars <- function(output_root) {
  if (is.null(output_root) || !length(output_root) || is.na(output_root) || !dir.exists(output_root)) {
    return(invisible(NULL))
  }

  cluster_files <- list.files(
    output_root,
    pattern = "^rrvgo_(BP|MF|CC)_clusters\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (!length(cluster_files)) return(invisible(NULL))

  for (cluster_file in cluster_files) {
    reduced <- tryCatch(
      utils::read.csv(cluster_file, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
    if (is.null(reduced) || !nrow(reduced)) next

    ont <- sub("^rrvgo_([A-Za-z]+)_clusters\\.csv$", "\\1", basename(cluster_file))
    root <- dirname(cluster_file)

    if (all(c("cluster", "score", "size", "term") %in% names(reduced))) {
      d <- reduced[order(reduced$score, decreasing = TRUE), , drop = FALSE]
      d <- utils::head(d, 120L)
      if (!"origin" %in% names(d)) d$origin <- "Evidence"
      p <- ggplot2::ggplot(
        d,
        ggplot2::aes(x = .data$cluster, y = .data$score, size = .data$size, colour = .data$origin)
      ) +
        ggplot2::geom_point(alpha = .72) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::labs(
          title = paste0("RRvGO semantic themes - ", ont),
          x = "Semantic cluster",
          y = "Representative order (descriptive)",
          colour = "Evidence profile",
          size = "Cluster size"
        )
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        lab <- utils::head(d, 24L)
        p <- p + ggrepel::geom_text_repel(
          data = lab,
          ggplot2::aes(label = .data$term),
          size = 3.2,
          max.overlaps = 30,
          show.legend = FALSE
        )
      }
      ggplot2::ggsave(
        file.path(root, paste0("rrvgo_", ont, "_bubbleplot.png")),
        p, width = 12, height = 8, dpi = 320, bg = "white", limitsize = FALSE
      )
    }

    if (requireNamespace("treemap", quietly = TRUE) &&
        all(c("term", "size") %in% names(reduced))) {
      top <- reduced[order(reduced$score, decreasing = TRUE), , drop = FALSE]
      top <- utils::head(top, 160L)
      top <- top[!is.na(top$size) & is.finite(top$size) & top$size > 0, , drop = FALSE]
      if (nrow(top)) {
        if (!"parentTerm" %in% names(top)) top$parentTerm <- "GO themes"
        top$scaled_size <- log1p(top$size)
        png_path <- file.path(root, paste0("rrvgo_", ont, "_treemap.png"))
        grDevices::png(png_path, width = 3000, height = 2200, res = 250, bg = "white")
        ok <- tryCatch({
          treemap::treemap(
            top,
            index = c("parentTerm", "term"),
            vSize = "scaled_size",
            type = "index",
            title = paste0("RRvGO semantic themes - ", ont),
            fontcolor.labels = c("#FFFFFFDD", "#00000080"),
            bg.labels = 0,
            border.col = "#00000080"
          )
          TRUE
        }, error = function(e) FALSE)
        grDevices::dev.off()
        if (!ok && file.exists(png_path)) unlink(png_path, force = TRUE)
      }
    }

    scatter <- file.path(root, paste0("rrvgo_", ont, "_scatterplot.pdf"))
    if (file.exists(scatter)) {
      size <- file.info(scatter)$size
      if (is.na(size) || size < 1000) {
        unlink(scatter, force = TRUE)
        writeLines(
          "RRvGO scatterplot was skipped because the graphics device did not produce a valid PDF. Exact semantic clusters remain available.",
          file.path(root, paste0("rrvgo_", ont, "_scatterplot_SKIPPED.txt"))
        )
      }
    }
  }

  invisible(NULL)
}

run_rrvgo_consensus_analysis <- function(...) {
  before <- grDevices::dev.list()
  before_ids <- if (is.null(before)) integer() else unname(before)
  on.exit(.echogo_rrvgo_close_new_devices(before_ids), add = TRUE)

  out <- .echogo_rrvgo_consensus_analysis_base(...)
  .echogo_rrvgo_close_new_devices(before_ids)
  .echogo_rrvgo_write_report_sidecars(out)
  out
}

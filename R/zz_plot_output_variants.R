# Report-facing plot output variants -----------------------------------------
#
# Preserve the established GOseq/g:Profiler lollipop design while adding PNG
# sidecars for reliable inline HTML display.

.plot_lollipop_core <- function(df,
                                y_col,
                                title, outfile,
                                x_lab = "Term", y_lab = NULL,
                                colour_col = NULL, size_col = NULL,
                                high_colour = NULL,
                                label_width = 35,
                                width = 12, height = 9,
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

  df$.__y__ <- df[[y_col]]
  has_col <- !is.null(colour_col) && (colour_col %in% names(df))
  has_size <- !is.null(size_col) && (size_col %in% names(df))
  if (has_col) df$.__col__ <- df[[colour_col]]
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
      title = title,
      subtitle = subtitle_text,
      x = x_lab,
      y = y_lab,
      colour = if (has_col) (legend_color_title %||% colour_col) else NULL,
      size = if (has_size) (legend_size_title %||% size_col) else NULL
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::expand_limits(y = 0)

  if (has_col && !is.null(high_colour)) {
    p <- p + ggplot2::scale_color_gradient(
      low = low_colour_override %||% .echogo_low_col(),
      high = high_colour
    )
  }

  ext <- tolower(tools::file_ext(outfile))
  stem <- tools::file_path_sans_ext(outfile)
  pdf_path <- if (identical(ext, "pdf")) outfile else paste0(stem, ".pdf")
  png_path <- paste0(stem, ".png")
  svg_path <- paste0(stem, ".svg")

  ggplot2::ggsave(pdf_path, p, width = width, height = height, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(png_path, p, width = width, height = height, dpi = 320, bg = "white", limitsize = FALSE)
  if (requireNamespace("svglite", quietly = TRUE)) {
    try(
      ggplot2::ggsave(
        svg_path, p, width = width, height = height,
        device = svglite::svglite, bg = "white", limitsize = FALSE
      ),
      silent = TRUE
    )
  }

  invisible(p)
}

.echogo_write_empty_pdf <- function(outfile, title, subtitle = NULL) {
  draw <- function() {
    grid::grid.newpage()
    grid::grid.text(title, y = 0.60, gp = grid::gpar(fontsize = 16, fontface = "bold"))
    if (!is.null(subtitle)) {
      grid::grid.text(subtitle, y = 0.48, gp = grid::gpar(fontsize = 12))
    }
  }

  grDevices::pdf(outfile, width = 10, height = 6)
  draw()
  grDevices::dev.off()

  png_path <- paste0(tools::file_path_sans_ext(outfile), ".png")
  grDevices::png(png_path, width = 2400, height = 1440, res = 240, bg = "white")
  draw()
  grDevices::dev.off()

  invisible(c(pdf = outfile, png = png_path))
}

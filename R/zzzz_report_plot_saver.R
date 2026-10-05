# Bounded report figure saver -------------------------------------------------
# Large validation datasets can contain hundreds of exact terms, while the
# report deliberately displays a capped subset. Bound physical device size so
# display-only calculations cannot create impractically large PNG/PDF files.

.echogo_save_plot_variants <- function(plot, directory, stem,
                                       width = 12, height = 8, dpi = 320) {
  if (is.null(plot)) {
    return(list(png = NA_character_, pdf = NA_character_, svg = NA_character_))
  }

  width <- max(6, min(as.numeric(width), 20))
  height <- max(4, min(as.numeric(height), 18))
  dpi <- max(150, min(as.numeric(dpi), 400))

  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  png <- file.path(directory, paste0(stem, ".png"))
  pdf <- file.path(directory, paste0(stem, ".pdf"))
  svg <- file.path(directory, paste0(stem, ".svg"))

  ggplot2::ggsave(
    png, plot, width = width, height = height, dpi = dpi,
    bg = "white", limitsize = FALSE
  )
  ggplot2::ggsave(
    pdf, plot, width = width, height = height,
    bg = "white", limitsize = FALSE
  )

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

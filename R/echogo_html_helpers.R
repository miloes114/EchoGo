#' Create an HTML details wrapper
#' @keywords internal
.echogo_html_details <- function(summary, inner_html) {
  paste0(
    "<details style='margin:8px 0;'>",
    "<summary style='cursor:pointer;font-weight:600;'>", summary, "</summary>",
    "<div style='margin-top:10px;'>", inner_html, "</div>",
    "</details>"
  )
}

.echogo_escape_html <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

#' Convert an output path to a browser-friendly URL
#' If target and report_dir share the same base_dir, return a relative URL.
#' Otherwise return a file:/// URL.
#' @keywords internal
.echogo_path_for_html <- function(p, report_dir, base_dir) {
  if (is.null(p) || is.na(p) || !nzchar(p)) return(NA_character_)

  norm <- function(x) {
    x <- normalizePath(x, winslash = "/", mustWork = FALSE)
    x <- gsub("/+$", "", x)
    x
  }

  p  <- norm(p)
  rd <- norm(report_dir)
  bd <- norm(base_dir)

  p_key  <- if (.Platform$OS.type == "windows") tolower(p)  else p
  rd_key <- if (.Platform$OS.type == "windows") tolower(rd) else rd
  bd_key <- if (.Platform$OS.type == "windows") tolower(bd) else bd

  if (nzchar(bd) &&
      startsWith(p_key,  paste0(bd_key, "/")) &&
      startsWith(rd_key, paste0(bd_key, "/"))) {

    rel_target <- sub(paste0("^", gsub("([\\W])", "\\\\\\1", bd), "/?"), "", p)
    rel_report <- sub(paste0("^", gsub("([\\W])", "\\\\\\1", bd), "/?"), "", rd)

    depth <- if (nzchar(rel_report)) length(strsplit(rel_report, "/", fixed = TRUE)[[1]]) else 0
    up <- if (depth > 0) paste(rep("..", depth), collapse = "/") else "."
    return(gsub("//+", "/", file.path(up, rel_target)))
  }

  paste0("file:///", p)
}

#' Embed a PDF with a collapsible toggle
#' @param title  Title shown in the toggle summary
#' @param path   Absolute path to the PDF
#' @param report_dir Directory where the HTML is written
#' @param base_dir   Root of all outputs (used for relative paths)
#' @param height CSS height (e.g., "650px")
#' @return HTML string (to be cat()'d with results='asis')
#' @export
echogo_embed_pdf_toggle <- function(title, path, report_dir, base_dir, height = "650px") {
  if (is.null(path) || is.na(path) || !nzchar(path)) {
    return(.echogo_html_details(title, "<em>Missing file:</em> (no path)"))
  }
  if (!file.exists(path)) {
    return(.echogo_html_details(title, sprintf("<em>Missing file:</em> %s", .echogo_escape_html(path))))
  }

  url <- .echogo_stage_asset(path, report_dir = report_dir, base_dir = base_dir)
  if (is.na(url)) {
    return(.echogo_html_details(title, sprintf("<em>Failed to stage file:</em> %s", .echogo_escape_html(path))))
  }

  inner <- sprintf("<object data='%s' type='application/pdf' width='100%%' height='%s'>", url, height)
  inner <- paste0(
    inner,
    "<p>PDF preview not supported. <a href='", url, "' target='_blank'>Open</a>.</p>",
    "</object>"
  )

  .echogo_html_details(title, inner)
}

#' Embed an HTML file (iframe) with a collapsible toggle
#' @param title Title shown in the toggle summary
#' @param path_abs Absolute path to the HTML file
#' @param report_dir Directory where the report HTML lives
#' @param base_dir Root of all outputs (used to compute relative URLs)
#' @param height CSS height (e.g., "650px")
#' @return HTML string (to be cat()'d with results='asis')
#' @export
embed_html_toggle_external_plus <- function(title, path_abs, report_dir, base_dir, height = "650px") {
  if (is.null(path_abs) || is.na(path_abs) || !nzchar(path_abs)) {
    return(.echogo_html_details(title, "<em>Missing HTML:</em> (no path)"))
  }
  if (!file.exists(path_abs)) {
    return(.echogo_html_details(title, "<em>Not found in index or not generated for this run.</em>"))
  }

  url <- .echogo_stage_asset(path_abs, report_dir = report_dir, base_dir = base_dir)
  if (is.na(url)) {
    return(.echogo_html_details(title, sprintf("<em>Failed to stage HTML:</em> %s", .echogo_escape_html(path_abs))))
  }

  iframe <- sprintf(
    "<iframe src='%s' data-external='1' style='width:100%%;height:%s;border:1px solid #ddd;border-radius:6px;'></iframe>",
    url, height
  )

  pdf_path <- file.path(
    dirname(path_abs),
    paste0(tools::file_path_sans_ext(basename(path_abs)), ".pdf")
  )
  fallback <- ""
  if (file.exists(pdf_path)) {
    pdf_url <- .echogo_stage_asset(pdf_path, report_dir = report_dir, base_dir = base_dir)
    if (!is.na(pdf_url)) {
      fallback <- paste0(
        "<p><a href='", pdf_url,
        "' target='_blank' rel='noopener'>Open the static PDF version</a></p>"
      )
    }
  }

  .echogo_html_details(title, paste0(iframe, fallback))
}

#' Embed an HTML file (iframe) with a collapsible toggle (deprecated alias)
#'
#' Deprecated. Use [embed_html_toggle_external_plus()] instead.
#' @inheritParams embed_html_toggle_external_plus
#' @param path Absolute path to the HTML file (deprecated name for `path_abs`).
#' @export
echogo_embed_html_toggle <- function(title, path, report_dir, base_dir, height = "650px") {
  .Deprecated("embed_html_toggle_external_plus", package = "EchoGO")
  embed_html_toggle_external_plus(
    title      = title,
    path_abs   = path,
    report_dir = report_dir,
    base_dir   = base_dir,
    height     = height
  )
}

.echogo_stage_html_dependencies <- function(path_abs, dst) {
  if (!identical(tolower(tools::file_ext(path_abs)), "html")) {
    return(invisible(TRUE))
  }

  dep_src <- paste0(tools::file_path_sans_ext(path_abs), "_files")
  if (!dir.exists(dep_src)) return(invisible(TRUE))

  dep_dst <- file.path(dirname(dst), basename(dep_src))
  dir.create(dep_dst, recursive = TRUE, showWarnings = FALSE)

  dep_files <- list.files(
    dep_src,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = FALSE,
    no.. = TRUE
  )
  if (!length(dep_files)) return(invisible(TRUE))

  dep_root <- normalizePath(dep_src, winslash = "/", mustWork = TRUE)
  rel <- substring(
    normalizePath(dep_files, winslash = "/", mustWork = TRUE),
    nchar(dep_root) + 2L
  )
  targets <- file.path(dep_dst, rel)
  invisible(vapply(unique(dirname(targets)), dir.create, logical(1),
                   recursive = TRUE, showWarnings = FALSE))

  copied <- file.copy(dep_files, targets, overwrite = TRUE, copy.date = TRUE)
  if (!all(copied)) {
    warning("Failed to stage one or more HTML widget dependencies for: ", path_abs)
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

#' Stage a report asset and return its relative URL
#' @keywords internal
.echogo_stage_asset <- function(path_abs, report_dir, base_dir) {
  if (is.null(path_abs) || is.na(path_abs) || !nzchar(path_abs)) return(NA_character_)
  path_abs <- normalizePath(path_abs, winslash = "/", mustWork = FALSE)
  if (!file.exists(path_abs)) return(NA_character_)

  assets_dir <- file.path(report_dir, "assets")
  dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)

  norm <- function(x) {
    x <- normalizePath(x, winslash = "/", mustWork = FALSE)
    gsub("/+$", "", x)
  }

  p  <- norm(path_abs)
  bd <- norm(base_dir)

  p_key  <- if (.Platform$OS.type == "windows") tolower(p)  else p
  bd_key <- if (.Platform$OS.type == "windows") tolower(bd) else bd

  prefix <- paste0(bd_key, "/")
  if (nzchar(bd_key) && startsWith(p_key, prefix)) {
    rel <- substr(p, nchar(prefix) + 1, nchar(p))
    dst <- file.path(assets_dir, rel)
    dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    copied <- file.copy(path_abs, dst, overwrite = TRUE, copy.date = TRUE)
    if (!isTRUE(copied)) return(NA_character_)
    if (!isTRUE(.echogo_stage_html_dependencies(path_abs, dst))) return(NA_character_)
    return(gsub("\\\\", "/", file.path("assets", rel)))
  }

  dst <- file.path(assets_dir, "external", basename(path_abs))
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  copied <- file.copy(path_abs, dst, overwrite = TRUE, copy.date = TRUE)
  if (!isTRUE(copied)) return(NA_character_)
  if (!isTRUE(.echogo_stage_html_dependencies(path_abs, dst))) return(NA_character_)
  gsub("\\\\", "/", file.path("assets", "external", basename(path_abs)))
}

# Biological-report helpers --------------------------------------------------

.echogo_report_metric_cards <- function(cards) {
  if (!is.data.frame(cards) || !all(c("value", "label") %in% names(cards))) {
    stop("cards must contain value and label columns", call. = FALSE)
  }
  items <- vapply(seq_len(nrow(cards)), function(i) {
    value <- cards$value[[i]]
    value <- if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
      format(value, big.mark = ",", scientific = FALSE, trim = TRUE)
    } else {
      as.character(value)
    }
    paste0(
      "<div class='metric-card'>",
      "<div class='value'>", .echogo_escape_html(value), "</div>",
      "<div class='label'>", .echogo_escape_html(cards$label[[i]]), "</div>",
      "</div>"
    )
  }, character(1))
  paste0("<div class='metric-grid'>", paste(items, collapse = ""), "</div>")
}

.echogo_report_callout <- function(kind = c("how-to-read", "takeaway", "meaning", "caution", "method-note"),
                                   title, body) {
  kind <- match.arg(kind)
  paste0(
    "<div class='", kind, "'><strong>", .echogo_escape_html(title),
    "</strong>", body, "</div>"
  )
}

.echogo_report_figure_card <- function(title, png, caption,
                                       report_dir, base_dir,
                                       pdf = NULL, svg = NULL,
                                       alt = title) {
  if (is.null(png) || is.na(png) || !nzchar(png) || !file.exists(png)) {
    return(.echogo_html_details(title, "<em>This figure was not generated for this run.</em>"))
  }

  png_url <- .echogo_stage_asset(png, report_dir = report_dir, base_dir = base_dir)
  if (is.na(png_url)) {
    return(.echogo_html_details(title, "<em>The figure was generated but could not be staged into the report.</em>"))
  }

  make_link <- function(path, label) {
    if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) return("")
    u <- .echogo_stage_asset(path, report_dir = report_dir, base_dir = base_dir)
    if (is.na(u)) return("")
    paste0("<a href='", u, "' target='_blank' rel='noopener'>", label, "</a>")
  }

  links <- c(
    make_link(png, "PNG"),
    make_link(pdf, "PDF"),
    make_link(svg, "SVG")
  )
  links <- links[nzchar(links)]

  paste0(
    "<div class='figure-card'>",
    "<img src='", png_url, "' alt='", .echogo_escape_html(alt), "'/>",
    "<div class='figure-caption'><strong>", .echogo_escape_html(title), ".</strong> ", caption, "</div>",
    if (length(links)) paste0(
      "<div class='download-bar'><strong>Full-resolution figure:</strong> ",
      paste(links, collapse = " "), "</div>"
    ) else "",
    "</div>"
  )
}

.echogo_report_table_details <- function(title, data, href = NULL, preview_rows = 12L) {
  if (is.null(data) || !is.data.frame(data) || !nrow(data)) {
    return(.echogo_html_details(title, "<em>No corresponding output was generated.</em>"))
  }

  preview <- utils::head(data, preview_rows)
  tab <- paste(
    capture.output(knitr::kable(preview, format = "html", escape = TRUE, digits = 4)),
    collapse = "\n"
  )

  dl <- if (!is.null(href) && !is.na(href) && nzchar(href)) {
    paste0(
      "<div class='table-download'><a href='", href,
      "' target='_blank'>Open complete table</a> - ",
      format(nrow(data), big.mark = ","), " rows.</div>"
    )
  } else {
    paste0(
      "<div class='table-download'>", format(nrow(data), big.mark = ","),
      " rows; preview shown.</div>"
    )
  }

  .echogo_html_details(title, paste0(dl, tab))
}

.echogo_report_semantic_badge <- function(orgdb, role, method = "Rel") {
  if (is.null(orgdb) || is.na(orgdb) || !nzchar(orgdb)) return("")
  role <- as.character(role %||% "")
  class_name <- if (identical(role, "proxy")) "proxy" else "target-reference"
  role_label <- if (identical(role, "proxy")) "proxy" else "target reference"
  paste0(
    "<span class='semantic-badge ", class_name, "'>Semantic reference: ",
    .echogo_escape_html(orgdb), " | role: ", role_label,
    " | method: ", .echogo_escape_html(method), "</span>"
  )
}

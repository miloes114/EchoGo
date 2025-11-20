#' Internal: HTML <details> wrapper
#' @keywords internal
.echogo_html_details <- function(summary, inner_html) {
  paste0(
    "<details style='margin:8px 0;'>",
    "<summary style='cursor:pointer;font-weight:600;'>", summary, "</summary>",
    "<div style='margin-top:10px;'>", inner_html, "</div>",
    "</details>"
  )
}

#' Internal: turn an absolute path into a browser-friendly URL
#' If target and report_dir share the same base_dir, return a relative URL.
#' Otherwise return a file:/// URL.
#' @keywords internal
.echogo_path_for_html <- function(p, report_dir, base_dir) {
  if (is.null(p) || is.na(p) || !nzchar(p)) return(NA_character_)
  p  <- normalizePath(p,      winslash = "/", mustWork = FALSE)
  rd <- normalizePath(report_dir, winslash = "/", mustWork = FALSE)
  bd <- normalizePath(base_dir,   winslash = "/", mustWork = FALSE)

  if (nzchar(bd) && startsWith(p, paste0(bd, "/")) && startsWith(rd, paste0(bd, "/"))) {
    # both under base_dir -> compute relative path from report_dir to p
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
  if (is.null(path) || is.na(path) || !file.exists(path)) {
    return(.echogo_html_details(title, sprintf("<em>Missing file:</em> %s", path)))
  }
  url <- .echogo_path_for_html(path, report_dir, base_dir)
  inner <- sprintf("<object data='%s' type='application/pdf' width='100%%' height='%s'>", url, height)
  inner <- paste0(inner,
                  "<p>PDF preview not supported. <a href='", url, "' target='_blank'>Open</a>.</p>",
                  "</object>")
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
  if (is.null(path_abs) || is.na(path_abs) || !file.exists(path_abs)) {
    return(.echogo_html_details(title, "<em>Not found in index or not generated for this run.</em>"))
  }
  url <- .echogo_path_for_html(path_abs, report_dir, base_dir)
  iframe <- sprintf(
    "<iframe src='%s' data-external='1' style='width:100%%;height:%s;border:1px solid #ddd;border-radius:6px;'></iframe>",
    url, height
  )
  .echogo_html_details(title, iframe)
}

#' Embed an HTML file (iframe) with a collapsible toggle (deprecated alias)
#'
#' Deprecated. Use [embed_html_toggle_external_plus()] instead.
#' @inheritParams embed_html_toggle_external_plus
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


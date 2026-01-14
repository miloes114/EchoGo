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

  norm <- function(x) {
    x <- normalizePath(x, winslash = "/", mustWork = FALSE)
    x <- gsub("/+$", "", x)
    x
  }

  p  <- norm(p)
  rd <- norm(report_dir)
  bd <- norm(base_dir)

  # Windows: compare case-insensitively
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

  # FIX: absolute existence check first
  if (is.null(path) || is.na(path) || !nzchar(path)) {
    return(.echogo_html_details(title, "<em>Missing file:</em> (no path)"))
  }
  if (!file.exists(path)) {
    return(.echogo_html_details(title, sprintf("<em>Missing file:</em> %s", path)))
  }

  url <- .echogo_stage_asset(path, report_dir = report_dir, base_dir = base_dir)
  if (is.na(url)) {
    return(.echogo_html_details(title, sprintf("<em>Failed to stage file:</em> %s", path)))
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

  # FIX: absolute existence check first
  if (is.null(path_abs) || is.na(path_abs) || !nzchar(path_abs)) {
    return(.echogo_html_details(title, "<em>Missing HTML:</em> (no path)</em>"))
  }
  if (!file.exists(path_abs)) {
    return(.echogo_html_details(title, "<em>Not found in index or not generated for this run.</em>"))
  }

  url <- .echogo_stage_asset(path_abs, report_dir = report_dir, base_dir = base_dir)
  if (is.na(url)) {
    return(.echogo_html_details(title, sprintf("<em>Failed to stage HTML:</em> %s", path_abs)))
  }

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

#' Internal: stage a file into report/assets and return a relative URL
#' - If file is under base_dir, preserve its relative subpath under assets/
#' - Otherwise stage into assets/external/
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

  # Windows: compare case-insensitively
  p_key  <- if (.Platform$OS.type == "windows") tolower(p)  else p
  bd_key <- if (.Platform$OS.type == "windows") tolower(bd) else bd

  prefix <- paste0(bd_key, "/")
  if (nzchar(bd_key) && startsWith(p_key, prefix)) {
    rel <- substr(p, nchar(prefix) + 1, nchar(p))
    dst <- file.path(assets_dir, rel)
    dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    file.copy(path_abs, dst, overwrite = TRUE)
    return(gsub("\\\\", "/", file.path("assets", rel)))
  }

  dst <- file.path(assets_dir, "external", basename(path_abs))
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  file.copy(path_abs, dst, overwrite = TRUE)
  gsub("\\\\", "/", file.path("assets", "external", basename(path_abs)))
}


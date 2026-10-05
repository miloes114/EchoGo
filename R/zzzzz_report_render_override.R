# Late-loaded report-render compatibility layer.
# It repairs RRvGO PNGs, applies report copy/print polish, and validates image
# assets used by the rendered report.

.render_echogo_report <- function(report_title,
                                  template = NULL,
                                  params   = list(),
                                  theme    = "flatly",
                                  outdir,
                                  keep_temp_rmd = FALSE,
                                  verbose_render = FALSE) {

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    warning("rmarkdown not installed; skipping report generation.")
    return(NA_character_)
  }
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  base_dir <- normalizePath(outdir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(base_dir)) { warning("Base dir not found: ", base_dir); return(NA_character_) }
  report_dir <- file.path(base_dir, "report")

  write_index_dir <- function(dir_path) {
    if (!nzchar(dir_path) || !dir.exists(dir_path)) return(invisible(FALSE))
    if (requireNamespace("EchoGO", quietly = TRUE) &&
        "echogo_write_index" %in% getNamespaceExports("EchoGO")) {
      try(EchoGO::echogo_write_index(dir_path), silent = TRUE)
      return(invisible(TRUE))
    }
    rels <- list.files(dir_path, recursive = TRUE, all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
    idx <- tibble::tibble(
      rel_path = gsub("\\\\", "/", rels),
      full_path = gsub("\\\\", "/", normalizePath(file.path(dir_path, rels), winslash = "/", mustWork = FALSE))
    )
    out_csv <- file.path(dir_path, "__file_index.csv")
    if (requireNamespace("readr", quietly = TRUE)) readr::write_csv(idx, out_csv) else utils::write.csv(idx, out_csv, row.names = FALSE)
    invisible(TRUE)
  }
  write_index_dir(base_dir)

  tpl_src <- template %||% system.file("reports", "echogo_report.Rmd", package = "EchoGO")
  if (!nzchar(tpl_src) || !file.exists(tpl_src)) {
    tpl_src <- file.path(base_dir, "echogo_report_minimal.Rmd")
    cat(
      '---\ntitle: "', report_title, '"\noutput:\n  html_document:\n    theme: ', theme,
      '\nparams:\n  dirs: NULL\n  top_n: 25\n  sections: ["overview","goseq","gprofiler","consensus","rrvgo","networks","evaluation"]\n---',
      file = tpl_src, sep = ""
    )
  }

  stable_rmd <- file.path(base_dir, "echogo_report.Rmd")
  ok_copy <- file.copy(tpl_src, stable_rmd, overwrite = TRUE, copy.date = TRUE)
  if (!ok_copy || !file.exists(stable_rmd)) {
    warning("Failed to stage Rmd in base: ", stable_rmd)
    return(NA_character_)
  }

  css_src <- file.path(dirname(tpl_src), "echogo_report.css")
  css_staged <- file.path(base_dir, "echogo_report.css")
  if (file.exists(css_src)) file.copy(css_src, css_staged, overwrite = TRUE, copy.date = TRUE)

  logo_source <- tryCatch(.echogo_logo_source_path(), error = function(e) NA_character_)
  logo_staged <- tryCatch(.echogo_stage_report_logo(base_dir, logo_source), error = function(e) NA_character_)
  if (is.na(logo_staged) || !file.exists(logo_staged)) {
    warning("Official EchoGO logo was not found; rendering without branding asset.")
  }
  on.exit({
    if (file.exists(css_staged)) unlink(css_staged, force = TRUE)
    if (file.exists(file.path(base_dir, "echogo-logo.png"))) unlink(file.path(base_dir, "echogo-logo.png"), force = TRUE)
  }, add = TRUE)

  # Repair or create RRvGO treemap PNG sidecars before the Rmd discovers its
  # primary semantic figures.  Invalid PNGs are removed if they cannot be
  # repaired, allowing the report to fall back to the valid PDF rather than
  # embedding a broken image.
  rrvgo_png_audit <- tryCatch(
    .echogo_ensure_rrvgo_report_pngs(file.path(base_dir, "rrvgo")),
    error = function(e) tibble::tibble()
  )
  if (is.data.frame(rrvgo_png_audit) && nrow(rrvgo_png_audit)) {
    bad <- which(!rrvgo_png_audit$valid_after)
    if (length(bad)) {
      for (i in bad) {
        p <- rrvgo_png_audit$png_file[[i]]
        if (!is.na(p) && file.exists(p)) unlink(p, force = TRUE)
      }
    }
  }

  render_params <- params %||% list()
  render_params$dirs <- render_params$dirs %||% list()
  render_params$dirs$base <- base_dir
  render_params$dirs$report <- report_dir

  yaml <- tryCatch(rmarkdown::yaml_front_matter(stable_rmd), error = function(e) NULL)
  if (!is.null(yaml) && is.list(yaml$params)) {
    allowed <- names(yaml$params)
    if (length(allowed)) render_params <- render_params[intersect(names(render_params), allowed)]
  }
  if ("sections" %in% names(render_params)) {
    s <- tolower(unlist(render_params$sections, use.names = FALSE))
    s <- unique(c(s, "evaluation"))
    render_params$sections <- s
  }

  old_opt <- options(echogo.outdir = base_dir); on.exit(options(old_opt), add = TRUE)
  old_env <- Sys.getenv("ECHOGO_OUTDIR", unset = NA_character_)
  Sys.setenv(ECHOGO_OUTDIR = base_dir)
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("ECHOGO_OUTDIR") else Sys.setenv(ECHOGO_OUTDIR = old_env)
  }, add = TRUE)

  ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
  safe_title <- gsub("[^A-Za-z0-9_-]+", "-", report_title)
  html_name <- paste0(safe_title, "_", ts, ".html")
  html_in_root <- file.path(base_dir, html_name)
  log_file <- file.path(base_dir, "__report_render.log")
  if (file.exists(log_file)) unlink(log_file, force = TRUE)

  old_wd <- setwd(base_dir); on.exit({ setwd(old_wd) }, add = TRUE)
  render_ok <- TRUE
  res_path <- tryCatch({
    log_txt <- capture.output(
      rmarkdown::render(
        input = basename(stable_rmd),
        output_file = basename(html_in_root),
        output_dir = base_dir,
        intermediates_dir = base_dir,
        knit_root_dir = base_dir,
        params = render_params,
        envir = new.env(parent = globalenv()),
        quiet = !isTRUE(verbose_render),
        encoding = "UTF-8"
      ),
      type = "output"
    )
    if (length(log_txt)) writeLines(log_txt, log_file)
    html_in_root
  }, error = function(e) {
    render_ok <<- FALSE
    writeLines(paste0("ERROR: ", conditionMessage(e)), log_file)
    warning("Report render failed: ", conditionMessage(e))
    NA_character_
  })

  if (!render_ok || is.na(res_path) || !file.exists(res_path)) {
    warning("Report HTML not produced. See log: ", log_file)
    return(NA_character_)
  }

  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(css_src)) file.copy(css_src, file.path(report_dir, "echogo_report.css"), overwrite = TRUE, copy.date = TRUE)
  if (!is.na(logo_source) && file.exists(logo_source)) file.copy(logo_source, file.path(report_dir, "echogo-logo.png"), overwrite = TRUE, copy.date = TRUE)

  final_html <- file.path(report_dir, basename(res_path))
  if (!file.rename(res_path, final_html)) {
    file.copy(res_path, final_html, overwrite = TRUE)
    unlink(res_path, force = TRUE)
  }

  final_rmd <- file.path(report_dir, "echogo_report.Rmd")
  if (!file.rename(stable_rmd, final_rmd)) {
    file.copy(stable_rmd, final_rmd, overwrite = TRUE)
    unlink(stable_rmd, force = TRUE)
  }

  if (file.exists(css_staged)) file.copy(css_staged, file.path(report_dir, "echogo_report.css"), overwrite = TRUE, copy.date = TRUE)
  if (!is.na(logo_staged) && file.exists(logo_staged)) file.copy(logo_staged, file.path(report_dir, "echogo-logo.png"), overwrite = TRUE, copy.date = TRUE)

  root_index <- file.path(base_dir, "__file_index.csv")
  if (file.exists(root_index)) file.copy(root_index, file.path(report_dir, "__file_index.csv"), overwrite = TRUE)

  if (file.exists(log_file)) {
    if (!file.rename(log_file, file.path(report_dir, "__report_render.log"))) {
      file.copy(log_file, file.path(report_dir, "__report_render.log"), overwrite = TRUE)
      unlink(log_file, force = TRUE)
    }
  }

  # Final visible-copy cleanup + print behavior, followed by a hard report
  # image-integrity gate.  A product marked/generated but referenced through a
  # missing or corrupt image is now a validation error rather than a broken
  # placeholder in the final report.
  .echogo_final_report_html_polish(final_html)
  .echogo_validate_report_image_assets(final_html, report_dir)

  normalizePath(final_html, winslash = "/", mustWork = FALSE)
}

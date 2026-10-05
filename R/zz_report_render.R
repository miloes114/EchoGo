# Final report renderer override ------------------------------------------------
# Loaded after report_render.R so the universal v0.1.4 biological report owns
# the current rendering contract without changing the exported API.

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
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(base_dir)) {
    warning("Base dir not found: ", base_dir)
    return(NA_character_)
  }

  report_dir <- file.path(base_dir, "report")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

  write_index_dir <- function(dir_path) {
    if (!nzchar(dir_path) || !dir.exists(dir_path)) return(invisible(FALSE))
    if (requireNamespace("EchoGO", quietly = TRUE) &&
        "echogo_write_index" %in% getNamespaceExports("EchoGO")) {
      try(EchoGO::echogo_write_index(dir_path), silent = TRUE)
      return(invisible(TRUE))
    }
    rels <- list.files(
      dir_path, recursive = TRUE, all.files = TRUE,
      include.dirs = FALSE, no.. = TRUE
    )
    idx <- tibble::tibble(
      rel_path = gsub("\\\\", "/", rels),
      full_path = gsub(
        "\\\\", "/",
        normalizePath(file.path(dir_path, rels), winslash = "/", mustWork = FALSE)
      )
    )
    out_csv <- file.path(dir_path, "__file_index.csv")
    if (requireNamespace("readr", quietly = TRUE)) {
      readr::write_csv(idx, out_csv)
    } else {
      utils::write.csv(idx, out_csv, row.names = FALSE)
    }
    invisible(TRUE)
  }
  write_index_dir(base_dir)

  tpl_src <- template %||% system.file("reports", "echogo_report.Rmd", package = "EchoGO")
  if (!nzchar(tpl_src) || !file.exists(tpl_src)) {
    warning("EchoGO report template was not found; skipping report generation.")
    return(NA_character_)
  }

  stable_rmd <- file.path(base_dir, "echogo_report.Rmd")
  if (!file.copy(tpl_src, stable_rmd, overwrite = TRUE, copy.date = TRUE)) {
    warning("Failed to stage Rmd in base: ", stable_rmd)
    return(NA_character_)
  }

  css_src <- system.file("reports", "echogo_report.css", package = "EchoGO")
  css_staged <- file.path(base_dir, "echogo_report.css")
  css_available <- nzchar(css_src) && file.exists(css_src) &&
    isTRUE(file.copy(css_src, css_staged, overwrite = TRUE, copy.date = TRUE))
  if (!css_available) {
    warning("EchoGO report CSS was not found; rendering without custom report styling.")
  }

  # Stage the exact official logo beside the temporary Rmd so self-contained
  # rendering embeds it without an absolute local path.
  logo_source <- tryCatch(.echogo_logo_source_path(), error = function(e) NA_character_)
  logo_staged <- tryCatch(
    .echogo_stage_report_logo(base_dir, logo_source),
    error = function(e) NA_character_
  )
  if (is.na(logo_staged) || !file.exists(logo_staged)) {
    warning("Official EchoGO logo was not found; rendering without branding asset.")
  }

  render_params <- params %||% list()
  render_params$dirs <- render_params$dirs %||% list()
  render_params$dirs$base <- base_dir
  render_params$dirs$report <- report_dir

  yaml <- tryCatch(rmarkdown::yaml_front_matter(stable_rmd), error = function(e) NULL)
  if (!is.null(yaml) && is.list(yaml$params)) {
    allowed <- names(yaml$params)
    if (length(allowed)) {
      render_params <- render_params[intersect(names(render_params), allowed)]
    }
  }

  old_opt <- options(echogo.outdir = base_dir)
  on.exit(options(old_opt), add = TRUE)
  old_env <- Sys.getenv("ECHOGO_OUTDIR", unset = NA_character_)
  Sys.setenv(ECHOGO_OUTDIR = base_dir)
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("ECHOGO_OUTDIR") else Sys.setenv(ECHOGO_OUTDIR = old_env)
  }, add = TRUE)

  ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
  safe_title <- gsub("[^A-Za-z0-9_-]+", "-", report_title)
  html_name <- paste0(safe_title, "_", ts, ".html")
  final_html <- file.path(report_dir, html_name)
  log_file <- file.path(base_dir, "__report_render.log")
  if (file.exists(log_file)) unlink(log_file, force = TRUE)

  old_wd <- setwd(base_dir)
  on.exit(setwd(old_wd), add = TRUE)

  render_ok <- TRUE
  res_path <- tryCatch({
    log_txt <- capture.output(
      rmarkdown::render(
        input = basename(stable_rmd),
        output_file = html_name,
        output_dir = report_dir,
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
    final_html
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

  final_rmd <- file.path(report_dir, "echogo_report.Rmd")
  if (!file.rename(stable_rmd, final_rmd)) {
    file.copy(stable_rmd, final_rmd, overwrite = TRUE)
    unlink(stable_rmd, force = TRUE)
  }

  if (file.exists(css_src)) {
    file.copy(css_src, file.path(report_dir, "echogo_report.css"), overwrite = TRUE, copy.date = TRUE)
  }
  if (!is.na(logo_source) && file.exists(logo_source)) {
    file.copy(logo_source, file.path(report_dir, "echogo-logo.png"), overwrite = TRUE, copy.date = TRUE)
  }

  root_index <- file.path(base_dir, "__file_index.csv")
  if (file.exists(root_index)) {
    file.copy(root_index, file.path(report_dir, "__file_index.csv"), overwrite = TRUE)
  }
  if (file.exists(log_file)) {
    target_log <- file.path(report_dir, "__report_render.log")
    if (!file.rename(log_file, target_log)) {
      file.copy(log_file, target_log, overwrite = TRUE)
      unlink(log_file, force = TRUE)
    }
  }

  if (file.exists(css_staged)) unlink(css_staged, force = TRUE)
  if (file.exists(file.path(base_dir, "echogo-logo.png"))) unlink(file.path(base_dir, "echogo-logo.png"), force = TRUE)
  write_index_dir(base_dir)

  normalizePath(final_html, winslash = "/", mustWork = FALSE)
}

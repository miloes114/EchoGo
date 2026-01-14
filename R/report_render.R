#' Internal: render the EchoGO HTML report (robust: knit in root, then move; clean root on success)
#' @keywords internal
.render_echogo_report <- function(report_title,
                                  template = NULL,
                                  params   = list(),
                                  theme    = "flatly",
                                  outdir,
                                  keep_temp_rmd = FALSE,   # ignored; we keep only stable Rmd in report/
                                  verbose_render = FALSE) {

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    warning("rmarkdown not installed; skipping report generation.")
    return(NA_character_)
  }
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # --- Base (root) & Report dir
  base_dir   <- normalizePath(outdir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(base_dir)) { warning("Base dir not found: ", base_dir); return(NA_character_) }
  report_dir <- file.path(base_dir, "report")  # will be created only after successful knit

  # --- Write root index (needed while knitting)
  write_index_dir <- function(dir_path) {
    if (!nzchar(dir_path) || !dir.exists(dir_path)) return(invisible(FALSE))
    if (requireNamespace("EchoGO", quietly = TRUE) &&
        "echogo_write_index" %in% getNamespaceExports("EchoGO")) {
      try(EchoGO::echogo_write_index(dir_path), silent = TRUE)
      return(invisible(TRUE))
    }
    rels <- list.files(dir_path, recursive = TRUE, all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
    idx  <- tibble::tibble(
      rel_path  = gsub("\\\\","/", rels),
      full_path = gsub("\\\\","/", normalizePath(file.path(dir_path, rels),
                                                 winslash="/", mustWork=FALSE))
    )
    out_csv <- file.path(dir_path, "__file_index.csv")
    if (requireNamespace("readr", quietly = TRUE)) readr::write_csv(idx, out_csv) else utils::write.csv(idx, out_csv, row.names = FALSE)
    invisible(TRUE)
  }
  write_index_dir(base_dir)

  # --- Resolve template
  tpl_src <- template %||% system.file("reports", "echogo_report.Rmd", package = "EchoGO")
  if (!nzchar(tpl_src) || !file.exists(tpl_src)) {
    tpl_src <- file.path(base_dir, "echogo_report_minimal.Rmd")
    cat(
      '---
title: "', report_title, '"
output:
  html_document:
    theme: ', theme, '
params:
  dirs: NULL
  top_n: 25
  sections: ["overview","goseq","gprofiler","consensus","rrvgo","networks","evaluation"]
---', file = tpl_src, sep = ""
    )
  }

  # --- Use a STABLE Rmd name in root (what the user sees & can re-knit)
  stable_rmd <- file.path(base_dir, "echogo_report.Rmd")
  ok_copy    <- file.copy(tpl_src, stable_rmd, overwrite = TRUE, copy.date = TRUE)
  if (!ok_copy || !file.exists(stable_rmd)) {
    warning("Failed to stage Rmd in base: ", stable_rmd)
    return(NA_character_)
  }

  # --- Params for the Rmd; tell it where the FINAL report dir will be
  render_params <- params %||% list()
  render_params$dirs <- render_params$dirs %||% list()
  render_params$dirs$base   <- base_dir
  render_params$dirs$report <- report_dir

  # NEW: keep only params that are declared in the Rmd's YAML
  yaml <- tryCatch(rmarkdown::yaml_front_matter(stable_rmd), error = function(e) NULL)
  if (!is.null(yaml) && is.list(yaml$params)) {
    allowed <- names(yaml$params)
    if (length(allowed)) {
      render_params <- render_params[intersect(names(render_params), allowed)]
    }
  }

  # MINIMAL FIX: if 'sections' is provided, normalize + force-include 'evaluation'
  if ("sections" %in% names(render_params)) {
    s <- tolower(unlist(render_params$sections, use.names = FALSE))
    s <- unique(c(s, "evaluation"))
    render_params$sections <- s
  }
  # (If 'sections' is not provided, the Rmd's default includes 'evaluation' already.)

  # --- Env exports used by helpers
  old_opt <- options(echogo.outdir = base_dir); on.exit(options(old_opt), add = TRUE)
  old_env <- Sys.getenv("ECHOGO_OUTDIR", unset = NA_character_)
  Sys.setenv(ECHOGO_OUTDIR = base_dir)
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("ECHOGO_OUTDIR") else Sys.setenv(ECHOGO_OUTDIR = old_env)
  }, add = TRUE)

  # --- Knit in root; HTML name is timestamped
  ts          <- format(Sys.time(), "%Y%m%d-%H%M%S")
  safe_title  <- gsub("[^A-Za-z0-9_-]+", "-", report_title)
  html_name   <- paste0(safe_title, "_", ts, ".html")
  html_in_root<- file.path(base_dir, html_name)
  log_file    <- file.path(base_dir, "__report_render.log")
  if (file.exists(log_file)) unlink(log_file, force = TRUE)

  old_wd <- setwd(base_dir); on.exit({ setwd(old_wd) }, add = TRUE)

  # Render (no fragile sinks; capture output to log)
  render_ok <- TRUE
  res_path <- tryCatch({
    log_txt <- capture.output(
      rmarkdown::render(
        input             = basename(stable_rmd),
        output_file       = basename(html_in_root),
        output_dir        = base_dir,
        intermediates_dir = base_dir,
        knit_root_dir     = base_dir,
        params            = render_params,
        envir             = new.env(parent = globalenv()),
        quiet             = !isTRUE(verbose_render),
        encoding          = "UTF-8"
      ),
      type = "output"
    )
    if (length(log_txt)) writeLines(log_txt, log_file)
    html_in_root
  }, error = function(e) {
    render_ok <<- FALSE
    writeLines(paste0("ERROR: ", conditionMessage(e)), log_file)
    warning("Report render failed: ", conditionMessage(e))
    return(NA_character_)
  })

  if (!render_ok || is.na(res_path) || !file.exists(res_path)) {
    warning("Report HTML not produced. See log: ", log_file)
    return(NA_character_)
  }

  # --- Move artifacts into report/, then clean root
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

  # Move HTML
  final_html <- file.path(report_dir, basename(res_path))
  if (!file.rename(res_path, final_html)) {
    file.copy(res_path, final_html, overwrite = TRUE)
    unlink(res_path, force = TRUE)
  }

  # Move Rmd (stable)
  final_rmd <- file.path(report_dir, "echogo_report.Rmd")
  if (!file.rename(stable_rmd, final_rmd)) {
    file.copy(stable_rmd, final_rmd, overwrite = TRUE)
    unlink(stable_rmd, force = TRUE)
  }

  # Copy index into report/ but KEEP the authoritative index in base_dir
  root_index <- file.path(base_dir, "__file_index.csv")
  if (file.exists(root_index)) {
    file.copy(root_index, file.path(report_dir, "__file_index.csv"), overwrite = TRUE)
  }

  if (file.exists(log_file)) {
    if (!file.rename(log_file, file.path(report_dir, "__report_render.log"))) {
      file.copy(log_file, file.path(report_dir, "__report_render.log"), overwrite = TRUE)
      unlink(log_file, force = TRUE)
    }
  }

  # Fresh indices (optional)
  #try(write_index_dir(report_dir), silent = TRUE)
  # Root is now clean of report artefacts

  normalizePath(final_html, winslash = "/", mustWork = FALSE)
}


#' EchoGO — Quickstart help
#'
#' Print a short overview of the main EchoGO entry points, vignettes,
#' and helper functions.
#'
#' @export
echogo_help <- function() {
  cat("
EchoGO — Quickstart

Documentation & vignettes:
  browseVignettes('EchoGO')                   # list all EchoGO vignettes
  vignette('EchoGO_workflow')                 # workflow & running the pipeline
  vignette('EchoGO_interpretation')           # interpreting consensus, RRvGO & networks

Run the built-in demo (fast: 3 species + mouse OrgDb):
  echogo_quickstart(run_demo = TRUE)

Use your own data:
  echogo_scaffold('my_project')
  # Put files into my_project/input (see _README.txt there)
  echogo_run(input_dir = 'my_project/input', outdir = 'my_project/results')

Choosing species codes (for g:Profiler):
  EchoGO uses g:Profiler organism IDs such as:
    hsapiens, mmusculus, drerio, dmelanogaster, celegans, rnorvegicus...
  Full organism list (external):
    https://biit.cs.ut.ee/gprofiler/page/organism-list

  EchoGO helpers:
    echogo_pick_species()                 # interactive table to browse/filter and copy IDs
    echogo_preflight_species(c('hsapiens','mmusculus'))
                                          # validate/suggest IDs before run_full_echogo()

Checking / installing OrgDb & GO.db:
  echogo_list_orgdb()                     # see which OrgDb packages are installed
  echogo_install_orgdb_instructions()     # print BiocManager::install() commands
  echogo_install_orgdb()                  # automatically install into user library
  echogo_require_orgdb()                  # ensure OrgDb is available before running

Demo data:
  EchoGO ships with example inputs and frozen results.
  Use:
    echogo_open_demo()                    # open demo input/results folders and list files

Run the full workflow:
  run_full_echogo(
    input_dir,
    species = getOption('EchoGO.default_species', c('hsapiens','mmusculus','drerio')),
    orgdb   = getOption('EchoGO.default_orgdb', 'org.Mm.eg.db'),
    outdir  = 'echogo_out',
    make_report = TRUE
  )

See also:
  ?echogo_pick_species
  ?echogo_preflight_species
  ?echogo_list_orgdb
  ?echogo_install_orgdb
  ?echogo_require_orgdb
  ?echogo_open_demo
")
}


#' Copy the built-in demo and (optionally) run it end-to-end
#'
#' Copies files from the package's `inst/extdata/echogo_demo` into `outdir/echogo_demo`.
#' If `run_demo = TRUE`, it runs `run_full_echogo(input_dir = <that folder>)`
#' using a small 3-species set (hsapiens, mmusculus, drerio) and OrgDb org.Mm.eg.db, with report ON.
#'
#' @param run_demo logical; run the pipeline after copying (default: FALSE)
#' @param outdir where to copy the demo (default: tempdir())
#' @param make_report logical; build the HTML report (default: TRUE)
#' @param ... passed through to `run_full_echogo()`
#' @export
echogo_quickstart <- function(run_demo = FALSE,
                              outdir = getOption("EchoGO.demo_root", default = tempdir()),
                              make_report = TRUE,
                              ...) {
  src <- system.file("extdata", "echogo_demo", package = "EchoGO")
  if (src == "") stop("Demo data not found in inst/extdata/echogo_demo of the package.")

  demo_dir <- file.path(outdir, "echogo_demo")
  dir.create(demo_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(src, full.names = TRUE)
  ok <- file.copy(files, demo_dir, recursive = TRUE, overwrite = TRUE)
  if (!all(ok)) warning("Some demo files could not be copied.")

  message("Demo data copied to: ", normalizePath(demo_dir, winslash = "/"))

  if (!isTRUE(run_demo)) return(invisible(demo_dir))

  # Fixed small demo set for speed
  species_demo <- c("hsapiens","mmusculus","drerio")
  orgdb_demo   <- "org.Mm.eg.db"

  # 🔹 Ensure required OrgDb is present (auto-install if missing)
  echogo_require_orgdb(c("GO.db", orgdb_demo), auto_install = TRUE)

  res <- run_full_echogo(
    input_dir   = demo_dir,
    species     = species_demo,
    orgdb       = orgdb_demo,
    outdir      = file.path(demo_dir, "results"),
    make_report = isTRUE(make_report),
    ...
  )

  if (interactive() && !is.null(res$files$report_html) && !is.na(res$files$report_html)) {
    utils::browseURL(res$files$report_html)
  }
  invisible(res)
}


#' Scaffold a minimal EchoGO project directory
#' @param path target directory
#' @export
echogo_scaffold <- function(path = "echogo_project") {
  input  <- file.path(path, "input")
  output <- file.path(path, "results")
  dir.create(input,  recursive = TRUE, showWarnings = FALSE)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)

  readme <- "
Place your files in input/ with these names (or set them in input/config.yml):
  - gene.counts.matrix.tsv
  - DE_*.tsv (one or more differential expression tables)
  - Trinotate.xls

Then run:
  echogo_run(input_dir = 'PATH/input', outdir = 'PATH/results')

File format notes:
  * counts: TSV with ID in first column and counts in remaining columns
  * DE_*: TSV with at least columns: id, log2FC, pvalue, padj
  * Trinotate.xls: standard Trinotate report

Species & annotation:
  EchoGO uses g:Profiler organism IDs (e.g., hsapiens, mmusculus, drerio).
  Full organism list:
    https://biit.cs.ut.ee/gprofiler/page/organism-list

  Helpful R helpers (run these in R after loading EchoGO):
    echogo_help()                        # overall quickstart
    echogo_list_orgdb()                  # check installed OrgDb packages
    echogo_install_orgdb_instructions()  # show BiocManager::install() commands
    echogo_install_orgdb()               # auto-install GO.db / OrgDb into user library
    echogo_require_orgdb()               # ensure OrgDb is available before running EchoGO
"

  writeLines(readme, file.path(input, "_README.txt"))

  yaml <- "species: [hsapiens, mmusculus, drerio]
orgdb: org.Mm.eg.db
report_title: EchoGO Report
"
  writeLines(yaml, file.path(input, "config.yml"))

  message("Scaffold created at: ", normalizePath(path, winslash = "/"))
  invisible(path)
}


#' Run EchoGO from a folder of inputs (and optional YAML)
#'
#' This wrapper prefers the modern `run_full_echogo(input_dir = ...)` interface.
#' It reads `input/config.yml` when present to pick species/orgdb.
#'
#' @param input_dir folder containing your inputs (see scaffold README)
#' @param outdir results directory (default sibling 'results')
#' @param config optional YAML (defaults to input/config.yml)
#' @param ... passed to `run_full_echogo()`
#' @export
echogo_run <- function(input_dir,
                       outdir  = file.path(dirname(normalizePath(input_dir)), "results"),
                       config  = file.path(input_dir, "config.yml"),
                       ...) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  cfg <- list()
  if (file.exists(config)) {
    if (!requireNamespace("yaml", quietly = TRUE)) stop("Install the 'yaml' package to read ", config)
    cfg <- yaml::read_yaml(config)
  }

  species <- cfg$species %||% getOption("EchoGO.default_species",
                                        c("hsapiens","mmusculus","drerio"))
  orgdb   <- cfg$orgdb   %||% getOption("EchoGO.default_orgdb", "org.Mm.eg.db")

  # 🔹 Ensure OrgDb + GO.db are available (auto-install to user library if needed)
  echogo_require_orgdb(c("GO.db", orgdb), auto_install = TRUE)

  run_full_echogo(input_dir   = input_dir,
                  species     = species,
                  orgdb       = orgdb,
                  outdir      = outdir,
                  make_report = TRUE,
                  ...)
}

#' Open or inspect the shipped EchoGO demo data
#'
#' Prints the paths to the packaged demo inputs and frozen results,
#' shows their contents, and attempts to open the folders in the
#' system file browser (Windows Explorer, macOS Finder, or Linux xdg-open).
#'
#' @export
echogo_open_demo <- function() {
  demo_in  <- echogo_demo_path()
  demo_out <- echogo_demo_results_path()

  cat("\nEchoGO demo inputs:\n  ", demo_in, "\n\n")
  print(list.files(demo_in, full.names = TRUE))

  cat("\nEchoGO frozen demo results:\n  ", demo_out, "\n\n")
  print(list.files(demo_out, full.names = TRUE))

  open_dir <- function(path) {
    if (.Platform$OS.type == "windows") {
      shell.exec(path)
    } else if (Sys.info()[["sysname"]] == "Darwin") {
      system2("open", path)
    } else {
      system2("xdg-open", path)
    }
  }

  cat("\nOpening demo directories...\n")
  try(open_dir(demo_in),  silent = TRUE)
  try(open_dir(demo_out), silent = TRUE)

  invisible(list(inputs = demo_in, results = demo_out))
}

#' Path to shipped demo inputs
#' @return A filesystem path to inst/extdata/echogo_demo
#' @export
echogo_demo_path <- function() {
  system.file("extdata", "echogo_demo", package = "EchoGO")
}

#' Path to shipped demo results (frozen snapshot)
#' @return A filesystem path to inst/extdata/echogo_demo_results
#' @export
echogo_demo_results_path <- function() {
  system.file("extdata", "echogo_demo_results", package = "EchoGO")
}

#' EchoGO - Quickstart help
#'
#' Prints offline-demo and own-data guidance. The standard workflow requires
#' precomputed target GOseq as its primary target evidence, DE results,
#' a complete tested-feature background and compatible annotation. Own-data
#' runs require researcher-selected annotation contexts and an explicit target
#' or no-target g:Profiler declaration; that declaration does not remove GOseq.
#' When RRvGO is enabled, declare its semantic reference and reference role.
#' @export
echogo_help <- function() {
  cat("
EchoGO - Quickstart

Run the built-in demo (light first-run validation; offline by default):
  echogo_quickstart(run_demo = TRUE)

Run the extended demo (default-domain exploration, RRvGO and diagnostics):
  echogo_quickstart(run_demo = TRUE, full = TRUE)

Use your own data:
  echogo_scaffold('my_project')
  # Populate inputs: DE results, complete tested-feature background,
  # annotation and REQUIRED precomputed target GOseq (see input/_README.txt).
  # Replace this illustrative zebrafish panel for your biological question.
  echogo_run(input_dir = 'my_project/input', outdir = 'my_project/results',
             species = c('drerio', 'omykiss'), target_context = 'drerio',
             run_exploratory_default_domain = FALSE, run_rrvgo = TRUE,
             semantic_reference_orgdb = 'org.Dr.eg.db',
             semantic_reference_role = 'target_reference')
  # Record optional context rationale with context_metadata or config.yml.
  # If no queried context is the target: target_context = NA_character_.
  # Target GOseq stays required; this declares only the g:Profiler context.
  # To skip semantic reduction, set run_rrvgo = FALSE.

Reference-based RNA-seq (DESeq2 + GOseq precomputed):
  If your contrast input folder already contains:
    - allcounts_table.txt
    - dge_<CONTRAST>.csv
    - dge_<CONTRAST>.GOseq.enriched.tsv
    - Trinotate_for_EchoGO.tsv (or <reference_label>_eggNOG_for_EchoGO.tsv)
  you can run:

    input_dir <- 'path/to/reference_project/echogo_input'
    outdir    <- 'path/to/reference_project/echogo_results'
    species   <- c('drerio','strutta','gaculeatus','olatipes','trubripes','amexicanus',
                   'oniloticus','ssalar','omykiss','okisutch','otshawytscha')

    run_full_echogo(input_dir = input_dir,
                    outdir = outdir,
                    species = species,
                    target_context = 'drerio',
                    semantic_reference_orgdb = 'org.Dr.eg.db',
                    semantic_reference_role = 'target_reference',
                    run_exploratory_default_domain = FALSE,
                    run_evaluation = TRUE,
                    make_report = TRUE,
                    verbose = TRUE)

  To build that folder from a GFF3/GTF and protein FASTA, start with:
    echogo_prepare_reference_inputs(root, gff_file, protein_fasta)
  Run the eggNOG-mapper command it writes, then repeat the same call.

  EchoGO resolves significant and tested genes with the same portable canonical-name contract.

Choosing species codes (for g:Profiler):
  EchoGO uses g:Profiler organism IDs such as:
    hsapiens, mmusculus, drerio, dmelanogaster, celegans, rnorvegicus...
  Full organism list (external):
    https://biit.cs.ut.ee/gprofiler/page/organism-list

  EchoGO helpers:
    echogo_pick_species()                 # interactive table to browse/filter and copy IDs
    echogo_preflight_species(c('hsapiens','mmusculus'))
                                          # validate/suggest species IDs before run_full_echogo()

Checking / installing OrgDb & GO.db:
  - echogo_list_orgdb()                 # see which OrgDb packages are installed
  - echogo_install_orgdb_instructions() # print BiocManager::install() commands
  - echogo_install_orgdb('org.Mm.eg.db')# install one OrgDb into the active library
  - echogo_require_orgdb(c('GO.db','org.Dr.eg.db')) # preflight a declared reference

Demo data:
  EchoGO ships with example inputs and frozen results.
  Use:
    echogo_open_demo()  # opens demo input/results folders and lists files

Run the full workflow:
  run_full_echogo(input_dir,
                  species = c('drerio','strutta'),
                  target_context = 'drerio',
                  semantic_reference_orgdb = 'org.Dr.eg.db',
                  semantic_reference_role = 'target_reference',
                  outdir  = 'echogo_out',
                  make_report = TRUE)

See also:
  ?echogo_pick_species
  ?echogo_preflight_species
  ?echogo_list_orgdb
  ?echogo_install_orgdb
  ?echogo_require_orgdb
  ?echogo_open_demo
  ?echogo_resolve_reference_inputs
  ?echogo_prepare_reference_inputs
")
}


.validate_quickstart_result <- function(res, require_report = FALSE) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  evidence_file <- res$files$exact_term_evidence %||% res$files$consensus_xlsx

  if (is.null(evidence_file) || !file.exists(evidence_file)) {
    stop("EchoGO demo did not produce its exact-term evidence table.", call. = FALSE)
  }

  evidence <- if (grepl("\\.xlsx$", evidence_file, ignore.case = TRUE)) openxlsx::read.xlsx(evidence_file) else readr::read_csv(evidence_file, show_col_types = FALSE)
  if (!nrow(evidence)) {
    stop(
      "EchoGO demo completed but returned an empty exact-term evidence table. ",
      "Inspect the GOseq parsing and g:Profiler input vectors.",
      call. = FALSE
    )
  }

  if (isTRUE(require_report)) {
    report_file <- res$files$report_html
    if (
      is.null(report_file) ||
      length(report_file) != 1L ||
      is.na(report_file) ||
      !file.exists(report_file)
    ) {
      stop("EchoGO demo did not produce its HTML report.", call. = FALSE)
    }
  }

  invisible(TRUE)
}

#' Copy the built-in demo and (optionally) run it end-to-end
#'
#' Copies files from the package's `inst/extdata/echogo_demo` into `outdir/echogo_demo`.
#' If `run_demo = TRUE`, it runs `run_full_echogo(input_dir = <that folder>)`
#' using a small fixed set of researcher-declared annotation contexts. The target
#' context is zebrafish (`drerio`); the full demonstration uses `org.Dr.eg.db`
#' as an explicitly declared target semantic reference.
#' The default is a light installation check. Set `full = TRUE` to include the
#' exploratory stream, RRvGO, and evaluation outputs.
#'
#' @param run_demo logical; run the pipeline after copying (default: FALSE)
#' @param outdir where to copy the demo (default: tempdir())
#' @param make_report logical; build the HTML report (default: TRUE)
#' @param clean logical; remove previous demo results before running (default: TRUE)
#' @param full Logical; run all exploratory, RRvGO, and evaluation stages. The
#'   default `FALSE` is intended as a quicker first-run validation.
#' @param live_gprofiler Logical; if `FALSE` (default), use deterministic cached
#'   demonstration responses. Set `TRUE` for an optional live integration test.
#' @param ... passed through to `run_full_echogo()`
#' @export
echogo_quickstart <- function(run_demo = FALSE,
                              outdir = getOption("EchoGO.demo_root", default = tempdir()),
                              make_report = TRUE,
                              clean = TRUE,
                              full = FALSE,
                              live_gprofiler = FALSE,
                              ...) {
  src <- system.file("extdata", "echogo_demo", package = "EchoGO")
  if (src == "") stop("Demo data not found in inst/extdata/echogo_demo of the package.")

  demo_dir <- file.path(outdir, "echogo_demo")
  .echogo_copy_tree(src, demo_dir, overwrite = TRUE)

  message("Demo data copied to: ", normalizePath(demo_dir, winslash = "/"))

  if (!isTRUE(run_demo)) return(invisible(demo_dir))

  # The built-in demo is explicitly offline unless live_gprofiler = TRUE.
  # Scope both species and taxonomy refresh controls to this invocation so a
  # missing or stale user cache cannot trigger an external request, while
  # preserving the caller's options on success or error.
  offline_options <- NULL
  if (!isTRUE(live_gprofiler)) {
    offline_options <- options(
      EchoGO.species_autoupdate = FALSE,
      EchoGO.taxonomy_online = FALSE
    )
    on.exit(options(offline_options), add = TRUE)
  }

  # This fixed panel is part of the documented demonstration, not an own-data
  # default. Names become context labels in all user-facing output.
  species_demo <- c(zebrafish = "drerio", mouse = "mmusculus", human = "hsapiens")
  context_metadata_demo <- data.frame(
    context_code = c("drerio", "mmusculus", "hsapiens"),
    context_rationale = c(
      "Experimental target annotation context.",
      "Selected mammalian annotation context for comparative follow-up.",
      "Selected human annotation context for comparative follow-up."
    ),
    stringsAsFactors = FALSE
  )
  semantic_reference_demo <- "org.Dr.eg.db"

  demo_files <- list(
    goseq = file.path(demo_dir, "GOseq_enrichment_demo.tsv"),
    trinotate = file.path(demo_dir, "Trinotate_demo.tsv"),
    de = file.path(demo_dir, "DE_results_demo.tsv"),
    counts = file.path(demo_dir, "counts_demo.tsv")
  )
  missing_demo_files <- names(demo_files)[
    !vapply(demo_files, file.exists, logical(1))
  ]
  if (length(missing_demo_files)) {
    stop(
      "The packaged EchoGO demo is incomplete. Missing: ",
      paste(missing_demo_files, collapse = ", "),
      call. = FALSE
    )
  }

  results_dir <- file.path(demo_dir, "results")
  if (isTRUE(clean) && dir.exists(results_dir)) {
    unlink(results_dir, recursive = TRUE, force = TRUE)
  }

  if (!isTRUE(live_gprofiler)) {
    cache <- file.path(demo_dir, "cache")
    if (!dir.exists(cache)) {
      stop("The deterministic demo g:Profiler cache is missing: ", cache, call. = FALSE)
    }
    cache_target <- file.path(results_dir, "gprofiler")
    dir.create(cache_target, recursive = TRUE, showWarnings = FALSE)
    cache_files <- list.files(cache, recursive = TRUE, full.names = TRUE, all.files = FALSE)
    for (source in cache_files) {
      if (dir.exists(source)) next
      relative <- substring(
        normalizePath(source, winslash = "/", mustWork = TRUE),
        nchar(normalizePath(cache, winslash = "/", mustWork = TRUE)) + 2L
      )
      destination <- file.path(cache_target, relative)
      dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(source, destination, overwrite = TRUE)) {
        stop("Could not copy cached demo response: ", basename(source), call. = FALSE)
      }
    }
    message("Using deterministic cached g:Profiler demo responses (offline mode).")
  } else {
    old_force <- getOption("EchoGO.force_rerun_gprofiler", NULL)
    on.exit(options(EchoGO.force_rerun_gprofiler = old_force), add = TRUE)
    options(EchoGO.force_rerun_gprofiler = TRUE)
    message("Using live g:Profiler integration mode.")
  }

  message("Demo GOseq input: ", basename(demo_files$goseq))

  run_args <- list(
    input_dir         = demo_dir,
    goseq_file        = demo_files$goseq,
    trinotate_file    = demo_files$trinotate,
    de_file           = demo_files$de,
    count_matrix_file = demo_files$counts,
    species           = species_demo,
    target_context    = "drerio",
    context_metadata  = context_metadata_demo,
    semantic_reference_orgdb = semantic_reference_demo,
    semantic_reference_role = "target_reference",
    outdir            = results_dir,
    make_report       = isTRUE(make_report),
    strict_only       = NULL,
    run_exploratory_default_domain = isTRUE(full),
    run_rrvgo         = isTRUE(full),
    run_evaluation    = isTRUE(full),
    use_trinotate_universe = FALSE
  )
  run_args <- utils::modifyList(run_args, list(...))
  res <- do.call(run_full_echogo, run_args)

  .validate_quickstart_result(res, require_report = isTRUE(make_report))
  `%||%` <- function(a, b) if (is.null(a)) b else a
  message(
    "EchoGO demo completed successfully.\n",
    "Start with the biological HTML report: ", res$files$report_html %||% "not rendered", "\n",
    "Exact-term evidence: ", res$files$exact_term_evidence %||% res$files$consensus_xlsx, "\n",
    "Representability audit: ", res$files$mapping_table %||% file.path(results_dir, "gprofiler", "submitted_vectors", "mapping_table.csv"), "\n",
    "Results directory: ", normalizePath(results_dir, winslash = "/", mustWork = FALSE)
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
EchoGO - Input folder guide
==========================

EchoGO supports TWO input layouts.

1) Classic / de novo transcriptome layout (Trinity + Trinotate)
--------------------------------------------------------------
Put these files into input/ (or set alternative names in input/config.yml):

  - gene.counts.matrix.tsv
      - TSV: first column = feature ID, remaining columns = raw counts per sample

  - DE_*.tsv (one or more contrasts)
      - TSV differential expression results
      - must contain at least: id, log2FC, pvalue, padj (or mappable equivalents)

  - Trinotate.xls (or Trinotate .tsv/.xlsx export)
      - standard Trinotate report with GO mappings / EggNOG fields if available

2) Reference-based RNA-seq layout (DESeq2 + GOseq precomputed)
--------------------------------------------------------------
EchoGO can run directly from a reference-based experiment input folder that already contains
the GOseq enrichment output and a Trinotate-like/EggNOG annotation table.

Expected files (example from a contrast folder):

  - allcounts_table.txt
      - background count matrix (required by the current input resolver)

  - dge_<CONTRAST>.csv
      - DESeq2 results for that contrast (required by the current input resolver)

  - dge_<CONTRAST>.GOseq.enriched.tsv
      - GOseq enriched categories table (REQUIRED for reference-based mode)

  - Trinotate_for_EchoGO.tsv or <reference_label>_eggNOG_for_EchoGO.tsv
      - Trinotate-like / EggNOG-mapped annotation table (REQUIRED)

  - config.yml
      - optional but recommended (organism contexts, OrgDb, report title)

Important note about gene IDs in reference-based mode:
  - If GOseq 'gene_ids' do not match annotation 'transcript_id', EchoGO assumes 'gene_ids'
    are already valid gene names and uses them as-is (this is intended for reference-based pipelines).

To create these files from a reference annotation:
  echogo_prepare_reference_inputs(root, gff_file, protein_fasta)
Then run the generated eggNOG-mapper command and repeat the same call.

--------------------------------------------------------------
Run EchoGO
--------------------------------------------------------------
After placing files here:

  # Replace the illustrative contexts for your biological question.
  echogo_run(input_dir = 'PATH/input', outdir = 'PATH/results',
             species = c('drerio', 'omykiss'), target_context = 'drerio',
             run_rrvgo = FALSE, run_exploratory_default_domain = FALSE)

You can also point echogo_run() directly to a reference-based contrast input folder, e.g.:
  echogo_run(input_dir = '.../dge_<CONTRAST>/input',
             outdir    = '.../dge_<CONTRAST>/results',
             species = c('drerio', 'omykiss'), target_context = 'drerio',
             run_rrvgo = FALSE, run_exploratory_default_domain = FALSE)

Target GOseq is required primary target evidence. Declaring NO_TARGET below
means no target g:Profiler context is queried; it does not remove GOseq.

--------------------------------------------------------------
Scientific configuration
------------------------
EchoGO does not choose biological annotation contexts for your experiment.
Set 'species' to researcher-selected g:Profiler organism IDs, for example:
  species: [drerio, strutta]

Declare exactly one of:
  target_context: drerio      # a queried target g:Profiler context
  target_context: NO_TARGET   # no target g:Profiler context is queried

Target GOseq remains the experimental anchor when NO_TARGET is declared.
Optional context_metadata records rationale/provenance and never affects statistics.

Default-domain g:Profiler exploration is opt-in:
  run_exploratory_default_domain: false

RRvGO is optional and uses a separate semantic reference:
  run_rrvgo: false
  semantic_reference_orgdb: null   # e.g. org.Dr.eg.db when RRvGO is enabled
  semantic_reference_role: null    # target_reference or proxy

The semantic OrgDb is not an enrichment context and EchoGO does not infer it.

Full organism list (external):
  https://biit.cs.ut.ee/gprofiler/page/organism-list

Helpful EchoGO helpers:
  echogo_help()
  echogo_pick_species()
  echogo_preflight_species(c('hsapiens','mmusculus'))

Checking / installing optional RRvGO dependencies:
  echogo_list_orgdb()
  echogo_install_orgdb_instructions()
  echogo_install_orgdb('org.Dr.eg.db')
  echogo_require_orgdb(c('GO.db','org.Dr.eg.db'))
"

  writeLines(readme, file.path(input, "_README.txt"))

  yaml <- "# REQUIRED: add researcher-selected g:Profiler organism codes.
species: []

# REQUIRED: replace null with a queried context code or NO_TARGET.
target_context: null

# Optional provenance. Add context_code/context_rationale rows if useful.
context_metadata: []

# Optional broader exploratory tier; never primary evidence.
run_exploratory_default_domain: false

# Optional RRvGO semantic summarization. When true, declare both fields.
run_rrvgo: false
semantic_reference_orgdb: null
semantic_reference_role: null

report_title: EchoGO biological interpretation
"
  writeLines(yaml, file.path(input, "config.yml"))

  message("Scaffold created at: ", normalizePath(path, winslash = "/"))
  invisible(path)
}


#' Run EchoGO from a folder of inputs (and optional YAML)
#'
#' This wrapper prefers the modern `run_full_echogo(input_dir = ...)` interface.
#' It reads `input/config.yml` when present to obtain researcher-declared
#' enrichment contexts, target/no-target state, and optional RRvGO semantic
#' reference configuration.
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

  dots <- list(...)
  species <- dots$species %||% cfg$species
  if (is.null(species) || !length(unlist(species, use.names = FALSE))) {
    stop(
      "echogo_run() requires researcher-selected enrichment contexts. Set ",
      "species in config.yml or supply species = c(...).",
      call. = FALSE
    )
  }

  target_from_dots <- "target_context" %in% names(dots)
  target_from_config <- "target_context" %in% names(cfg) && !is.null(cfg$target_context)
  if (target_from_dots && is.null(dots$target_context)) {
    stop(
      "target_context = NULL is not an explicit current declaration. Use a ",
      "queried context code/label or target_context = NA_character_ for no target.",
      call. = FALSE
    )
  }
  if (!target_from_dots && !target_from_config) {
    stop(
      "echogo_run() requires a target-context declaration. Supply a queried ",
      "context code/label, or use target_context = NA_character_ in R / ",
      "target_context: NO_TARGET in config.yml.",
      call. = FALSE
    )
  }

  context_metadata <- cfg$context_metadata
  if (is.list(context_metadata) && !is.data.frame(context_metadata) && length(context_metadata)) {
    context_metadata <- dplyr::bind_rows(lapply(context_metadata, as.data.frame))
  }

  base_args <- list(
    input_dir    = input_dir,
    species      = unlist(species, use.names = TRUE),
    outdir       = outdir,
    make_report  = TRUE,
    report_title = cfg$report_title %||% NULL,
    run_rrvgo = cfg$run_rrvgo %||% TRUE,
    run_exploratory_default_domain = cfg$run_exploratory_default_domain %||% FALSE
  )
  if (target_from_config) base_args$target_context <- cfg$target_context
  if (!is.null(context_metadata) && length(context_metadata)) {
    base_args$context_metadata <- context_metadata
  }
  if (!is.null(cfg$semantic_reference_orgdb)) {
    base_args$semantic_reference_orgdb <- cfg$semantic_reference_orgdb
  }
  if (!is.null(cfg$semantic_reference_role)) {
    base_args$semantic_reference_role <- cfg$semantic_reference_role
  }
  if (!is.null(cfg$orgdb)) base_args$orgdb <- cfg$orgdb
  args <- utils::modifyList(base_args, dots)

  do.call(run_full_echogo, args)
}

#' Open or inspect the shipped EchoGO demo data
#'
#' Prints the paths to the packaged demo inputs and the v0.1.4 demo snapshot,
#' shows their contents, and, in an interactive R session, opens the folders in
#' the system file browser (Windows Explorer, macOS Finder, or Linux xdg-open).
#'
#' @export
echogo_open_demo <- function() {
  demo_in  <- echogo_demo_path()
  demo_out <- echogo_demo_results_path()

  cat("\nEchoGO demo inputs:\n  ", demo_in, "\n\n")
  print(list.files(demo_in, full.names = TRUE))

  cat("\nEchoGO v0.1.4 demo snapshot:\n  ", demo_out, "\n\n")
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

  if (interactive()) {
    cat("\nOpening demo directories...\n")
    try(open_dir(demo_in),  silent = TRUE)
    try(open_dir(demo_out), silent = TRUE)
  } else {
    cat("\nNon-interactive session: demo directories were listed but not opened.\n")
  }

  invisible(list(inputs = demo_in, results = demo_out))
}

#' Path to shipped demo inputs
#' @return A filesystem path to inst/extdata/echogo_demo
#' @export
echogo_demo_path <- function() {
  system.file("extdata", "echogo_demo", package = "EchoGO")
}

#' Path to shipped demo results (frozen v0.1.4 snapshot)
#' @return A filesystem path to the v0.1.4 product-facing demo snapshot.
#' @export
echogo_demo_results_path <- function() {
  versioned <- system.file("extdata", "echogo_demo_results_v0.1.4", package = "EchoGO")
  if (nzchar(versioned)) versioned else system.file("extdata", "echogo_demo_results", package = "EchoGO")
}

#' @title Run EchoGO pipeline for a reference-based RNA-seq scaffold
#' @description
#' Convenience wrapper that:
#'   - detects standard scaffold inputs inside \code{input_dir}
#'   - calls \code{run_echogo_pipeline()} with the right files
#'   - uses config.yml (if present) to set species/orgdb, otherwise defaults.
#'
#' Expected files inside \code{input_dir}:
#'   - *_GOseq.enriched.tsv         (required)
#'   - Trout_eggNOG_for_EchoGO.tsv  (or Trinotate.*) (required)
#'   - dge_*.csv                    (required)
#'   - allcounts_table.txt          (required)
#'   - config.yml                   (optional: species, orgdb, report_title)
#'
#' @param input_dir Path to the scaffold input/ folder for one contrast.
#' @param outdir Base output directory (EchoGO will create its standard layout here).
#' @param strict_only,run_evaluation,verbose Passed to \code{run_echogo_pipeline()}.
#' @return Invisible list returned by \code{run_echogo_pipeline()}.
#' @export
echogo_run_reference_rnaseq <- function(
    input_dir,
    outdir,
    strict_only    = FALSE,
    run_evaluation = TRUE,
    verbose        = TRUE
) {
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  files     <- list.files(input_dir, full.names = TRUE)

  # ---- Detect required files in scaffold input/ -----------------------------
  goseq_file <- grep("GOseq\\.enriched\\.", files, value = TRUE)
  if (!length(goseq_file))
    stop("No GOseq enriched file (*GOseq.enriched.*) found in: ", input_dir)
  goseq_file <- goseq_file[1]

  trinotate_file <- grep("Trout_eggNOG_for_EchoGO|Trinotate\\.", files, value = TRUE)
  if (!length(trinotate_file))
    stop("No Trinotate-like file (Trout_eggNOG_for_EchoGO*.tsv or Trinotate.*) found in: ", input_dir)
  trinotate_file <- trinotate_file[1]

  de_file <- grep("dge_.*\\.csv$", files, value = TRUE)
  if (!length(de_file))
    stop("No DESeq2 file (dge_*.csv) found in: ", input_dir)
  de_file <- de_file[1]

  count_matrix_file <- grep("allcounts_table|gene.counts.matrix", files, value = TRUE)
  if (!length(count_matrix_file))
    stop("No background count matrix (allcounts_table.txt / gene.counts.matrix.*) found in: ", input_dir)
  count_matrix_file <- count_matrix_file[1]

  # Optional config.yml for species/orgdb overrides
  config_file <- grep("config\\.yml$", files, value = TRUE)
  species_cfg <- NULL
  orgdb_cfg   <- NULL
  if (length(config_file) == 1L && requireNamespace("yaml", quietly = TRUE)) {
    cfg <- try(yaml::read_yaml(config_file), silent = TRUE)
    if (!inherits(cfg, "try-error")) {
      if (!is.null(cfg$species)) species_cfg <- unlist(cfg$species)
      if (!is.null(cfg$orgdb))   orgdb_cfg   <- cfg$orgdb
    }
  }

  # ---- Decide species / orgdb (config.yml > options > defaults) ------------
  species <- species_cfg %||%
    getOption("EchoGO.default_species",
              c("hsapiens","mmusculus","drerio"))
  orgdb   <- orgdb_cfg %||%
    getOption("EchoGO.default_orgdb", "org.Dr.eg.db")

  # ---- Run the main EchoGO pipeline ----------------------------------------
  outdir <- normalizePath(outdir, winslash = "/", mustWork = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  run_echogo_pipeline(
    goseq_file        = goseq_file,
    trinotate_file    = trinotate_file,
    de_file           = de_file,
    count_matrix_file = count_matrix_file,
    species           = species,
    orgdb             = orgdb,
    outdir            = outdir,
    strict_only       = strict_only,
    run_evaluation    = run_evaluation,
    verbose           = verbose
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

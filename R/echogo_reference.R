#' @title Run EchoGO pipeline for a reference-based RNA-seq scaffold
#' @description
#' Convenience wrapper that:
#'   - detects standard scaffold inputs inside \code{input_dir}
#'   - calls \code{run_echogo_pipeline()} with the right files
#'   - uses config.yml (if present) for explicit context and optional RRvGO
#'     semantic-reference declarations.
#'
#' Expected files inside \code{input_dir}:
#'   - *_GOseq.enriched.tsv         (required)
#'   - *_eggNOG_for_EchoGO.tsv or Trinotate.* (required)
#'   - dge_*.csv                    (required)
#'   - allcounts_table.txt          (required)
#'   - config.yml                   (recommended: contexts, target declaration,
#'                                    and optional RRvGO configuration)
#'
#' @param input_dir Path to the scaffold input/ folder for one contrast.
#' @param outdir Base output directory (EchoGO will create its standard layout here).
#' @param strict_only Deprecated compatibility alias passed to \code{run_echogo_pipeline()}.
#' @param run_exploratory_default_domain,target_context,context_metadata Passed to \code{run_echogo_pipeline()}.
#' @param species Researcher-selected g:Profiler organism contexts. If omitted,
#'   they must be present in config.yml.
#' @param run_rrvgo Enable optional RRvGO semantic reduction.
#' @param semantic_reference_orgdb Preferred RRvGO semantic-reference OrgDb.
#' @param semantic_reference_role Declared `target_reference` or `proxy` role.
#' @param orgdb Deprecated compatibility alias for `semantic_reference_orgdb`.
#' @param run_evaluation Run scoreless evaluation/diagnostic summaries.
#' @param verbose Print progress messages.
#' @return Invisible list returned by \code{run_echogo_pipeline()}.
#' @export
echogo_run_reference_rnaseq <- function(
    input_dir,
    outdir,
    strict_only    = NULL,
    run_exploratory_default_domain = FALSE,
    target_context = NULL,
    context_metadata = NULL,
    run_evaluation = TRUE,
    verbose        = TRUE,
    species = NULL,
    run_rrvgo = FALSE,
    semantic_reference_orgdb = NULL,
    semantic_reference_role = NULL,
    orgdb = NULL
) {
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  files     <- list.files(input_dir, full.names = TRUE)

  # ---- Detect required files in scaffold input/ -----------------------------
  goseq_file <- grep("GOseq\\.enriched\\.", files, value = TRUE)
  if (!length(goseq_file))
    stop("No GOseq enriched file (*GOseq.enriched.*) found in: ", input_dir)
  if (length(goseq_file) > 1L)
    stop("Multiple GOseq enriched files found. Use one contrast per input directory.")
  goseq_file <- goseq_file[1]

  trinotate_file <- grep(
    "eggNOG_for_EchoGO\\.tsv$|Trinotate\\.",
    files,
    value = TRUE,
    ignore.case = TRUE
  )
  if (!length(trinotate_file))
    stop("No Trinotate-like file (*_eggNOG_for_EchoGO.tsv or Trinotate.*) found in: ", input_dir)
  trinotate_file <- trinotate_file[1]

  de_file <- grep("dge_.*\\.csv$", files, value = TRUE)
  if (!length(de_file))
    stop("No DESeq2 file (dge_*.csv) found in: ", input_dir)
  if (length(de_file) > 1L)
    stop("Multiple DESeq2 files found. Use one contrast per input directory.")
  de_file <- de_file[1]

  count_matrix_file <- grep("allcounts_table|gene.counts.matrix", files, value = TRUE)
  if (!length(count_matrix_file))
    stop("No background count matrix (allcounts_table.txt / gene.counts.matrix.*) found in: ", input_dir)
  count_matrix_file <- count_matrix_file[1]

  # Optional config.yml for species/orgdb overrides
  config_file <- grep("config\\.yml$", files, value = TRUE)
  species_cfg <- NULL
  semantic_orgdb_cfg <- NULL
  semantic_role_cfg <- NULL
  target_cfg <- NULL
  run_rrvgo_cfg <- NULL
  if (length(config_file) == 1L && requireNamespace("yaml", quietly = TRUE)) {
    cfg <- try(yaml::read_yaml(config_file), silent = TRUE)
    if (!inherits(cfg, "try-error")) {
      if (!is.null(cfg$species)) species_cfg <- unlist(cfg$species)
      if (!is.null(cfg$semantic_reference_orgdb)) semantic_orgdb_cfg <- cfg$semantic_reference_orgdb
      if (!is.null(cfg$semantic_reference_role)) semantic_role_cfg <- cfg$semantic_reference_role
      if (!is.null(cfg$target_context)) target_cfg <- cfg$target_context
      if (!is.null(cfg$run_rrvgo)) run_rrvgo_cfg <- cfg$run_rrvgo
      if (is.null(semantic_orgdb_cfg) && !is.null(cfg$orgdb)) orgdb <- cfg$orgdb
    }
  }

  species <- species %||% species_cfg
  if (is.null(species) || !length(species)) {
    stop(
      "echogo_run_reference_rnaseq() requires researcher-selected species in ",
      "config.yml or species = c(...).",
      call. = FALSE
    )
  }
  if ((missing(target_context) || is.null(target_context)) && !is.null(target_cfg)) {
    target_context <- target_cfg
  }
  if (is.null(target_context)) {
    stop(
      "Declare target_context as a queried context or NA_character_/NO_TARGET.",
      call. = FALSE
    )
  }
  semantic_reference_orgdb <- semantic_reference_orgdb %||% semantic_orgdb_cfg
  semantic_reference_role <- semantic_reference_role %||% semantic_role_cfg
  if (!is.null(run_rrvgo_cfg)) run_rrvgo <- isTRUE(run_rrvgo_cfg)

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
    run_exploratory_default_domain = run_exploratory_default_domain,
    target_context = target_context,
    context_metadata = context_metadata,
    run_rrvgo         = run_rrvgo,
    semantic_reference_orgdb = semantic_reference_orgdb,
    semantic_reference_role = semantic_reference_role,
    run_evaluation    = run_evaluation,
    verbose           = verbose
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

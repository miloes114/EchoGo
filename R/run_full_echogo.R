#' Run the full EchoGO workflow from an input folder or explicit file paths
#'
#' Runs the canonical EchoGO pipeline end-to-end:
#' GOseq enrichment loading/annotation -> multi-context g:Profiler enrichment ->
#' scoreless exact-term evidence assembly (primary + optional default-domain exploratory layer) -> RRvGO reduction (if enabled) ->
#' GO-term networks (if enabled) -> optional evaluation and HTML report.
#'
#' Target GOseq provides the primary target evidence in the standard EchoGO
#' evidence model. Supply its precomputed results explicitly or through the
#' input-folder resolver. Declaring no target g:Profiler context does not
#' remove this GOseq experimental anchor.
#'
#' EchoGO supports two common input styles:
#'
#' \strong{1) Classic / de novo transcriptome layout} (Trinity/Trinotate-style)
#' \itemize{
#'   \item a count matrix (e.g. \code{gene.counts.matrix.tsv})
#'   \item one DE table per contrast (e.g. \code{DE_\\*.tsv})
#'   \item a Trinotate report (e.g. \code{Trinotate.xls} or \code{Trinotate\\*.tsv})
#'   \item a precomputed target GOseq enrichment table (required target evidence)
#' }
#'
#' \strong{2) Reference-based RNA-seq layout} (DESeq2 + GOseq precomputed)
#' \itemize{
#'   \item \code{allcounts_table.txt} (count matrix; required by current input resolver)
#'   \item \code{dge_CONTRAST.csv} (DESeq2 table; required by current input resolver)
#'   \item \code{CONTRAST.GOseq.enriched.tsv} (required target GOseq evidence)
#'   \item \code{Trinotate_for_EchoGO.tsv} or \code{<reference_label>_eggNOG_for_EchoGO.tsv} (annotation; required)
#' }
#'
#' In reference-based workflows, GOseq \code{gene_ids} may already be gene symbols rather
#' than transcript IDs. EchoGO preserves those values for GOseq display while the g:Profiler
#' vectors still pass through the portable canonical-name resolver.
#'
#' @param input_dir Optional. Folder containing inputs. If provided, EchoGO will attempt to
#'   auto-detect required files using filename patterns.
#' @param goseq_file Path to the required precomputed target GOseq enrichment
#'   table (e.g. \code{\\*GOseq\\*.tsv}). The argument may be omitted when
#'   \code{input_dir} enables detection of that table; the target GOseq evidence
#'   itself is required for the standard workflow.
#' @param trinotate_file Optional but typically required. Path to a Trinotate or Trinotate-like/EggNOG
#'   annotation table used to resolve portable canonical names and annotate GOseq output.
#' @param de_file Optional. Path to a differential expression table (classic mode: \code{DE_\\*.tsv};
#'   reference-based: \code{dge_\\*.csv}). If \code{input_dir} is provided and \code{de_file} is NULL,
#'   EchoGO attempts to detect it.
#' @param count_matrix_file Optional. Path to a count matrix file (classic or reference-based).
#'   If \code{input_dir} is provided and \code{count_matrix_file} is NULL, EchoGO attempts to detect it.
#' @param species Researcher-selected character vector of g:Profiler organism
#'   IDs (e.g. \code{"hsapiens"}, \code{"mmusculus"}, \code{"drerio"}).
#'   Own-data runs have no built-in biological context panel. Values are
#'   validated via \code{echogo_preflight_species()}.
#' @param species_expr Optional. A tag/taxonomy expression resolved via \code{echogo_resolve()} into
#'   g:Profiler organism IDs; merged with \code{species} when both are provided.
#' @param orgdb Deprecated compatibility alias for
#'   \code{semantic_reference_orgdb}. It is not an enrichment context.
#' @param outdir Output directory for all EchoGO results.
#' @param make_report Logical; if TRUE, render the HTML report into \code{outdir}.
#' @param report_sections Character vector selecting report sections to include.
#' @param report_top_n Integer; top-N terms to display in report tables/figures where applicable.
#' @param report_theme Character; Bootswatch theme name for the report (passed to R Markdown).
#' @param report_template Optional; path to a custom Rmd template for report rendering.
#' @param report_title Optional report title. Defaults to the output-directory name.
#' @param strict_only Deprecated compatibility alias. TRUE disables the optional
#'   default-domain exploratory layer.
#' @param run_exploratory_default_domain Logical; run the optional default-domain
#'   exploratory tier. Defaults to FALSE.
#' @param target_context Researcher-declared target g:Profiler context code or
#'   label. Use \code{NA_character_} to declare explicitly that no target
#'   g:Profiler context is queried; target GOseq remains the experimental anchor.
#' @param context_metadata Optional data frame of researcher-supplied context rationale
#'   provenance. It is not used for statistical inference or ranking.
#' @param run_rrvgo Logical; if TRUE (default), run RRvGO semantic reduction.
#' @param run_evaluation Logical; if TRUE, run evaluation/diagnostic summaries (when enabled).
#' @param semantic_reference_role Required when RRvGO is enabled: either
#'   \code{"target_reference"} or \code{"proxy"}. EchoGO does not infer it.
#' @param semantic_reference_orgdb Preferred Bioconductor OrgDb package name for
#'   RRvGO semantic similarity. Required when RRvGO is enabled and independent
#'   of the g:Profiler organism contexts.
#' @param use_trinotate_universe Deprecated compatibility flag retained in run
#'   metadata. It never redefines the experiment-derived tested universe or
#'   permits raw transcript-ID fallback.
#' @param tested_gene_ids Optional explicit tested-gene vector or file.
#' @param de_id_column,de_significant_column,de_padj_column,de_lfc_column Optional DE columns.
#' @param count_id_column,annotation_id_column Optional count/annotation ID columns.
#' @param padj_threshold,log2fc_threshold Default DE significance thresholds.
#' @param de_table_significant_only Declare that the DE table contains significant rows only.
#' @param verbose Logical; if TRUE, print progress messages.
#'
#' @return An (invisible) list returned by \code{run_echogo_pipeline()}, augmented with:
#' \itemize{
#'   \item \code{$files$species_used}: path to \code{__species_used.txt}
#'   \item \code{$files$report_html}: path to the rendered HTML report (or \code{NA_character_})
#' }
#'
#' @examples
#' \dontrun{
#' ## Scaffold workflow: populate the experimental input files first.
#' ## Required: DE results, complete tested-feature background, annotation,
#' ## and precomputed target GOseq. Replace these illustrative contexts.
#' echogo_scaffold("my_project")
#' echogo_run(
#'   input_dir = "my_project/input",
#'   outdir = "my_project/results",
#'   species = c("drerio", "omykiss"),
#'   target_context = "drerio",
#'   context_metadata = data.frame(
#'     context_code = c("drerio", "omykiss"),
#'     context_rationale = c("Target organism resource",
#'                           "Comparative teleost resource for the study question")
#'   ),
#'   run_exploratory_default_domain = FALSE,
#'   run_rrvgo = TRUE,
#'   semantic_reference_orgdb = "org.Dr.eg.db",
#'   semantic_reference_role = "target_reference"
#' )
#' ## If no queried context is the target, use target_context = NA_character_.
#' ## Target GOseq is still required. Disable RRvGO with run_rrvgo = FALSE
#' ## when no justified semantic reference is available.
#'
#' ## Reference-based RNA-seq (DESeq2 + GOseq precomputed)
#' input_dir <- "path/to/reference_project/echogo_input"
#' outdir    <- "path/to/reference_project/echogo_results"
#' fish_species <- c("drerio","strutta","gaculeatus","olatipes","trubripes","amexicanus",
#'                   "oniloticus","ssalar","omykiss","okisutch","otshawytscha")
#'
#' res <- run_full_echogo(
#'   input_dir              = input_dir,
#'   species                = fish_species,
#'   target_context         = "drerio",
#'   semantic_reference_orgdb = "org.Dr.eg.db",
#'   semantic_reference_role = "target_reference",
#'   outdir                 = outdir,
#'   run_exploratory_default_domain = FALSE,
#'   run_evaluation         = TRUE,
#'   make_report            = TRUE,
#'   verbose                = TRUE
#' )
#' }
#'
#' @seealso
#' \code{\link{echogo_help}}, \code{\link{echogo_scaffold}}, \code{\link{echogo_run}},
#' \code{\link{echogo_preflight_species}}, \code{\link{echogo_pick_species}},
#' \code{\link{echogo_resolve_reference_inputs}}
#'
#' @export


run_full_echogo <- function(
    input_dir = NULL,
    goseq_file = NULL,
    trinotate_file = NULL,
    de_file = NULL,
    count_matrix_file = NULL,
    species = NULL,
    species_expr = NULL,   # allow tag/taxonomy expressions
    orgdb   = NULL,
    outdir  = "echogo_out",
    make_report = TRUE,
    report_sections = c("overview","goseq","gprofiler","consensus","rrvgo","networks"),
    report_top_n = 25,
    report_theme = "flatly",
    report_template = NULL,
    report_title = NULL,
    strict_only = NULL,
    run_exploratory_default_domain = FALSE,
    target_context = NULL,
    context_metadata = NULL,
    run_rrvgo = TRUE,
    run_evaluation = TRUE,
    semantic_reference_role = NULL,
    use_trinotate_universe = FALSE,
    tested_gene_ids = NULL,
    de_id_column = NULL,
    de_significant_column = NULL,
    de_padj_column = NULL,
    de_lfc_column = NULL,
    count_id_column = NULL,
    annotation_id_column = NULL,
    padj_threshold = 0.05,
    log2fc_threshold = 1,
    de_table_significant_only = FALSE,
    verbose = TRUE,
    semantic_reference_orgdb = NULL
) {
  `%||%` <- function(a,b) if (!is.null(a)) a else b
  species_was_missing <- missing(species)
  target_was_missing <- missing(target_context)

  # Keep compatibility mirrors disabled unless explicitly requested.
  if (is.null(getOption("EchoGO.legacy_aliases", NULL))) {
    options(EchoGO.legacy_aliases = FALSE)
  }

  # ---- species: resolve expression to codes, then validate ------------------
  if (!is.null(species_expr) && nzchar(paste(species_expr, collapse = ""))) {
    expr_ids <- tryCatch(echogo_resolve(species_expr, refresh = FALSE),
                         error = function(e) character(0))
    if (length(expr_ids)) {
      if (species_was_missing || is.null(species) || !length(species)) {
        species <- expr_ids
        species_was_missing <- FALSE
        if (isTRUE(verbose)) message("Resolved species_expr to ", paste(species, collapse = ", "))
      } else {
        species <- unique(c(as.character(species), expr_ids))
        if (isTRUE(verbose)) message("Merged species and species_expr to ", paste(species, collapse = ", "))
      }
    } else if (isTRUE(verbose)) {
      message("species_expr did not resolve to any species; falling back to 'species' argument.")
    }
  }

  # ---- help users find & validate -------------------------------------------
  if (species_was_missing || is.null(species) || !length(species)) {
    if (interactive()) {
      message("No species provided. Opening interactive chooser...")
      try(echogo_pick_species(refresh = FALSE), silent = TRUE)
    }
  }

  species <- .echogo_resolve_enrichment_contexts(
    species,
    argument_missing = species_was_missing,
    caller = "run_full_echogo()"
  )

  species <- echogo_preflight_species(species, refresh = FALSE, error_if_unknown = TRUE)
  if (isTRUE(verbose)) message("Using species: ", paste(species, collapse = ", "))

  target_contract <- .echogo_resolve_target_declaration(
    target_context,
    argument_missing = target_was_missing,
    caller = "run_full_echogo()"
  )
  semantic_contract <- .echogo_resolve_semantic_reference(
    semantic_reference_orgdb = semantic_reference_orgdb,
    orgdb = orgdb,
    semantic_reference_role = semantic_reference_role,
    run_rrvgo = run_rrvgo,
    caller = "run_full_echogo()"
  )
  if (isTRUE(run_rrvgo)) {
    .echogo_preflight_semantic_dependencies(semantic_contract$orgdb)
  }

  # ---- resolve files if input_dir is used -----------------------------------
  if (!is.null(input_dir)) {
    find_one <- function(patterns, must = TRUE, label = "input") {
      .echogo_find_one(
        patterns = patterns,
        must = must,
        label = label,
        verbose = verbose
      )
    }
    goseq_file <- goseq_file %||% find_one(
      file.path(input_dir, c(
        "*.GOseq.enriched.tsv",
        "*GOseq*enrich*.tsv",
        "*GOseq*.tsv",
        "*.GOseq.enriched.csv",
        "*GOseq*enrich*.csv",
        "*GOseq*.csv"
      )),
      must = FALSE,
      label = "GOseq"
    )
    trinotate_file <- trinotate_file %||% find_one(
      file.path(input_dir, c(
        "Trinotate*.xls",
        "Trinotate*.xlsx",
        "Trinotate*.tsv",
        "Trinotate*.txt",
        "*eggNOG*EchoGO*.tsv",
        "*eggNOG*EchoGO*.txt"
      )),
      label = "Trinotate/annotation"
    )

    count_matrix_file <- count_matrix_file %||% find_one(
      file.path(input_dir, c("*count*matrix*", "*counts*.tsv", "*counts*.csv", "*counts*.txt")),
      label = "count matrix"
    )
    if (is.null(de_file)) {
      de_file <- find_one(
        file.path(input_dir, c("DE_*.tsv", "*DE_results*subset*", "*DE_results*.tsv", "dge_*.csv")),
        label = "differential-expression"
      )
    }
  }

  # ---- Save species provenance (one per line) ----
  outdir_norm <- normalizePath(outdir, winslash = "/", mustWork = FALSE)
  dir.create(outdir_norm, showWarnings = FALSE, recursive = TRUE)
  species_file <- file.path(outdir_norm, "__species_used.txt")
  header <- sprintf("# EchoGO species used\n# Generated: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  writeLines(c(header, species), con = species_file)

  # Expose active results root for downstream helpers (e.g., RRvGO writer hints)
  options(EchoGO.active_results_dir = outdir_norm)

  # ---- run the canonical pipeline ----
  res <- run_echogo_pipeline(
    goseq_file        = goseq_file,
    trinotate_file    = trinotate_file,
    de_file           = de_file,
    count_matrix_file = count_matrix_file,
    species           = species,
    orgdb             = NULL,
    outdir            = outdir_norm,
    strict_only       = strict_only,
    run_exploratory_default_domain = run_exploratory_default_domain,
    target_context    = target_contract$configured_value,
    context_metadata  = context_metadata,
    run_rrvgo         = run_rrvgo,
    run_evaluation    = run_evaluation,
    semantic_reference_role = semantic_contract$role,
    use_trinotate_universe = use_trinotate_universe,
    tested_gene_ids = tested_gene_ids,
    de_id_column = de_id_column,
    de_significant_column = de_significant_column,
    de_padj_column = de_padj_column,
    de_lfc_column = de_lfc_column,
    count_id_column = count_id_column,
    annotation_id_column = annotation_id_column,
    padj_threshold = padj_threshold,
    log2fc_threshold = log2fc_threshold,
    de_table_significant_only = de_table_significant_only,
    verbose           = verbose,
    semantic_reference_orgdb = semantic_contract$orgdb
  )

  # ---- render report ----
  report_html <- NA_character_
  if (isTRUE(make_report)) {
    res_dirs <- if (!is.null(res$dirs) && is.list(res$dirs)) res$dirs else list()
    report_html <- tryCatch(
      .render_echogo_report(
        report_title   = report_title %||% basename(outdir_norm),
        template       = report_template,
        params         = list(
          top_n    = report_top_n,
          sections = report_sections,
          dirs     = utils::modifyList(
            res_dirs,
            list(base = outdir_norm)
          )
        ),
        theme          = report_theme,
        outdir         = outdir_norm,
        keep_temp_rmd  = FALSE,
        verbose_render = isTRUE(verbose)
      ),
      error = function(e) {
        warning("Report rendering failed: ", conditionMessage(e))
        NA_character_
      }
    )
  }

  # --- Normalize nested compatibility output under the requested directory ----
  run_id <- basename(normalizePath(dirname(outdir_norm), winslash = "/", mustWork = FALSE))
  nested_root <- file.path(outdir_norm, run_id)
  if (dir.exists(nested_root)) {
    lift <- c("rrvgo","Similarity_based_consensus","networks","Network_analysis",
              "consensus","goseq","gprofiler","report","diagnostics","evaluation",
              "consensus_enrichment","consensus_plots_all_exploratory","consensus_plots_strict_true_consensus",
              "orthology_based_enrichment_support","cross_species_gprofiler")
    for (d in lift) {
      src <- file.path(nested_root, d)
      if (dir.exists(src)) {
        dest <- file.path(outdir_norm, d)
        dir.create(dest, recursive = TRUE, showWarnings = FALSE)
        files <- list.files(src, all.files = TRUE, full.names = TRUE, no.. = TRUE)
        if (length(files)) file.copy(files, dest, recursive = TRUE, overwrite = TRUE)
      }
    }
    try(unlink(nested_root, recursive = TRUE, force = TRUE), silent = TRUE)
  }

  # --- store paths explicitly as list fields & return
  if (is.null(res$files) || !is.list(res$files)) res$files <- list()
  res$files$species_used <- species_file
  res$files$report_html  <- if (is.character(report_html) && length(report_html) == 1 && nzchar(report_html)) {
    normalizePath(report_html, winslash = "/", mustWork = FALSE)
  } else NA_character_

  if (isTRUE(verbose)) {
    message("\nEchoGO completed. Start with the biological HTML report: ", res$files$report_html)
    message("Exact-term evidence: ", res$files$exact_term_evidence)
    message("Representability audit: ", res$files$mapping_table)
    if (length(res$files$semantic_products) == 1L && !is.na(res$files$semantic_products)) {
      message("Semantic products (if terms qualified): ", res$files$semantic_products)
    }
  }

  return(invisible(res))
}

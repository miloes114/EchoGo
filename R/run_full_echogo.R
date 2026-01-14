#' Run the full EchoGO workflow from an input folder or explicit file paths
#'
#' Runs the canonical EchoGO pipeline end-to-end:
#' GOseq enrichment loading/annotation -> cross-species g:Profiler enrichment ->
#' consensus scoring (strict + optional exploratory) -> RRvGO reduction (if enabled) ->
#' GO-term networks (if enabled) -> optional evaluation and HTML report.
#'
#' EchoGO supports two common input styles:
#'
#' \strong{1) Classic / de novo transcriptome layout} (Trinity/Trinotate-style)
#' \itemize{
#'   \item a count matrix (e.g. \code{gene.counts.matrix.tsv})
#'   \item one DE table per contrast (e.g. \code{DE_\\*.tsv})
#'   \item a Trinotate report (e.g. \code{Trinotate.xls} or \code{Trinotate\\*.tsv})
#'   \item (optional) a GOseq enrichment table if you already computed it
#' }
#'
#' \strong{2) Reference-based RNA-seq layout} (DESeq2 + GOseq precomputed)
#' \itemize{
#'   \item \code{allcounts_table.txt} (count matrix; required by current input resolver)
#'   \item \code{dge_CONTRAST.csv} (DESeq2 table; required by current input resolver)
#'   \item \code{CONTRAST.GOseq.enriched.tsv} (GOseq enriched categories; used when present)
#'   \item \code{Trinotate_for_EchoGO.tsv} or \code{Trout_eggNOG_for_EchoGO.tsv} (annotation; required)
#' }
#'
#' In reference-based workflows, GOseq \code{gene_ids} may already be gene symbols/names rather
#' than transcript IDs. If GOseq \code{gene_ids} do not match the annotation \code{transcript_id},
#' EchoGO will treat them as already-usable names (intended behavior for reference-based pipelines).
#'
#' @param input_dir Optional. Folder containing inputs. If provided, EchoGO will attempt to
#'   auto-detect required files using filename patterns.
#' @param goseq_file Optional. Path to a GOseq enrichment file (e.g. \code{\\*GOseq\\*.tsv}).
#'   If \code{input_dir} is provided, EchoGO will attempt to detect it.
#' @param trinotate_file Optional but typically required. Path to a Trinotate or Trinotate-like/EggNOG
#'   annotation table used to map gene IDs to names and derive background universe when requested.
#' @param de_file Optional. Path to a differential expression table (classic mode: \code{DE_\\*.tsv};
#'   reference-based: \code{dge_\\*.csv}). If \code{input_dir} is provided and \code{de_file} is NULL,
#'   EchoGO attempts to detect it.
#' @param count_matrix_file Optional. Path to a count matrix file (classic or reference-based).
#'   If \code{input_dir} is provided and \code{count_matrix_file} is NULL, EchoGO attempts to detect it.
#' @param species Character vector of g:Profiler organism IDs (e.g. \code{"hsapiens"}, \code{"mmusculus"},
#'   \code{"drerio"}). Validated via \code{echogo_preflight_species()}.
#' @param species_expr Optional. A tag/taxonomy expression resolved via \code{echogo_resolve()} into
#'   g:Profiler organism IDs; merged with \code{species} when both are provided.
#' @param orgdb Character scalar. Bioconductor OrgDb package name for semantic similarity steps
#'   (e.g. \code{"org.Dr.eg.db"}, \code{"org.Mm.eg.db"}).
#' @param outdir Output directory for all EchoGO results.
#' @param make_report Logical; if TRUE, render the HTML report into \code{outdir}.
#' @param report_sections Character vector selecting report sections to include.
#' @param report_top_n Integer; top-N terms to display in report tables/figures where applicable.
#' @param report_theme Character; Bootswatch theme name for the report (passed to R Markdown).
#' @param report_template Optional; path to a custom Rmd template for report rendering.
#' @param strict_only Logical; if TRUE, only compute the strict consensus stream (skip exploratory stream).
#' @param run_evaluation Logical; if TRUE, run evaluation/diagnostic summaries (when enabled).
#' @param use_trinotate_universe Logical; if TRUE, prefer a background universe derived from Trinotate/EggNOG
#'   primary names (when available) for cross-species steps (useful for reference-based workflows).
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
#' ## Classic scaffold workflow
#' echogo_scaffold("my_project")
#' echogo_run(input_dir = "my_project/input", outdir = "my_project/results")
#'
#' ## Reference-based RNA-seq (DESeq2 + GOseq precomputed)
#' input_dir <- "E:/.../dge_RRf_BBf/input"
#' outdir    <- "E:/.../dge_RRf_BBf/results"
#' fish_species <- c("drerio","strutta","gaculeatus","olatipes","trubripes","amexicanus",
#'                   "oniloticus","ssalar","omykiss","okisutch","otshawytscha")
#'
#' res <- run_full_echogo(
#'   input_dir              = input_dir,
#'   species                = fish_species,
#'   orgdb                  = "org.Dr.eg.db",
#'   outdir                 = outdir,
#'   strict_only            = FALSE,
#'   run_evaluation         = TRUE,
#'   use_trinotate_universe = TRUE,
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
    species = getOption("EchoGO.default_species", c("hsapiens","mmusculus","drerio")),
    species_expr = NULL,   # allow tag/taxonomy expressions
    orgdb   = getOption("EchoGO.default_orgdb", "org.Dr.eg.db"),
    outdir  = "echogo_out",
    make_report = TRUE,
    report_sections = c("overview","goseq","gprofiler","consensus","rrvgo","networks"),
    report_top_n = 25,
    report_theme = "flatly",
    report_template = NULL,
    strict_only = FALSE,
    run_evaluation = TRUE,
    use_trinotate_universe = FALSE,   # <-- NEW: explicit reference-mode switch
    verbose = TRUE
) {
  `%||%` <- function(a,b) if (!is.null(a)) a else b

  # --- ensure legacy mirrors are OFF by default (user can still opt-in)
  if (is.null(getOption("EchoGO.legacy_aliases", NULL))) {
    options(EchoGO.legacy_aliases = FALSE)
  }

  # ---- species: resolve expr → codes, then validate -------------------------
  if (!is.null(species_expr) && nzchar(paste(species_expr, collapse = ""))) {
    expr_ids <- tryCatch(echogo_resolve(species_expr, refresh = FALSE),
                         error = function(e) character(0))
    if (length(expr_ids)) {
      if (missing(species) || is.null(species) || !length(species)) {
        species <- expr_ids
        if (isTRUE(verbose)) message("Resolved species_expr → ", paste(species, collapse = ", "))
      } else {
        species <- unique(c(as.character(species), expr_ids))
        if (isTRUE(verbose)) message("Merged species & species_expr → ", paste(species, collapse = ", "))
      }
    } else if (isTRUE(verbose)) {
      message("species_expr did not resolve to any species; falling back to 'species' argument.")
    }
  }

  # ---- help users find & validate -------------------------------------------
  if (missing(species) || is.null(species) || !length(species)) {
    if (interactive()) {
      message("No species provided. Opening interactive chooser…")
      try(echogo_choose_species(refresh = FALSE), silent = TRUE)
    }
    stop("Please provide species=... (e.g., c('hsapiens','mmusculus')) or species_expr='...'. ",
         "Use echogo_list_species(view=TRUE) or echogo_species_lookup('human') to discover IDs.", call. = FALSE)
  }

  species <- echogo_preflight_species(species, refresh = FALSE, error_if_unknown = TRUE)
  if (isTRUE(verbose)) message("Using species: ", paste(species, collapse = ", "))

  # ---- resolve files if input_dir is used -----------------------------------
  if (!is.null(input_dir)) {
    find_one <- function(glob_vec, must = TRUE) {
      hits <- unlist(lapply(glob_vec, Sys.glob))
      if (!length(hits)) {
        if (must) stop("Could not find any of: ", paste(glob_vec, collapse = " | "), " in ", normalizePath(input_dir))
        return(NULL)
      }
      hits[[1]]
    }
    goseq_file        <- goseq_file        %||% find_one(file.path(input_dir, c("*GOseq*enrich*","*GOseq*.tsv","*GOseq*.csv")), must = FALSE)
    trinotate_file <- trinotate_file %||% find_one(
      file.path(input_dir, c(
        "Trinotate*.xls",
        "Trinotate*.xlsx",
        "Trinotate*.tsv",
        "Trinotate*.txt",
        "*eggNOG*EchoGO*.tsv",      # <- your Trout_eggNOG_for_EchoGO.tsv
        "*eggNOG*EchoGO*.txt"
      ))
    )

    count_matrix_file <- count_matrix_file %||% find_one(file.path(input_dir, c("*count*matrix*","*counts*.tsv","*counts*.csv","*counts*.txt")))
    if (is.null(de_file)) {
      de_file <- find_one(file.path(input_dir, c("DE_*.tsv","*DE_results*subset*","*DE_results*.tsv","dge_*.csv")))
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
    orgdb             = orgdb,
    outdir            = outdir_norm,
    strict_only       = strict_only,
    run_evaluation    = run_evaluation,
    use_trinotate_universe = use_trinotate_universe,  # <-- pass reference-mode flag
    verbose           = verbose
  )

  # ---- render report ----
  report_html <- NA_character_
  if (isTRUE(make_report)) {
    res_dirs <- if (!is.null(res$dirs) && is.list(res$dirs)) res$dirs else list()
    report_html <- tryCatch(
      .render_echogo_report(
        report_title   = basename(outdir_norm),
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

  # --- de-nest accidental '<outdir>/<run_id>/*' created by older helpers ----
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

  return(invisible(res))
}

run_full_echogo <- function(
    input_dir = NULL,
    goseq_file = NULL,
    trinotate_file = NULL,
    de_file = NULL,
    count_matrix_file = NULL,
    species = getOption("EchoGO.default_species", c("hsapiens","mmusculus","drerio")),
    species_expr = NULL,   # <-- NEW: allow tag/taxonomy expressions
    orgdb   = getOption("EchoGO.default_orgdb", "org.Dr.eg.db"),
    outdir  = "echogo_out",
    make_report = TRUE,
    report_sections = c("overview","goseq","gprofiler","consensus","rrvgo","networks"),
    report_top_n = 25,
    report_theme = "flatly",
    report_template = NULL,
    strict_only = FALSE,
    run_evaluation = TRUE,
    verbose = TRUE
) {
  `%||%` <- function(a,b) if (!is.null(a)) a else b

  # --- ensure legacy mirrors are OFF by default (user can still opt-in)
  if (is.null(getOption("EchoGO.legacy_aliases", NULL))) {
    options(EchoGO.legacy_aliases = FALSE)
  }

  # ---- species: resolve expr → codes, then validate -------------------------
  # If species_expr provided, resolve it and merge/override sensibly.
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
    trinotate_file    <- trinotate_file    %||% find_one(file.path(input_dir, c("Trinotate*.xls","Trinotate*.xlsx","Trinotate*.tsv","Trinotate*.txt")))
    count_matrix_file <- count_matrix_file %||% find_one(file.path(input_dir, c("*count*matrix*","*counts*.tsv","*counts*.csv")))
    if (is.null(de_file)) {
      de_file <- find_one(file.path(input_dir, c("DE_*.tsv","*DE_results*subset*","*DE_results*.tsv")))
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

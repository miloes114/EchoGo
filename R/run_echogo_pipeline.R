#' @name run_echogo_pipeline
#' @title Run the full EchoGO enrichment pipeline
#' @description Loads and processes enrichment results using GOseq and g:Profiler,
#' integrates them into a consensus table, and generates semantic, network, and
#' evaluation outputs.
#' @param goseq_file Path to GOseq enrichment results
#' @param trinotate_file Path to Trinotate annotation file
#' @param de_file Path to the full contrast-specific DE results table.
#' @param count_matrix_file Path to full count matrix
#' @param species Character vector of g:Profiler organism codes (unlimited). You may pass
#'   a named vector like c(hsapiens="human", mmusculus="mouse"); if unnamed, labels default
#'   to the codes themselves.
#' @param orgdb Character|OrgDb|list for RRvGO (optional; passed through if supported)
#' @param outdir Base output directory (standard layout will be created here)
#' @param strict_only If TRUE, skip exploratory no-background runs.
#' @param run_rrvgo Whether to run RRvGO semantic reduction (default: TRUE)
#' @param run_evaluation Whether to run enrichment evaluation (default: TRUE)
#' @param verbose Print progress
#' @param use_trinotate_universe Deprecated compatibility flag retained in run
#'   metadata. It does not intersect or redefine the tested experiment universe.
#' @param tested_gene_ids Optional explicit tested-gene vector or file.
#' @param de_id_column,de_significant_column,de_padj_column,de_lfc_column Optional DE columns.
#' @param count_id_column,annotation_id_column Optional count/annotation ID columns.
#' @param padj_threshold,log2fc_threshold Default DE significance thresholds.
#' @param de_table_significant_only Declare that the DE table contains significant rows only.
#' @return Named list with inputs, dirs, objects, files
#' @export
#' @importFrom dplyr %>% filter pull
#' @importFrom stringr str_detect
run_echogo_pipeline <- function(
    goseq_file,
    trinotate_file,
    de_file,
    count_matrix_file,
    species = getOption("EchoGO.default_species", c("hsapiens","mmusculus","drerio")),
    orgdb   = getOption("EchoGO.default_orgdb", "org.Dr.eg.db"),
    outdir  = "echogo_out",
    strict_only = FALSE,
    run_rrvgo = TRUE,
    run_evaluation = TRUE,
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
    verbose = TRUE
) {
  `%||%` <- function(a,b) if (!is.null(a)) a else b
  .mk <- function(...) { p <- file.path(...); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
  .mirror_tree <- function(src, dest) {
    if (!dir.exists(src)) return(invisible(FALSE))
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    files <- list.files(src, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
    ok <- TRUE
    for (f in files) {
      if (dir.exists(f)) next
      rel <- sub(paste0("^", gsub("\\\\","\\\\\\\\", normalizePath(src, winslash="/", mustWork=FALSE))), "", normalizePath(f, winslash="/", mustWork=FALSE))
      rel <- sub("^[/\\\\]", "", rel)
      tgt <- file.path(dest, rel)
      dir.create(dirname(tgt), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(f, tgt, overwrite = TRUE)) ok <- FALSE
    }
    cat("This legacy folder mirrors: ", basename(src), "/\n", file = file.path(dest, "__moved_to.txt"))
    invisible(ok)
  }

  # Respect user preference for legacy alias folders (default off)
  legacy_on <- isTRUE(getOption("EchoGO.legacy_aliases", FALSE))
  if (verbose && !legacy_on) message("Legacy alias folders are disabled; writing only to canonical layout.")

  # ---- Canonical standard folders under outdir ----
  dirs <- list(
    base        = .mk(outdir),
    goseq       = .mk(outdir, "goseq"),
    gprof       = .mk(outdir, "gprofiler"),
    consensus   = .mk(outdir, "consensus"),
    rrvgo       = .mk(outdir, "rrvgo"),
    network     = .mk(outdir, "networks"),
    evaluation  = .mk(outdir, "evaluation"),
    diagnostics = .mk(outdir, "diagnostics"),
    report      = .mk(outdir, "report")
  )
  # canonical consensus plot subfolders
  dirs$consensus_strict_plots <- .mk(dirs$consensus, "plots_strict")
  dirs$consensus_all_plots    <- .mk(dirs$consensus, "plots_exploratory")

  # expose active results root so relative writers (e.g., RRvGO) can anchor themselves
  options(EchoGO.active_results_dir = dirs$base)

  # ---- Legacy alias folder names (Windows-safe mirroring; only if enabled) ----
  legacy <- NULL
  if (legacy_on) {
    legacy <- list(
      orthology_based_enrichment_support    = file.path(outdir, "orthology_based_enrichment_support"),
      cross_species_gprofiler               = file.path(outdir, "cross_species_gprofiler"),
      consensus_enrichment                  = file.path(outdir, "consensus_enrichment"),
      Similarity_based_consensus            = file.path(outdir, "Similarity_based_consensus"),
      Network_analysis                      = file.path(outdir, "Network_analysis"),
      consensus_plots_strict_true_consensus = file.path(outdir, "consensus_plots_strict_true_consensus"),
      consensus_plots_all_exploratory       = file.path(outdir, "consensus_plots_all_exploratory")
    )
    invisible(lapply(legacy, dir.create, recursive = TRUE, showWarnings = FALSE))
  }

  # Build a stable species map: codes -> labels
  species_map <- if (!is.null(names(species)) && any(nzchar(names(species)))) {
    species
  } else {
    stats::setNames(as.character(species), as.character(species))
  }

  gene_sets <- prepare_gprofiler_gene_sets(
    de_results = de_file,
    tested_gene_ids = tested_gene_ids,
    annotation = trinotate_file,
    count_matrix = count_matrix_file,
    de_id_column = de_id_column,
    significant_column = de_significant_column,
    padj_column = de_padj_column,
    lfc_column = de_lfc_column,
    count_id_column = count_id_column,
    annotation_id_column = annotation_id_column,
    padj_threshold = padj_threshold,
    log2fc_threshold = log2fc_threshold,
    de_table_significant_only = de_table_significant_only,
    use_trinotate_universe = use_trinotate_universe
  )
  vectors_dir <- .mk(dirs$gprof, "submitted_vectors")
  readr::write_csv(gene_sets$mapping_table, file.path(vectors_dir, "mapping_table.csv"))
  writeLines(gene_sets$foreground_original, file.path(vectors_dir, "foreground_original.txt"))
  writeLines(gene_sets$background_original, file.path(vectors_dir, "background_original.txt"))
  writeLines(gene_sets$foreground_resolved, file.path(vectors_dir, "foreground_resolved.txt"))
  writeLines(gene_sets$background_resolved, file.path(vectors_dir, "background_resolved.txt"))
  jsonlite::write_json(
    list(
      significance_rule = gene_sets$significance_rule,
      tested_universe_source = gene_sets$tested_universe_source,
      resolver_definition = gene_sets$resolver_definition,
      counts = list(
        original_significant = length(gene_sets$foreground_original),
        original_tested = length(gene_sets$background_original),
        resolved_foreground = length(gene_sets$foreground_resolved),
        resolved_background = length(gene_sets$background_resolved),
        foreground_unmapped = length(gene_sets$foreground_unmapped),
        background_unmapped = length(gene_sets$background_unmapped)
      )
    ),
    file.path(vectors_dir, "gene_set_contract.json"),
    auto_unbox = TRUE,
    pretty = TRUE,
    na = "null"
  )

  # ---- Step 1: Load GOseq + annotations (canonical goseq/) ----
  if (verbose) message("Loading and annotating GOseq results...")
  goseq_df <- load_and_annotate_goseq(
    goseq_file        = goseq_file,
    trinotate_file    = trinotate_file,
    de_file           = de_file,
    count_matrix_file = count_matrix_file,
    output_dir        = dirs$goseq,
    gene_sets         = gene_sets
  )

  # Ensure the canonical GOseq CSV exists at the expected location
  goseq_output <- file.path(dirs$goseq, "GOseq_enrichment_full_annotated.csv")
  if (!file.exists(goseq_output)) {
    # fallback to legacy location if user reuses an old run (do not create legacy dirs)
    alt <- file.path(outdir, "orthology_based_enrichment_support", "GOseq_enrichment_full_annotated.csv")
    if (file.exists(alt)) file.copy(alt, goseq_output, overwrite = TRUE)
  }

  # GOseq lollipops
  if (exists("goseq_make_lollipops")) {
    try(goseq_make_lollipops(goseq_csv = goseq_output,
                             outdir_goseq = dirs$goseq),
        silent = TRUE)
  }
  # mirror canonical to legacy (only if enabled)
  if (legacy_on) .mirror_tree(dirs$goseq, legacy$orthology_based_enrichment_support)

  # ---- Step 2: prepare the matched portable foreground/background vectors ----
  # The DE foreground is never derived from GOseq category membership.
  de_eggnog_filtered <- gene_sets$foreground_resolved

  # The custom background is the matched tested experiment universe resolved
  # through exactly the same canonical-name procedure.
  bg_eggnog_filtered <- gene_sets$background_resolved
  if (verbose) {
    message("   - Resolved significant foreground: ", length(de_eggnog_filtered), " names.")
    message("   - Resolved tested background: ", length(bg_eggnog_filtered), " names.")
  }


  if (verbose) {
    message("Running g:Profiler across organism contexts: ", paste(names(species_map), collapse = ", "))
    message("   - resolved foreground names available for g:Profiler: ", length(de_eggnog_filtered))
  }

  # ---- Step 3: g:Profiler organism contexts (canonical gprofiler/) ----
  gprof_results <- NULL

  if (!length(de_eggnog_filtered)) {
    # No usable gene symbols; skip g:Profiler gracefully
    dir.create(dirs$gprof, showWarnings = FALSE, recursive = TRUE)
    cat(
      "g:Profiler was skipped because no portable canonical names were resolved.\n",
      "Check the mapping audit and the SwissProt, eggNOG preferred-name, or native-symbol fields.\n",
      file = file.path(dirs$gprof, "__SKIPPED_no_gene_symbols.txt")
    )
    if (verbose) message("Warning: skipping g:Profiler because no portable canonical names are available.")
  } else {
    # Skip re-run if cached CSVs already exist unless forced
    have_cached <- length(list.files(
      dirs$gprof, pattern = "^gprofiler_.*_(with_bg|nobg).*\\.csv$", recursive = TRUE
    )) > 0
    if (have_cached && !file.exists(file.path(dirs$gprof, "run_manifest.json")) &&
        !isTRUE(getOption("EchoGO.force_rerun_gprofiler", FALSE))) {
      stop(
        "Cached pre-v0.1.3 g:Profiler tables lack exact vectors and response metadata. ",
        "Use a new versioned output directory or explicitly force a complete rerun.",
        call. = FALSE
      )
    }
    if (have_cached && !isTRUE(getOption("EchoGO.force_rerun_gprofiler", FALSE))) {
      cached_queries <- list.files(
        file.path(dirs$gprof, "with_custom_background"),
        pattern = "_query\\.txt$", full.names = TRUE
      )
      cached_backgrounds <- list.files(
        file.path(dirs$gprof, "with_custom_background"),
        pattern = "_background\\.txt$", full.names = TRUE
      )
      if (!length(cached_queries) || !length(cached_backgrounds)) {
        stop("Cached g:Profiler results lack exact per-run query/background files.", call. = FALSE)
      }
      query_matches <- vapply(
        cached_queries,
        function(path) identical(readLines(path, warn = FALSE), de_eggnog_filtered),
        logical(1)
      )
      background_matches <- vapply(
        cached_backgrounds,
        function(path) identical(readLines(path, warn = FALSE), bg_eggnog_filtered),
        logical(1)
      )
      if (!all(query_matches) || !all(background_matches)) {
        stop(
          "Cached g:Profiler vectors do not match the current foreground/background gene sets. ",
          "Regenerate the cache or use a new live output directory.",
          call. = FALSE
        )
      }
    }
    if (!have_cached || isTRUE(getOption("EchoGO.force_rerun_gprofiler", FALSE))) {
      gprof_results <- run_gprofiler_cross_species(
        de_genes = de_eggnog_filtered,
        bg_genes = bg_eggnog_filtered,
        species  = species_map,
        outdir   = dirs$gprof,
        do_no_bg = !strict_only,
        significance_rule = gene_sets$significance_rule,
        resolver_definition = gene_sets$resolver_definition,
        verbose  = verbose
      )
    } else {
      if (verbose) message("   - using cached g:Profiler CSVs in: ", dirs$gprof)
    }
  }

  # g:Profiler lollipops (keep, but allow skipping via option)
  if (!isTRUE(getOption("EchoGO.skip_gprofiler_lollipops", FALSE)) && exists("gprofiler_make_lollipops")) {
    try(gprofiler_make_lollipops(gprof_base_dir = dirs$gprof,
                                 species_map    = species_map),
        silent = TRUE)
  }
  # mirror canonical to legacy (only if enabled)
  if (legacy_on) .mirror_tree(dirs$gprof, file.path(outdir, "cross_species_gprofiler"))

  # ---- Step 4: Build consensus table (canonical consensus/) ----
  if (verbose) message("Building consensus enrichment table...")
  consensus_df <- tryCatch({
    build_consensus_table(
      goseq_file    = goseq_output,
      gprofiler_dir = dirs$gprof,
      species_map   = species_map,
      output_dir    = dirs$consensus
    )
  }, error = function(e) stop("build_consensus_table() failed: ", e$message))

  if (!"consensus_score" %in% names(consensus_df)) {
    consensus_df <- .echogo_add_consensus_scores(consensus_df)
  }

  # Consensus lollipops to canonical subfolders
  if (exists("consensus_make_lollipops")) {
    try(consensus_make_lollipops(
      consensus_df = consensus_df,
      base_outdir  = dirs$consensus  # should write into plots_strict/ & plots_exploratory/
    ),
    silent = TRUE)
  }
  if (legacy_on) {
    .mirror_tree(dirs$consensus, file.path(outdir, "consensus_enrichment"))
    .mirror_tree(dirs$consensus_strict_plots, file.path(outdir, "consensus_plots_strict_true_consensus"))
    .mirror_tree(dirs$consensus_all_plots,    file.path(outdir, "consensus_plots_all_exploratory"))
  }

  # ---- Step 5: RRvGO semantic clustering (canonical rrvgo/, plus optional legacy mirror) ----
  if (isTRUE(run_rrvgo) && verbose) message("Running semantic clustering (RRVGO)...")
  rr_fun <- if (exists("run_rrvgo_consensus_analysis")) run_rrvgo_consensus_analysis else NULL
  if (isTRUE(run_rrvgo) && !is.null(rr_fun)) {
    true_consensus_df <- consensus_df %>%
      dplyr::filter(
        significant_in_any == TRUE,
        ontology %in% c("BP", "MF", "CC"),
        in_goseq == TRUE |
          origin %in% c("GO terms - g:Profiler only (with BG)", "GO terms - Consensus (with BG)")
      )

    rr_formals <- names(formals(rr_fun))
    rr_call <- list(
      df_input = true_consensus_df,
      label    = "true_consensus_with_bg"
    )
    if ("orgdb" %in% rr_formals)       rr_call$orgdb       <- orgdb
    if ("output_base" %in% rr_formals) rr_call$output_base <- dirs$rrvgo
    if ("outdir" %in% rr_formals)      rr_call$outdir      <- file.path(dirs$rrvgo, "true_consensus_with_bg")
    do.call(rr_fun, rr_call)

    extra_terms <- setdiff(
      subset(consensus_df, significant_in_any)$term_id,
      subset(consensus_df, origin %in% c("GO terms - GOseq only",
                                         "GO terms - g:Profiler only (with BG)",
                                         "GO terms - Consensus (with BG)"))$term_id
    )
    if (length(extra_terms) > 0) {
      rr2_call <- list(
        df_input = subset(consensus_df, significant_in_any),
        label    = "exploratory_all_significant"
      )
      if ("orgdb" %in% rr_formals)       rr2_call$orgdb       <- orgdb
      if ("output_base" %in% rr_formals) rr2_call$output_base <- dirs$rrvgo
      if ("outdir" %in% rr_formals)      rr2_call$outdir      <- file.path(dirs$rrvgo, "exploratory_all_significant")
      do.call(rr_fun, rr2_call)
    }
    if (legacy_on) .mirror_tree(dirs$rrvgo, file.path(outdir, "Similarity_based_consensus"))
  } else if (isTRUE(run_rrvgo)) {
    warning("run_rrvgo_consensus_analysis() not found; skipping RRvGO.")
  } else if (verbose) {
    message("Skipping RRvGO semantic clustering (run_rrvgo = FALSE).")
  }

  # ---- Step 6: GO networks (writes into canonical networks/; optional legacy mirror) ----
  if (verbose) message("Building GO term networks...")
  if (exists("run_all_networks")) {
    run_all_networks(consensus_df, outdir = dirs$base)  # keep call as-is
    if (legacy_on) {
      legacy_net <- file.path(outdir, "Network_analysis")
      if (dir.exists(legacy_net)) .mirror_tree(legacy_net, dirs$network)
      .mirror_tree(dirs$network, legacy_net)
    }
  } else {
    warning("run_all_networks() not found; skipping networks.")
  }

  # ---- Step 7: Optional evaluation (canonical evaluation/) ----
  if (run_evaluation && exists("evaluate_consensus_vs_goseq")) {
    eval_consensus_file <- file.path(dirs$consensus, "consensus_enrichment_results_with_and_without_bg.xlsx")
    eval_goseq_file     <- file.path(dirs$goseq,     "GOseq_enrichment_full_annotated.csv")
    if (!file.exists(eval_goseq_file)) {
      # try legacy location without creating it
      legacy_goseq <- file.path(outdir, "orthology_based_enrichment_support", "GOseq_enrichment_full_annotated.csv")
      if (file.exists(legacy_goseq)) eval_goseq_file <- legacy_goseq
    }
    eval_output_dir <- .mk(dirs$evaluation)
    evaluate_consensus_vs_goseq(
      consensus_file = eval_consensus_file,
      goseq_file     = eval_goseq_file,
      output_dir     = eval_output_dir
    )
  }

  # ---- return ----
  res <- list(
    inputs  = list(
      goseq_file        = goseq_file,
      trinotate_file    = trinotate_file,
      de_file           = de_file,
      count_matrix_file = count_matrix_file,
      species           = species_map,
      orgdb             = orgdb
    ),
    dirs    = dirs,
    objects = list(
      goseq_df     = goseq_df,
      consensus_df = consensus_df,
      gprofiler    = gprof_results,
      gene_sets    = gene_sets
    ),
    files   = list(
      consensus_xlsx = file.path(
        dirs$consensus,
        "consensus_enrichment_results_with_and_without_bg.xlsx"
      )
    )
  )
  invisible(res)
}


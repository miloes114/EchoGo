#' Run RRvGO Semantic Clustering on Consensus Terms (multi-OrgDb, option-aware)
#'
#' Applies RRvGO-based semantic similarity reduction to GO terms from either the strict
#' consensus set or exploratory enrichment set. Produces annotated cluster tables, bubble plots,
#' heatmaps, scatter plots, treemaps, and wordclouds per ontology.
#'
#' @param df_input A consensus enrichment data frame filtered for one mode (true consensus or exploratory).
#' @param label A label to use for the output subfolder (e.g. "true_consensus_with_bg").
#' @param output_base Directory where output will be saved (default: "similarity_based_consensus").
#'        Tip: pass a canonical path like file.path(outdir, "rrvgo") from the pipeline.
#' @param ontologies Vector of GO ontologies to process (default: c("BP", "MF", "CC")).
#' @param similarity_threshold Similarity cutoff for clustering (default: 0.7).
#' @param orgdb Character vector of OrgDb package names. If NULL or missing,
#'   falls back to getOption("EchoGO.default_orgdb", "org.Mm.eg.db").
#' @export
run_rrvgo_consensus_analysis <- function(
    df_input,
    label = "with_bg_only",
    output_base = "similarity_based_consensus",
    ontologies = c("BP", "MF", "CC"),
    orgdb = NULL,
    similarity_threshold = 0.7
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  .mk <- function(...) { p <- file.path(...); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
  .mirror_legacy <- function(src_root, legacy_root) {
    if (!nzchar(src_root) || !nzchar(legacy_root) || !dir.exists(src_root)) return(invisible(NULL))
    dir.create(legacy_root, recursive = TRUE, showWarnings = FALSE)
    src_files <- list.files(src_root, recursive = TRUE, full.names = TRUE)
    for (f in src_files) {
      if (dir.exists(f)) next
      rel <- sub(paste0("^", gsub("\\\\","\\\\\\\\", normalizePath(src_root, winslash="/", mustWork=FALSE))), "", normalizePath(f, winslash="/", mustWork=FALSE))
      rel <- sub("^[/\\\\]", "", rel)
      dest <- file.path(legacy_root, rel)
      dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
      file.copy(f, dest, overwrite = TRUE)
    }
    cat("This legacy folder mirrors: ", basename(src_root), "/\n", file = file.path(legacy_root, "__moved_to.txt"))
    invisible(NULL)
  }
  .legacy_on <- function() isTRUE(getOption("EchoGO.legacy_aliases", FALSE))

  # If output_base is relative, allow the caller to hint the active results dir:
  # options(EchoGO.active_results_dir = path_to_results)
  if (!grepl("^([A-Za-z]:)?[\\/]", output_base)) {
    hinted <- getOption("EchoGO.active_results_dir", NULL)
    if (!is.null(hinted) && nzchar(hinted)) {
      # If the legacy default name is used, canonicalize to rrvgo/ under the hinted results dir
      if (identical(output_base, "similarity_based_consensus")) {
        output_base <- file.path(hinted, "rrvgo")
      } else {
        output_base <- file.path(hinted, output_base)
      }
    }
  }

  # Root for this RRvGO run (all writes stay inside here)
  output_base <- normalizePath(.mk(output_base), winslash = "/", mustWork = FALSE)
  sub_dir     <- .mk(output_base, paste0("rrvgo_", label))

  # --- default from package option (mouse by default)
  orgdb_pkgs <- orgdb %||% getOption("EchoGO.default_orgdb", "org.Mm.eg.db")
  orgdb_pkgs <- as.character(orgdb_pkgs)
  orgdb_pkgs <- orgdb_pkgs[nzchar(orgdb_pkgs)]

  # Require GO.db (fail fast with helpful tip)
  if (!requireNamespace("GO.db", quietly = TRUE)) {
    stop("Package 'GO.db' is required. Install with: ",
         "BiocManager::install('GO.db') or EchoGO::echogo_install_orgdb('GO.db')")
  }

  # --- graceful skip of missing OrgDb
  has_pkg <- vapply(orgdb_pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (any(!has_pkg)) {
    warning("Skipping missing OrgDb: ", paste(orgdb_pkgs[!has_pkg], collapse = ", "),
            "\nTip: EchoGO::echogo_install_orgdb(c(",
            paste(sprintf("'%s'", orgdb_pkgs[!has_pkg]), collapse = ", "), "))")
  }
  orgdb_pkgs <- orgdb_pkgs[has_pkg]
  if (!length(orgdb_pkgs)) {
    message("No valid OrgDb packages available; RRvGO step skipped.")
    return(invisible(NULL))
  }

  # --- ontology guard
  ontologies <- toupper(ontologies)
  ontologies <- intersect(ontologies, c("BP","MF","CC"))
  if (!length(ontologies)) {
    message("No valid ontologies; skipping RRvGO.")
    return(invisible(NULL))
  }

  df_rrvgo <- df_input %>%
    dplyr::filter(!is.na(term_id), grepl("^GO:\\d{7}$", term_id)) %>%
    dplyr::mutate(
      go_term = trimws(term_id),
      score_source = dplyr::case_when(
        !is.na(min_pval_goseq) ~ "GOseq",
        !is.na(min_pval_gprof_bg) ~ "gprof_bg",
        !is.na(min_pval_gprof_nobg) ~ "gprof_nobg",
        TRUE ~ "unknown"
      ),
      score = -log10(pmin(min_pval_goseq, min_pval_gprof_bg, min_pval_gprof_nobg, na.rm = TRUE)),
      origin = origin
    )

  for (odb in orgdb_pkgs) {
    message("▶ RRvGO with OrgDb = ", odb)
    odb_dir <- .mk(sub_dir, paste0("OrgDb=", odb))

    for (ont in ontologies) {
      df_sub <- dplyr::filter(df_rrvgo, ontology == ont)
      scores <- df_sub$score; names(scores) <- df_sub$go_term
      scores <- scores[!is.na(scores) & is.finite(scores)]

      if (length(scores) < 2 || length(unique(scores)) < 2) {
        message("⏩ Skipping ", label, " [", odb, "]: ", ont, " — too few valid or unique scores")
        next
      }

      simMatrix <- tryCatch(
        rrvgo::calculateSimMatrix(names(scores), orgdb = odb, ont = ont, method = "Rel"),
        error = function(e) { message("⚠️ calculateSimMatrix [", odb, ":", ont, "]: ", e$message); NULL }
      )
      if (is.null(simMatrix) || nrow(simMatrix) < 2 || all(is.na(simMatrix))) {
        message("⏩ Skipping ", label, " [", odb, "]: ", ont, " — similarity matrix too sparse")
        next
      }

      reducedTerms <- tryCatch(
        rrvgo::reduceSimMatrix(simMatrix, scores, threshold = similarity_threshold, orgdb = odb),
        error = function(e) { message("⚠️ reduceSimMatrix [", odb, ":", ont, "]: ", e$message); NULL }
      )
      if (is.null(reducedTerms)) next

      reducedTerms <- reducedTerms %>%
        dplyr::mutate(go = trimws(go)) %>%
        dplyr::left_join(
          df_input %>%
            dplyr::filter(ontology == ont, grepl("^GO:\\d{7}$", term_id)) %>%
            dplyr::transmute(go_term = trimws(term_id), origin = origin) %>%
            dplyr::distinct(),
          by = c("go" = "go_term")
        ) %>%
        dplyr::mutate(origin = dplyr::coalesce(origin, "Unmatched"))

      utils::write.csv(reducedTerms, file.path(odb_dir, paste0("rrvgo_", ont, "_clusters.csv")), row.names = FALSE)

      # ---- plotting (unchanged logic; only target dir = odb_dir) ----
      tryCatch({
        pdf(file.path(odb_dir, paste0("rrvgo_", ont, "_bubbleplot.pdf")), width = 12, height = 8)
        print(
          ggplot2::ggplot(head(reducedTerms[order(reducedTerms$score, decreasing = TRUE), ], 300),
                          ggplot2::aes(x = cluster, y = score, size = size, label = term, color = origin)) +
            ggrepel::geom_text_repel(max.overlaps = 25, size = 3.5) +
            ggplot2::geom_point(alpha = 0.7) +
            ggplot2::theme_minimal(base_size = 14) +
            ggplot2::labs(title = paste("RRVGO Semantic Clusters -", ont, "[", label, "] —", odb),
                          x = "Cluster", y = "-log10(p-value)", color = "Origin")
        )
        dev.off()
      }, error = function(e) message("⚠️ Bubble plot [", odb, ":", ont, "]: ", e$message))

      tryCatch({
        pdf(file.path(odb_dir, paste0("rrvgo_", ont, "_heatmap.pdf")), width = 12, height = 10)
        simMatrix_jittered <- simMatrix + matrix(runif(length(simMatrix), -1e-6, 1e-6), nrow = nrow(simMatrix))
        rrvgo::heatmapPlot(simMatrix_jittered, reducedTerms, annotateParent = TRUE,
                           annotationLabel = "parentTerm", fontsize = 7)
        dev.off()
      }, error = function(e) message("⚠️ Heatmap [", odb, ":", ont, "]: ", e$message))

      if (nrow(simMatrix) >= 3 && any(simMatrix != 0, na.rm = TRUE)) {
        tryCatch({
          pdf(file.path(odb_dir, paste0("rrvgo_", ont, "_scatterplot.pdf")), width = 12, height = 10)
          rrvgo::scatterPlot(simMatrix, head(reducedTerms[order(reducedTerms$score, decreasing = TRUE), ], 300))
          dev.off()
        }, error = function(e) message("⚠️ Scatter [", odb, ":", ont, "]: ", e$message))
      }

      tryCatch({
        topTerms_clean <- reducedTerms %>%
          dplyr::arrange(dplyr::desc(score)) %>% head(300) %>%
          dplyr::filter(!is.na(size), is.finite(size), size > 0) %>%
          dplyr::mutate(scaled_size = log1p(size))
        pdf(file.path(odb_dir, paste0("rrvgo_", ont, "_treemap.pdf")), width = 12, height = 10)
        treemap::treemap(
          topTerms_clean, index = c("parentTerm", "term"), vSize = "scaled_size", type = "index",
          title = paste("RRVGO Treemap -", ont, "[", label, "] —", odb),
          palette = scales::hue_pal()(length(unique(topTerms_clean$parentTerm))),
          fontcolor.labels = c("#FFFFFFDD", "#00000080"), bg.labels = 0, border.col = "#00000080"
        )
        dev.off()
      }, error = function(e) message("⚠️ Treemap [", odb, ":", ont, "]: ", e$message))

      tryCatch({
        pdf(file.path(odb_dir, paste0("rrvgo_", ont, "_wordcloud.pdf")), width = 12, height = 10)
        rrvgo::wordcloudPlot(head(reducedTerms[order(reducedTerms$score, decreasing = TRUE), ], 300),
                             min.freq = 1, colors = "darkblue")
        dev.off()
      }, error = function(e) message("⚠️ Wordcloud [", odb, ":", ont, "]: ", e$message))
    }
  }

  # ---- Optional legacy mirror (keeps old scripts happy) ----
  # If you're using canonical '.../rrvgo', also mirror to sibling 'Similarity_based_consensus'.
  base_name <- basename(output_base)
  parent_dir <- dirname(output_base)
  if (.legacy_on() && tolower(base_name) %in% c("rrvgo", "similarity_based_consensus")) {
    legacy_dir <- file.path(parent_dir, "Similarity_based_consensus")
    .mirror_legacy(src_root = output_base, legacy_root = legacy_dir)
  }

  invisible(output_base)
}

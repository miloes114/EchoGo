.echogo_require_rrvgo <- function(
    is_available = requireNamespace("rrvgo", quietly = TRUE)
) {
  if (!isTRUE(is_available)) {
    stop(
      "RRvGO semantic reduction was requested, but package 'rrvgo' ",
      "is not installed.\nInstall it with:\n",
      "  BiocManager::install('rrvgo')",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.echogo_invalid_similarity_matrix <- function(x) {
  if (is.null(x)) return(TRUE)
  dims <- dim(x)
  if (length(dims) != 2L || anyNA(dims) || any(dims < 2L)) return(TRUE)
  isTRUE(all(is.na(x)))
}

# Partition primary evidence into the two Decision 0008 semantic products.
# This is deliberately a pure evidence-profile partition: it does not alter
# any source-specific values, vectors, roles, or exact-term membership.
.echogo_rrvgo_partition_inputs <- function(df_input) {
  if (!is.data.frame(df_input)) {
    stop("df_input must be a data frame", call. = FALSE)
  }
  required <- c("primary_evidence", "evidence_profile", "ontology")
  missing <- setdiff(required, names(df_input))
  if (length(missing)) {
    stop("RRvGO partition requires columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  ontologies <- c("BP", "MF", "CC")
  list(
    target_supported = dplyr::filter(
      df_input,
      .data$primary_evidence %in% TRUE,
      .data$evidence_profile %in% c("TARGET_ONLY", "TARGET_PLUS_CONTEXT"),
      .data$ontology %in% ontologies
    ),
    alternative_context_hypothesis = dplyr::filter(
      df_input,
      .data$primary_evidence %in% TRUE,
      .data$evidence_profile %in% "ALTERNATIVE_CONTEXT",
      .data$ontology %in% ontologies
    )
  )
}

.echogo_prepare_rrvgo_primary_terms <- function(df_input) {
  if ("primary_evidence" %in% names(df_input)) {
    df_input <- dplyr::filter(df_input, .data$primary_evidence %in% TRUE)
  }
  if (!"representative_order" %in% names(df_input)) {
    df_input <- .echogo_add_display_order(df_input)
  }
  df_input$.echogo_rrvgo_origin <- if ("evidence_profile" %in% names(df_input)) {
    as.character(df_input$evidence_profile)
  } else {
    "PRIMARY_EVIDENCE"
  }
  df_input %>%
    dplyr::filter(!is.na(.data$term_id), grepl("^GO:\\d{7}$", .data$term_id)) %>%
    dplyr::mutate(
      go_term = trimws(.data$term_id),
      representative_order_source = "non_inferential_representative_order",
      # RRvGO requires numeric values where larger ranks win. This is not a
      # p-value, effect size, confidence score, or integrated statistic.
      rrvgo_numeric_order = -as.numeric(.data$representative_order),
      origin = .data$.echogo_rrvgo_origin
    )
}

#' Run RRvGO semantic clustering on evidence terms
#'
#' Applies RRvGO-based semantic similarity reduction to primary custom-background
#' evidence. Representative selection uses deterministic non-inferential order,
#' never a p-value or composite score. Produces annotated cluster tables, bubble plots,
#' heatmaps, scatter plots, treemaps, and wordclouds per ontology. The exported
#' function name and default output-folder name retain the historical word
#' `consensus` for API compatibility; they do not denote consensus scoring.
#'
#' @param ... Arguments forwarded to the underlying RRvGO analysis:
#'   `df_input` (evidence data frame for one semantic product), `label` (output
#'   subfolder label), `output_base` (output directory), `ontologies` (GO
#'   ontologies to process), `orgdb` (deprecated compatibility alias),
#'   `similarity_threshold` (semantic-similarity cutoff), `semantic_product`,
#'   `semantic_reference_role`, and `semantic_reference_orgdb`. The historical
#'   score-era word `consensus` in argument values and the exported function
#'   name is retained only for API compatibility.
#' @export
run_rrvgo_consensus_analysis <- function(
    df_input,
    label = "with_bg_only",
    output_base = "similarity_based_consensus",
    ontologies = c("BP", "MF", "CC"),
    orgdb = NULL,
    similarity_threshold = 0.7,
    semantic_product = "target_supported",
    semantic_reference_role = NULL,
    semantic_reference_orgdb = NULL
) {
  semantic_contract <- .echogo_resolve_semantic_reference(
    semantic_reference_orgdb = semantic_reference_orgdb,
    orgdb = orgdb,
    semantic_reference_role = semantic_reference_role,
    run_rrvgo = TRUE,
    caller = "run_rrvgo_consensus_analysis()"
  )
  .echogo_require_rrvgo()

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

  orgdb_pkgs <- semantic_contract$orgdb

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

  semantic_product <- match.arg(
    semantic_product,
    c("target_supported", "alternative_context_hypothesis")
  )
  semantic_reference_role <- semantic_contract$role

  df_rrvgo <- .echogo_prepare_rrvgo_primary_terms(df_input)
  status_rows <- list()
  record_status <- function(ontology, status, reason = NA_character_, input_term_n = NA_integer_,
                            valid_term_n = NA_integer_, cluster_file = NA_character_) {
    status_rows[[paste(semantic_product, ontology, sep = "_")]] <<- .echogo_rrvgo_status_row(
      semantic_product = semantic_product,
      ontology = ontology,
      status = status,
      reason = reason,
      input_term_n = input_term_n,
      valid_term_n = valid_term_n,
      semantic_reference_orgdb = semantic_contract$orgdb,
      semantic_reference_role = semantic_reference_role,
      semantic_method = "Rel",
      cluster_file = cluster_file
    )
  }

  for (odb in orgdb_pkgs) {
    message("RRvGO with OrgDb = ", odb)
    odb_dir <- .mk(sub_dir, paste0("OrgDb=", odb))

    for (ont in ontologies) {
      df_sub <- dplyr::filter(df_rrvgo, ontology == ont)
      scores <- df_sub$rrvgo_numeric_order; names(scores) <- df_sub$go_term
      scores <- scores[!is.na(scores) & is.finite(scores)]

      if (length(scores) < 2 || length(unique(scores)) < 2) {
        message("Skipping ", label, " [", odb, "]: ", ont, " - too few valid or unique scores")
        record_status(ont, "SKIPPED_TOO_FEW_TERMS",
                      reason = "Fewer than two valid terms or unique representative orders were available.",
                      input_term_n = nrow(df_sub), valid_term_n = length(scores))
        next
      }

      similarity_error <- NULL
      simMatrix <- tryCatch(
        rrvgo::calculateSimMatrix(names(scores), orgdb = odb, ont = ont, method = "Rel"),
        error = function(e) { similarity_error <<- e$message; message("Warning: calculateSimMatrix [", odb, ":", ont, "]: ", e$message); NULL }
      )
      if (!is.null(similarity_error)) {
        record_status(ont, "ERROR_SIMILARITY", reason = similarity_error,
                      input_term_n = nrow(df_sub), valid_term_n = length(scores))
        next
      }
      if (.echogo_invalid_similarity_matrix(simMatrix)) {
        message("Skipping ", label, " [", odb, "]: ", ont, " - similarity matrix too sparse")
        record_status(ont, "SKIPPED_SPARSE_SIMILARITY",
                      reason = "The semantic similarity matrix was unavailable or too sparse.",
                      input_term_n = nrow(df_sub), valid_term_n = length(scores))
        next
      }

      reduction_error <- NULL
      reducedTerms <- tryCatch(
        rrvgo::reduceSimMatrix(simMatrix, scores, threshold = similarity_threshold, orgdb = odb),
        error = function(e) { reduction_error <<- e$message; message("Warning: reduceSimMatrix [", odb, ":", ont, "]: ", e$message); NULL }
      )
      if (!is.null(reduction_error)) {
        record_status(ont, "ERROR_REDUCTION", reason = reduction_error,
                      input_term_n = nrow(df_sub), valid_term_n = length(scores))
        next
      }
      if (is.null(reducedTerms)) {
        record_status(ont, "ERROR_REDUCTION", reason = "RRvGO returned no reduced term table.",
                      input_term_n = nrow(df_sub), valid_term_n = length(scores))
        next
      }

      reducedTerms <- reducedTerms %>%
        dplyr::mutate(go = trimws(go)) %>%
        dplyr::left_join(
          df_sub %>%
            dplyr::transmute(go_term = .data$go_term, origin = .data$.echogo_rrvgo_origin) %>%
            dplyr::distinct(),
          by = c("go" = "go_term")
        ) %>%
        dplyr::mutate(
          origin = dplyr::coalesce(origin, "Unmatched"),
          semantic_method = "Rel",
          semantic_reference_orgdb = odb,
          semantic_reference_role = semantic_reference_role,
          semantic_reference_configuration_source = semantic_contract$source,
          semantic_product = semantic_product
        )

      cluster_path <- file.path(odb_dir, paste0("rrvgo_", ont, "_clusters.csv"))
      utils::write.csv(reducedTerms, cluster_path, row.names = FALSE)
      record_status(ont, "GENERATED", input_term_n = nrow(df_sub), valid_term_n = length(scores),
                    cluster_file = normalizePath(cluster_path, winslash = "/", mustWork = FALSE))

      # ---- plotting (unchanged logic; only target dir = odb_dir) ----
      tryCatch({
        pdf(file.path(odb_dir, paste0("rrvgo_", ont, "_bubbleplot.pdf")), width = 12, height = 8)
        print(
          ggplot2::ggplot(head(reducedTerms[order(reducedTerms$score, decreasing = TRUE), ], 300),
                          ggplot2::aes(x = cluster, y = score, size = size, label = term, color = origin)) +
            ggrepel::geom_text_repel(max.overlaps = 25, size = 3.5) +
            ggplot2::geom_point(alpha = 0.7) +
            ggplot2::theme_minimal(base_size = 14) +
            ggplot2::labs(title = paste("RRVGO Semantic Clusters -", ont, "[", label, "] -", odb),
                          x = "Cluster", y = "Representative order (non-inferential)", color = "Evidence profile")
        )
        dev.off()
      }, error = function(e) message("Warning: bubble plot [", odb, ":", ont, "]: ", e$message))

      tryCatch({
        pdf(file.path(odb_dir, paste0("rrvgo_", ont, "_heatmap.pdf")), width = 12, height = 10)
        simMatrix_jittered <- simMatrix + matrix(runif(length(simMatrix), -1e-6, 1e-6), nrow = nrow(simMatrix))
        rrvgo::heatmapPlot(simMatrix_jittered, reducedTerms, annotateParent = TRUE,
                           annotationLabel = "parentTerm", fontsize = 7)
        dev.off()
      }, error = function(e) message("Warning: heatmap [", odb, ":", ont, "]: ", e$message))

      if (nrow(simMatrix) >= 3 && any(simMatrix != 0, na.rm = TRUE)) {
        tryCatch({
          pdf(file.path(odb_dir, paste0("rrvgo_", ont, "_scatterplot.pdf")), width = 12, height = 10)
          rrvgo::scatterPlot(simMatrix, head(reducedTerms[order(reducedTerms$score, decreasing = TRUE), ], 300))
          dev.off()
        }, error = function(e) message("Warning: scatter [", odb, ":", ont, "]: ", e$message))
      }

      tryCatch({
        topTerms_clean <- reducedTerms %>%
          dplyr::arrange(dplyr::desc(score)) %>% head(300) %>%
          dplyr::filter(!is.na(size), is.finite(size), size > 0) %>%
          dplyr::mutate(scaled_size = log1p(size))
        pdf(file.path(odb_dir, paste0("rrvgo_", ont, "_treemap.pdf")), width = 12, height = 10)
        treemap::treemap(
          topTerms_clean, index = c("parentTerm", "term"), vSize = "scaled_size", type = "index",
          title = paste("RRVGO Treemap -", ont, "[", label, "] -", odb),
          palette = scales::hue_pal()(length(unique(topTerms_clean$parentTerm))),
          fontcolor.labels = c("#FFFFFFDD", "#00000080"), bg.labels = 0, border.col = "#00000080"
        )
        dev.off()
      }, error = function(e) message("Warning: treemap [", odb, ":", ont, "]: ", e$message))

      tryCatch({
        pdf(file.path(odb_dir, paste0("rrvgo_", ont, "_wordcloud.pdf")), width = 12, height = 10)
        rrvgo::wordcloudPlot(head(reducedTerms[order(reducedTerms$score, decreasing = TRUE), ], 300),
                             min.freq = 1, colors = "darkblue")
        dev.off()
      }, error = function(e) message("Warning: wordcloud [", odb, ":", ont, "]: ", e$message))
    }
  }

  if (exists(".echogo_upsert_rrvgo_status", mode = "function") && length(status_rows)) {
    .echogo_upsert_rrvgo_status(status_rows, output_base)
  }

  # Optionally mirror canonical output to the legacy directory name.
  base_name <- basename(output_base)
  parent_dir <- dirname(output_base)
  if (.legacy_on() && tolower(base_name) %in% c("rrvgo", "similarity_based_consensus")) {
    legacy_dir <- file.path(parent_dir, "Similarity_based_consensus")
    .mirror_legacy(src_root = output_base, legacy_root = legacy_dir)
  }

  invisible(output_base)
}

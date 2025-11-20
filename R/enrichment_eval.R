#' @name evaluate_consensus_vs_goseq
#' @title Evaluate Consensus vs. GOseq Enrichment
#' @description
#' Compares EchoGO consensus enrichment results against GOseq-only enrichment
#' using multiple complementary metrics:
#' \itemize{
#'   \item Jaccard Index for term overlap (overall, per ontology, by depth)
#'   \item Enrichment Quality Index (EQI) distributions and summary statistics
#'   \item Fold Enrichment distributions
#'   \item Term origin distributions
#'   \item Rarefaction curves showing recovery of new terms from multi-species integration
#'   \item Network complexity comparison between strict (True Consensus) and exploratory modes
#' }
#' Generates Venn diagrams, density plots, cumulative and permutation rarefaction curves,
#' and summary tables. Designed to mirror the functionality and output structure of the
#' corresponding .Rmd analysis step in the EchoGO workflow.
#'
#' @param consensus_file Path to consensus enrichment Excel file.
#' @param goseq_file Path to GOseq enrichment results CSV.
#' @param output_dir Directory to store evaluation outputs
#'   (recommended canonical: "<run>/evaluation"; legacy accepted).
#' @param network_dir Directory where network outputs live. If omitted, this function
#'   prefers "<run>/networks" (canonical) and falls back to "<run>/Network_analysis" (legacy).
#' @return Invisibly returns NULL; writes plots and tables to `output_dir`.
#' @export
#' @importFrom dplyr %>% filter pull mutate select group_by summarise bind_rows arrange desc n_distinct across all_of rename_with transmute slice_head
#' @importFrom stringr str_detect str_remove_all
#' @importFrom ggplot2 ggplot aes geom_density facet_wrap labs scale_fill_manual theme_minimal ggsave geom_line geom_point geom_ribbon theme element_text
#' @importFrom openxlsx read.xlsx write.xlsx
#' @importFrom readr read_csv write_csv
#' @importFrom ggvenn ggvenn
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom igraph read_graph degree vcount ecount
evaluate_consensus_vs_goseq <- function(
    consensus_file,
    goseq_file,
    output_dir,
    network_dir = file.path(dirname(dirname(normalizePath(consensus_file, winslash = "/", mustWork = FALSE))),
                            "Network_analysis")   # legacy default kept for API compatibility
) {
  # ---- Basic I/O guards ----
  if (!file.exists(consensus_file)) stop("evaluate_consensus_vs_goseq(): consensus_file not found: ", consensus_file)
  if (!file.exists(goseq_file))     stop("evaluate_consensus_vs_goseq(): goseq_file not found: ", goseq_file)

  # Canonical vs legacy resolution (prefer 'networks/', tolerate 'Network_analysis/')
  if (missing(network_dir) || is.null(network_dir)) {
    run_root <- dirname(dirname(normalizePath(consensus_file, winslash = "/", mustWork = FALSE)))
    candidates <- c(file.path(run_root, "networks"),
                    file.path(run_root, "Network_analysis"))
    hit <- candidates[dir.exists(candidates)][1]
    network_dir <- if (length(hit)) hit else candidates[1]
  } else {
    # If a legacy path was synthesized by the default, try to upgrade to canonical if present
    run_root <- dirname(dirname(normalizePath(consensus_file, winslash = "/", mustWork = FALSE)))
    canon_try <- file.path(run_root, "networks")
    if (dir.exists(canon_try)) network_dir <- canon_try
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  exploratory_dir <- file.path(output_dir, "exploratory_no_bg")
  dir.create(exploratory_dir, recursive = TRUE, showWarnings = FALSE)

  consensus_all <- openxlsx::read.xlsx(consensus_file) %>%
    dplyr::filter(.data$ontology %in% c("BP", "MF", "CC"),
                  .data$significant_in_any == TRUE)

  goseq <- readr::read_csv(goseq_file, show_col_types = FALSE) %>%
    dplyr::filter(!is.na(.data$term),
                  .data$over_represented_FDR <= 0.05) %>%
    dplyr::mutate(term_id = .data$clean_go_term,
                  label   = .data$term)

  # Step 1: Venn diagrams + Jaccard comparison
  compare_top_terms(consensus_all, goseq, output_dir,      mode = "with_bg")
  compare_top_terms(consensus_all, goseq, exploratory_dir, mode = "no_bg")

  # Step 2: EQI distributions + summary
  summarize_eqi_distributions(consensus_all, goseq, output_dir,      mode = "with_bg")
  summarize_eqi_distributions(consensus_all, goseq, exploratory_dir, mode = "no_bg")

  # Step 3: Rarefaction (cumulative + permutation)
  run_rarefaction_curves(consensus_file, goseq_file, output_dir)

  # Step 4: Network complexity (read from canonical/legacy; write comparison to both)
  # helper to normalize suffixes from summarize_network_complexity()
  rename_from_label <- function(df, label, alias) {
    pat <- paste0("_", label, "$")
    dplyr::rename_with(df, ~ sub(pat, paste0("_", alias), .x))
  }

  net_bg  <- summarize_network_complexity("with_bg",           network_dir = network_dir) %>%
    rename_from_label("with_bg", "bg")
  net_exp <- summarize_network_complexity("with_bg_and_nobg",  network_dir = network_dir) %>%
    rename_from_label("with_bg_and_nobg", "exploratory")

  comp <- dplyr::full_join(net_bg, net_exp, by = "ontology") %>%
    dplyr::mutate(
      delta_terms      = .data$total_terms_exploratory - .data$total_terms_bg,
      delta_edges      = .data$total_edges_exploratory - .data$total_edges_bg,
      delta_avg_degree = .data$avg_degree_exploratory - .data$avg_degree_bg
    ) %>%
    dplyr::select(
      ontology,
      total_terms_bg, total_edges_bg, avg_degree_bg,
      total_terms_exploratory, total_edges_exploratory, avg_degree_exploratory,
      delta_terms, delta_edges, delta_avg_degree
    )

  # Write inside evaluation folder; also drop a canonical copy into <run>/networks
  readr::write_csv(comp, file.path(output_dir, "network_complexity_comparison.csv"))

  # Canonical target
  net_canon <- if (grepl("[/\\\\]networks$", network_dir)) network_dir
  else file.path(dirname(network_dir), "networks")
  dir.create(net_canon, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(comp, file.path(net_canon, "network_complexity_comparison.csv"))

  invisible(NULL)
}  # end evaluate_consensus_vs_goseq()


#' Compute Jaccard index between two sets
#' @keywords internal
jaccard_index <- function(set1, set2) {
  inter <- length(intersect(set1, set2))
  uni   <- length(union(set1, set2))
  if (uni == 0) return(NA_real_)
  round(inter / uni, 3)
}

#' Compare top N terms and draw Venn diagrams + Jaccard CSVs
#' Mirrors Rmd "venn_top50_*" + Jaccard summaries (with_bg / no_bg)
#' @keywords internal
compare_top_terms <- function(consensus_all, goseq, output_dir, mode = "with_bg", top_n = 50) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  score_col <- if (identical(mode, "with_bg")) "consensus_score" else "consensus_score_all"
  if (!score_col %in% names(consensus_all)) {
    warning("compare_top_terms(): missing ", score_col, " in consensus table; skipping.")
    return(invisible(NULL))
  }

  # choose consensus slice for strict vs exploratory (and build robust label/id)
  consensus_use <- if (identical(mode, "with_bg")) {
    consensus_all %>%
      dplyr::filter(.data$origin %in% c("GO terms - Consensus (with BG)")) %>%
      dplyr::mutate(score = .data[[score_col]])
  } else {
    consensus_all %>%
      dplyr::mutate(score = .data[[score_col]])
  } %>%
    dplyr::mutate(
      plot_label = dplyr::coalesce(.data$term_name, .data$term_id),
      plot_id    = .data$term_id
    )

  # Defensive column checks for GOseq
  req_cols <- c("clean_go_term","over_represented_FDR","foldEnrichment","ontology","term")
  if (!all(req_cols %in% names(goseq))) {
    warning("compare_top_terms(): GOseq is missing required columns: ",
            paste(setdiff(req_cols, names(goseq)), collapse = ", "),
            ". Skipping.")
    return(invisible(NULL))
  }

  depth_vec <- if ("depth" %in% names(goseq)) unique(stats::na.omit(as.numeric(goseq$depth))) else numeric(0)

  pick_top <- function(df, score_var, n = top_n) {
    df %>% dplyr::arrange(dplyr::desc(.data[[score_var]])) %>% dplyr::slice_head(n = n)
  }

  # Overall
  g_top_all <- pick_top(goseq, "foldEnrichment")
  c_top_all <- pick_top(consensus_use, "score")

  # Save Venn (overall) — only if both sides non-empty
  jacc_rows <- list()
  if (nrow(g_top_all) > 0 && nrow(c_top_all) > 0) {
    if (requireNamespace("ggvenn", quietly = TRUE)) {
      p_venn <- ggvenn::ggvenn(
        list(GOseq = g_top_all$term, Consensus = c_top_all$plot_label),
        fill_color = c("steelblue","orange"), stroke_size = 0.5, set_name_size = 5
      )
      ggplot2::ggsave(file.path(output_dir, "venn_top50_all.pdf"), p_venn, width = 7, height = 6)
    } else {
      message("ggvenn not installed; skipping Venn plot (overall).")
    }

    jacc_rows[[length(jacc_rows) + 1]] <- tibble::tibble(
      mode = mode, comparison = "top50_all",
      jaccard = jaccard_index(g_top_all$clean_go_term, c_top_all$plot_id),
      n_goseq = dplyr::n_distinct(g_top_all$clean_go_term),
      n_consensus = dplyr::n_distinct(c_top_all$plot_id)
    )
  }

  # By ontology (all significant)
  for (ont in c("BP","MF","CC")) {
    g_ont <- goseq %>% dplyr::filter(.data$ontology == ont, .data$over_represented_FDR <= 0.05)
    c_ont <- consensus_use %>% dplyr::filter(.data$ontology == ont, .data$significant_in_any == TRUE)

    if (nrow(g_ont) > 0 && nrow(c_ont) > 0) {
      if (requireNamespace("ggvenn", quietly = TRUE)) {
        p_venn_ont <- ggvenn::ggvenn(
          list(GOseq = g_ont$term, Consensus = c_ont$plot_label),
          fill_color = c("steelblue","orange"), stroke_size = 0.5, set_name_size = 5
        )
        ggplot2::ggsave(file.path(output_dir, paste0("venn_all_significant_", ont, ".pdf")), p_venn_ont, width = 7, height = 6)
      } else {
        message("ggvenn not installed; skipping Venn plot for ontology: ", ont)
      }

      jacc_rows[[length(jacc_rows) + 1]] <- tibble::tibble(
        mode = mode,
        comparison = paste0("all_significant_", ont),
        jaccard = jaccard_index(g_ont$clean_go_term, c_ont$plot_id),
        n_goseq = dplyr::n_distinct(g_ont$clean_go_term),
        n_consensus = dplyr::n_distinct(c_ont$plot_id)
      )
    }
  }

  # Depth-specific top-50 venns/Jaccard (if depth exists)
  if (length(depth_vec)) {
    for (d in depth_vec) {
      g_d <- goseq        %>% dplyr::filter(.data$depth == d) %>% pick_top("foldEnrichment")
      c_d <- consensus_use %>% dplyr::filter(.data$depth == d) %>% pick_top("score")

      if (nrow(g_d) && nrow(c_d)) {
        if (requireNamespace("ggvenn", quietly = TRUE)) {
          p_venn_d <- ggvenn::ggvenn(
            list(GOseq = g_d$term, Consensus = c_d$plot_label),
            fill_color = c("steelblue","orange"), stroke_size = 0.5, set_name_size = 5
          )
          ggplot2::ggsave(file.path(output_dir, paste0("venn_top50_depth", d, ".pdf")), p_venn_d, width = 7, height = 6)
        } else {
          message("ggvenn not installed; skipping Venn plot for depth: ", d)
        }

        jacc_rows[[length(jacc_rows) + 1]] <- tibble::tibble(
          mode = mode,
          comparison = paste0("top50_depth", d),
          jaccard = jaccard_index(g_d$clean_go_term, c_d$plot_id),
          n_goseq = dplyr::n_distinct(g_d$clean_go_term),
          n_consensus = dplyr::n_distinct(c_d$plot_id)
        )
      }
    }
  }

  if (length(jacc_rows)) {
    jacc_tab <- dplyr::bind_rows(jacc_rows)
    readr::write_csv(jacc_tab, file.path(output_dir, paste0("top", top_n, "_jaccard_summary_", mode, ".csv")))
  }

  invisible(NULL)
}

#' Summarize EQI distributions and enrichment strength (PDF + Excel)
#' Mirrors Rmd "EQI_distribution_by_ontology.pdf" and FE distribution
#' @keywords internal
summarize_eqi_distributions <- function(consensus_all, goseq, output_dir, mode = "with_bg") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  score_col <- if (identical(mode, "with_bg")) "consensus_score" else "consensus_score_all"
  if (!score_col %in% names(consensus_all)) {
    warning("summarize_eqi_distributions(): missing ", score_col, "; skipping.")
    return(invisible(NULL))
  }

  safe_mlog10 <- function(p) -log10(pmax(p, .Machine$double.xmin))

  cons_use <- if (identical(mode, "with_bg")) {
    consensus_all[consensus_all$origin %in% "GO terms - Consensus (with BG)", , drop = FALSE]
  } else consensus_all

  cons_use$score <- cons_use[[score_col]]
  cons_use$score[is.na(cons_use$score)] <- 0

  cons_use$min_pval_any <- pmin(
    ifelse(is.na(cons_use$min_pval_goseq),      1, cons_use$min_pval_goseq),
    ifelse(is.na(cons_use$min_pval_gprof_bg),   1, cons_use$min_pval_gprof_bg),
    ifelse(is.na(cons_use$min_pval_gprof_nobg), 1, cons_use$min_pval_gprof_nobg),
    na.rm = TRUE
  )

  cons_use$enrichment <- pmax(
    ifelse(is.na(cons_use$fold_enrichment_goseq), 0, cons_use$fold_enrichment_goseq),
    ifelse(is.na(cons_use$avg_fold_gprof_bg),     0, cons_use$avg_fold_gprof_bg),
    ifelse(is.na(cons_use$avg_fold_gprof_nobg),   0, cons_use$avg_fold_gprof_nobg),
    na.rm = TRUE
  )

  cons_use$EQI <- log2(1 + cons_use$score) * safe_mlog10(cons_use$min_pval_any)

  goseq_plot <- goseq %>%
    dplyr::filter(ontology %in% c("BP","MF","CC"),
                  over_represented_FDR <= 0.05) %>%
    dplyr::transmute(
      term_id    = clean_go_term,
      ontology   = ontology,
      method     = "GOseq",
      EQI        = log2(1 + foldEnrichment) * safe_mlog10(over_represented_FDR),
      enrichment = foldEnrichment
    )

  cons_plot <- cons_use %>%
    dplyr::mutate(method = "Consensus") %>%
    dplyr::select(term_id, ontology, method, EQI, enrichment)

  combined <- dplyr::bind_rows(goseq_plot, cons_plot)

  p_eqi <- ggplot2::ggplot(combined, ggplot2::aes(x = EQI, fill = method)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::facet_wrap(~ontology, scales = "free") +
    ggplot2::labs(title = "EQI Distributions by Method and Ontology",
                  x = "Enrichment Quality Index", y = "Density") +
    ggplot2::scale_fill_manual(values = c("GOseq" = "steelblue", "Consensus" = "darkgreen")) +
    ggplot2::theme_minimal(base_size = 13)

  p_fe <- ggplot2::ggplot(combined, ggplot2::aes(x = enrichment, fill = method)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::facet_wrap(~ontology, scales = "free") +
    ggplot2::labs(title = "Fold Enrichment Distributions by Method and Ontology",
                  x = "Fold Enrichment", y = "Density") +
    ggplot2::scale_fill_manual(values = c("GOseq" = "steelblue", "Consensus" = "darkgreen")) +
    ggplot2::theme_minimal(base_size = 13)

  ggplot2::ggsave(file.path(output_dir, "EQI_distribution_by_ontology.pdf"),       p_eqi, width = 10, height = 6)
  ggplot2::ggsave(file.path(output_dir, "FoldEnrichment_distribution_by_ontology.pdf"), p_fe, width = 10, height = 6)

  eq_summary <- combined %>%
    dplyr::group_by(ontology, method) %>%
    dplyr::summarise(
      n                 = dplyr::n(),
      median_EQI        = stats::median(EQI, na.rm = TRUE),
      mean_EQI          = mean(EQI, na.rm = TRUE),
      median_enrichment = stats::median(enrichment, na.rm = TRUE),
      .groups = "drop"
    )
  openxlsx::write.xlsx(eq_summary, file.path(output_dir, "summary_consensus_vs_goseq.xlsx"), rowNames = FALSE)
  invisible(NULL)
}

#' Rarefaction curves (true consensus vs exploratory) + cumulative tables
#' Mirrors the Rmd chunk "rarefaction_binary_matrix" + permutation rarefaction
#' @keywords internal
run_rarefaction_curves <- function(consensus_file, goseq_file, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  comparison_dir <- output_dir
  suppressWarnings(dir.create(file.path(comparison_dir, "exploratory_no_bg"), recursive = TRUE, showWarnings = FALSE))

  consensus <- openxlsx::read.xlsx(consensus_file) %>%
    dplyr::filter(.data$significant_in_any == TRUE, .data$ontology %in% c("BP","MF","CC"))

  species_cols <- grep("^in_", names(consensus), value = TRUE)
  species_cols <- setdiff(species_cols, "in_goseq")

  # presence/absence transform
  consensus_pa <- consensus %>%
    dplyr::select(.data$term_id, .data$ontology, dplyr::all_of(species_cols), .data$in_goseq) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(species_cols), ~ as.integer(. == TRUE), .names = "pa_{.col}")) %>%
    dplyr::mutate(pa_in_goseq = as.integer(.data$in_goseq == TRUE))

  pa_species_cols <- grep("^pa_in_", names(consensus_pa), value = TRUE)
  consensus_pa[pa_species_cols] <- lapply(consensus_pa[pa_species_cols], function(x) ifelse(is.na(x), 0, x))

  species_info <- tibble::tibble(col = pa_species_cols) %>%
    dplyr::mutate(
      species      = stringr::str_remove_all(.data$col, "^pa_in_|_nobg$"),
      nobg         = stringr::str_detect(.data$col, "_nobg$"),
      species_mode = .data$col
    ) %>%
    dplyr::filter(.data$species != "goseq")

  build_mode <- function(bg_mode = TRUE) {
    tag <- if (bg_mode) "true_consensus" else "exploratory"
    rarefaction_all <- list()

    # Per-ontology cumulative (ordered by contribution)
    for (ont in c("BP","MF","CC")) {
      base_terms <- consensus_pa %>%
        dplyr::filter(.data$ontology == ont, .data$pa_in_goseq == 1) %>%
        dplyr::pull(.data$term_id) %>% unique()

      entries <- species_info %>% dplyr::filter(.data$nobg == !bg_mode)
      term_sets <- purrr::map(entries$species_mode, ~ {
        colname <- .
        consensus_pa %>%
          dplyr::filter(.data$ontology == ont, .data[[colname]] == 1, .data$pa_in_goseq == 0) %>%
          dplyr::pull(.data$term_id) %>% unique()
      })
      names(term_sets) <- entries$species

      species_order <- names(sort(sapply(term_sets, length)))

      cumulative <- base_terms
      curve <- tibble::tibble(
        species = "GOseq", ontology = ont,
        n_terms = length(cumulative), n_new_terms = length(base_terms),
        step = 0L
      )
      for (i in seq_along(species_order)) {
        sp <- species_order[i]
        new_terms <- setdiff(term_sets[[sp]], cumulative)
        cumulative <- union(cumulative, new_terms)
        curve <- dplyr::bind_rows(curve, tibble::tibble(
          species = sp, ontology = ont,
          n_terms = length(cumulative), n_new_terms = length(new_terms),
          step = i
        ))
      }
      rarefaction_all[[ont]] <- curve
    }

    # Combined ALL
    base_combined <- consensus_pa %>%
      dplyr::filter(.data$pa_in_goseq == 1) %>%
      dplyr::pull(.data$term_id) %>% unique()

    entries <- species_info %>% dplyr::filter(.data$nobg == !bg_mode)
    term_sets_all <- purrr::map(entries$species_mode, ~ {
      colname <- .
      consensus_pa %>%
        dplyr::filter(.data[[colname]] == 1, .data$pa_in_goseq == 0) %>%
        dplyr::pull(.data$term_id) %>% unique()
    })
    names(term_sets_all) <- entries$species

    species_order <- names(sort(sapply(term_sets_all, length)))

    cumulative <- base_combined
    curve_all <- tibble::tibble(
      species = "GOseq", ontology = "ALL",
      n_terms = length(cumulative), n_new_terms = length(base_combined),
      step = 0L
    )
    for (i in seq_along(species_order)) {
      sp <- species_order[i]
      new_terms <- setdiff(term_sets_all[[sp]], cumulative)
      cumulative <- union(cumulative, new_terms)
      curve_all <- dplyr::bind_rows(curve_all, tibble::tibble(
        species = sp, ontology = "ALL",
        n_terms = length(cumulative), n_new_terms = length(new_terms),
        step = i
      ))
    }

    cumulative_df <- dplyr::bind_rows(rarefaction_all$BP, rarefaction_all$MF, rarefaction_all$CC, curve_all)

    if (nrow(cumulative_df) > 0) {
      openxlsx::write.xlsx(
        cumulative_df,
        file.path(comparison_dir, paste0("cumulative_by_origin_", tag, ".xlsx")),
        asTable = TRUE
      )

      ontology_cols <- c(BP = "#1b9e77", MF = "#d95f02", CC = "#7570b3", ALL = "#e7298a")

      mk_panel <- function(df, ont) {
        df_ont <- df %>% dplyr::filter(.data$ontology == ont)
        if (all(df_ont$species == "GOseq")) {
          df_ont <- df_ont %>% dplyr::mutate(species = factor(.data$species, levels = "GOseq"))
        } else {
          species_levels <- df_ont %>%
            dplyr::filter(.data$species != "GOseq") %>%
            dplyr::arrange(.data$step, .data$species) %>%
            dplyr::pull(.data$species) %>% unique()
          species_levels <- c("GOseq", species_levels)
          df_ont <- df_ont %>%
            dplyr::mutate(species = factor(.data$species, levels = species_levels)) %>%
            dplyr::arrange(.data$step)
        }
        ggplot2::ggplot(df_ont, ggplot2::aes(x = .data$species, y = .data$n_terms)) +
          ggplot2::geom_line(group = 1, linewidth = 1.2, color = ontology_cols[[ont]]) +
          ggplot2::geom_point(size = 3, color = ontology_cols[[ont]]) +
          ggplot2::labs(
            title = paste("Cumulative GO Term Curve -", ont),
            x = "Species Added (addition order)", y = "Cumulative Unique GO Terms"
          ) +
          ggplot2::theme_minimal(base_size = 13) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      }

      plots <- list(
        BP  = mk_panel(cumulative_df, "BP"),
        MF  = mk_panel(cumulative_df, "MF"),
        CC  = mk_panel(cumulative_df, "CC"),
        ALL = mk_panel(cumulative_df, "ALL")
      )

      combo <- plots$BP + plots$MF + plots$CC + plots$ALL +
        patchwork::plot_layout(ncol = 2) +
        patchwork::plot_annotation(
          title = if (bg_mode) "Cumulative GO Term Curves (True Consensus)"
          else                 "Cumulative GO Term Curves (Exploratory)"
        )

      ggplot2::ggsave(file.path(comparison_dir, paste0("cumulative_curve_combined_", tag, ".pdf")),
                      combo, width = 14, height = 10)
    } else {
      message("No cumulative rarefaction data for mode '", tag, "'. Skipping cumulative plot and export.")
    }

    # Permutation rarefaction (1000)
    set.seed(1L)
    n_perm <- 1000L
    rarefied_all <- list()
    for (ont in c("BP","MF","CC","ALL")) {
      base_terms <- if (identical(ont, "ALL")) base_combined else {
        consensus_pa %>%
          dplyr::filter(.data$ontology == ont, .data$pa_in_goseq == 1) %>%
          dplyr::pull(.data$term_id) %>% unique()
      }
      term_sets <- if (identical(ont, "ALL")) term_sets_all else {
        entries$species_mode %>%
          stats::setNames(stringr::str_remove_all(entries$species_mode, "^pa_in_|_nobg$")) %>%
          purrr::map(function(colname) {
            consensus_pa %>%
              dplyr::filter(if (ont != "ALL") .data$ontology == ont else TRUE,
                            .data[[colname]] == 1, .data$pa_in_goseq == 0) %>%
              dplyr::pull(.data$term_id) %>% unique()
          })
      }
      species_list <- names(term_sets)
      if (!length(species_list)) next

      rarem <- replicate(n_perm, {
        sp_order <- sample(species_list)
        cumulative <- base_terms
        cumsum <- numeric(length(sp_order))
        for (i in seq_along(sp_order)) {
          new_terms <- setdiff(term_sets[[sp_order[i]]], cumulative)
          cumulative <- union(cumulative, new_terms)
          cumsum[i] <- length(cumulative)
        }
        cumsum
      }, simplify = "matrix")

      stats_df <- tibble::tibble(
        species_index = seq_along(species_list) + 1L, # +1 for GOseq baseline
        mean_terms = rowMeans(rarem),
        ci_lower   = apply(rarem, 1, stats::quantile, probs = 0.025),
        ci_upper   = apply(rarem, 1, stats::quantile, probs = 0.975),
        ontology   = ont
      )
      rarefied_all[[ont]] <- stats_df
    }

    if (length(rarefied_all)) {
      rarefied_df <- dplyr::bind_rows(rarefied_all)
      openxlsx::write.xlsx(rarefied_df, file.path(comparison_dir, paste0("rarefaction_by_origin_", tag, ".xlsx")), asTable = TRUE)

      if (nrow(rarefied_df) > 0) {
        p2 <- ggplot2::ggplot(rarefied_df, ggplot2::aes(x = .data$species_index, y = .data$mean_terms, color = .data$ontology, group = .data$ontology)) +
          ggplot2::geom_line(linewidth = 1.2) +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper, fill = .data$ontology), alpha = 0.2, color = NA) +
          ggplot2::facet_wrap(~ontology, scales = "free_y") +
          ggplot2::labs(
            title = if (bg_mode) "Rarefaction Curve (1000 Perms): GOseq + g:Profiler (True Consensus)"
            else                 "Rarefaction Curve (1000 Perms): GOseq + All Significant Terms (Exploratory)",
            x = "Species Randomly Added (Index)", y = "Mean Cumulative Unique GO Terms"
          ) +
          ggplot2::theme_minimal(base_size = 13)
        ggplot2::ggsave(file.path(comparison_dir, paste0("rarefaction_curve_by_origin_", tag, ".pdf")), p2, width = 10, height = 6)
      }
    } else {
      message("No species contributed to permutation rarefaction for mode '", tag, "'. Skipping permutation outputs.")
    }

    invisible(NULL)
  }

  build_mode(TRUE)   # strict
  build_mode(FALSE)  # exploratory

  openxlsx::write.xlsx(consensus_pa, file = file.path(comparison_dir, "presence_absence_table.xlsx"), asTable = TRUE)
  invisible(NULL)
}

#' Summarize graph stats per ontology by reading network outputs
#' Tries (in order):
#'   1) root CSV: <network_dir>/summary_<label>.csv
#'   2) legacy CSV: <network_dir>/<label>/network_summary_gene_count.csv (from the Rmd)
#'   3) compute on the fly from GraphML: <network_dir>/<label>/network_<ONT>_gene_count.graphml
#' If computed, it also writes the root summary_<label>.csv to stabilize future runs.
#' @keywords internal
summarize_network_complexity <- function(label, size_by = "gene_count", network_dir = ".") {
  # prefer provided network_dir; tolerate both canonical ('networks') and legacy ('Network_analysis')
  cand_dirs <- c(
    network_dir,
    sub("Network_analysis$", "network_analysis", network_dir),
    sub("network_analysis$", "Network_analysis", network_dir),
    sub("networks$", "Network_analysis", network_dir),
    sub("Network_analysis$", "networks", network_dir)
  )
  base <- cand_dirs[which(dir.exists(cand_dirs))[1]]
  if (is.na(base)) base <- cand_dirs[1]

  root_csv  <- file.path(base, paste0("summary_", label, ".csv"))
  rmd_csv   <- file.path(base, label, paste0("network_summary_", size_by, ".csv"))
  graph_dir <- file.path(base, label)

  # 1) root summary_<label>.csv
  if (file.exists(root_csv)) {
    df <- readr::read_csv(root_csv, show_col_types = FALSE)
    return(
      df %>%
        dplyr::transmute(
          ontology    = .data$ontology,
          total_terms = .data$total_terms,
          total_edges = .data$total_edges,
          avg_degree  = .data$avg_degree
        ) %>%
        dplyr::rename(
          !!paste0("total_terms_", label) := .data$total_terms,
          !!paste0("total_edges_", label) := .data$total_edges,
          !!paste0("avg_degree_",  label) := .data$avg_degree
        )
    )
  }

  # 2) legacy Rmd summary
  if (file.exists(rmd_csv)) {
    df <- readr::read_csv(rmd_csv, show_col_types = FALSE)
    out <- df %>%
      dplyr::transmute(
        ontology    = .data$ontology,
        total_terms = .data$total_terms,
        total_edges = .data$total_edges,
        avg_degree  = .data$avg_degree
      )
    try(readr::write_csv(out, root_csv), silent = TRUE)
    return(
      out %>%
        dplyr::rename(
          !!paste0("total_terms_", label) := .data$total_terms,
          !!paste0("total_edges_", label) := .data$total_edges,
          !!paste0("avg_degree_",  label) := .data$avg_degree
        )
    )
  }

  # 3) compute from GraphML if present
  onts <- c("BP","MF","CC")
  rows <- lapply(onts, function(ont) {
    gml <- file.path(graph_dir, paste0("network_", ont, "_", size_by, ".graphml"))
    if (!file.exists(gml)) {
      return(tibble::tibble(ontology = ont, total_terms = NA_integer_, total_edges = NA_integer_, avg_degree = NA_real_))
    }
    g   <- igraph::read_graph(gml, format = "graphml")
    deg <- igraph::degree(g)
    tibble::tibble(
      ontology    = ont,
      total_terms = igraph::vcount(g),
      total_edges = igraph::ecount(g),
      avg_degree  = if (length(deg)) mean(as.numeric(deg), na.rm = TRUE) else NA_real_
    )
  })
  out <- dplyr::bind_rows(rows)

  if (any(!is.na(out$total_terms))) {
    try(readr::write_csv(out, root_csv), silent = TRUE)
  } else {
    message("ℹ No GraphML found for label '", label, "'. Returning NA skeleton.")
  }

  out %>%
    dplyr::rename(
      !!paste0("total_terms_", label) := .data$total_terms,
      !!paste0("total_edges_", label) := .data$total_edges,
      !!paste0("avg_degree_",  label) := .data$avg_degree
    )
}

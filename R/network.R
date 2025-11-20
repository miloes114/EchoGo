#' @name run_all_networks
#' @title Run GO Term Network Construction for EchoGO Results (bulletproof)
#' @description
#' Builds two network modes from the consensus table:
#' \itemize{
#'   \item \strong{with_bg (True Consensus)} = GOseq (with BG) + g:Profiler (with BG) + Consensus (with BG)
#'   \item \strong{with_bg_and_nobg (Exploratory)} = all of the above \emph{plus} the no-background sources
#' }
#' Then constructs per-ontology (BP, MF, CC) GO-term overlap networks.
#'
#' @param consensus_df Data frame from \code{build_consensus_table()} (must include
#'   \code{term_id}, \code{term_name}, \code{ontology}, \code{all_genes}, \code{origin},
#'   \code{significant_in_any}; optionally \code{consensus_score}, \code{consensus_score_all}, \code{in_goseq}, and GOseq p-values).
#' @param min_shared_genes Minimum shared genes to draw an edge when \code{edge_metric="shared"}. Default \code{2}.
#' @param size_by Node size: \code{"gene_count"} (default) or \code{"consensus_score"}.
#' @param outdir Base output directory (writes under \code{outdir/networks}).
#' @param sep_regex Regex used to split \code{all_genes}. Default \code{"[,;]"}.
#' @param min_gene_count Drop terms with < this many genes before building. Default \code{3}.
#' @param max_terms_per_ontology (kept for backward-compat; not used)
#' @param edge_metric \code{"shared"} (fast; uses \code{min_shared_genes}) or \code{"jaccard"} (uses \code{min_jaccard}). Default \code{"shared"}.
#' @param min_jaccard Minimum Jaccard when \code{edge_metric="jaccard"}. Default \code{0.1}.
#' @param use_sparse Use \pkg{Matrix} sparse path with \code{tcrossprod}. Default \code{TRUE}.
#' @param analysis_cap_per_ontology Cap used for graph construction + stats. Default \code{500}.
#' @param plotting_cap_per_ontology Cap used only for plotting. Default \code{200}.
#'
#' @return Writes GraphML, PDF/SVG, HTML and summary CSVs; invisibly returns \code{NULL}.
#' @export
run_all_networks <- function(consensus_df,
                             min_shared_genes = 2,
                             size_by = "gene_count",
                             outdir = getwd(),
                             sep_regex = "[,;]",
                             min_gene_count = 3,
                             max_terms_per_ontology = 200,  # retained but unused
                             edge_metric = "shared",
                             min_jaccard = 0.1,
                             use_sparse = TRUE,
                             analysis_cap_per_ontology = 500,
                             plotting_cap_per_ontology = 200) {

  # --- helpers (local scope) ---
  `%||%` <- function(a,b) if (!is.null(a)) a else b
  .truthy <- function(x) {
    if (is.logical(x)) return(!is.na(x) & x)
    if (is.numeric(x)) return(!is.na(x) & x != 0)
    if (is.character(x)) return(tolower(trimws(x)) %in% c("true","1","yes","y"))
    rep(FALSE, length(x))
  }
  .has_genes <- function(x) {
    z <- trimws(as.character(x))
    nzchar(z) & !is.na(z)
  }
  .norm_ont <- function(x) {
    x <- toupper(substr(as.character(x), 1, 2))
    x[x %in% c("BP","MF","CC")] <- x[x %in% c("BP","MF","CC")]
    x
  }
  .mk <- function(...) { p <- file.path(...); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
  .mirror_tree <- function(src, dest) {
    if (!dir.exists(src)) return(invisible(FALSE))
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    files <- list.files(src, recursive = TRUE, full.names = TRUE, no.. = TRUE)
    ok <- TRUE
    for (f in files) {
      if (dir.exists(f)) next
      rel <- sub(paste0("^", gsub("\\\\","\\\\\\\\", normalizePath(src, winslash="/", mustWork=FALSE))),
                 "", normalizePath(f, winslash="/", mustWork=FALSE))
      rel <- sub("^[/\\\\]", "", rel)
      tgt <- file.path(dest, rel)
      dir.create(dirname(tgt), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(f, tgt, overwrite = TRUE)) ok <- FALSE
    }
    cat("This legacy folder mirrors: ", basename(src), "/\n", file = file.path(dest, "__moved_to.txt"))
    invisible(ok)
  }

  legacy_on <- isTRUE(getOption("EchoGO.legacy_aliases", FALSE))

  # Ensure scores & non-empty
  if (!"consensus_score" %in% names(consensus_df)) {
    consensus_df <- .echogo_add_consensus_scores(consensus_df)
  }
  if (!nrow(consensus_df)) { message("No consensus rows; skipping networks."); return(invisible(NULL)) }

  canon_root  <- .mk(outdir, "networks")
  legacy_root <- file.path(outdir, "Network_analysis")

  # ---------- TRUE CONSENSUS (with background only) ----------
  message("🚀 Generating GO term networks (True Consensus — all with background)…")
  network_df_with_bg <- consensus_df %>%
    dplyr::mutate(
      ontology = .norm_ont(.data$ontology),
      sig_any  = .truthy(.data$significant_in_any),
      in_gsq   = if ("in_goseq" %in% names(.)) .truthy(.data$in_goseq) else FALSE,
      origin   = as.character(.data$origin %||% "")
    ) %>%
    dplyr::filter(
      .data$ontology %in% c("BP","MF","CC"),
      sig_any,
      .has_genes(.data$all_genes),
      in_gsq |
        (grepl("g:?profiler", origin, ignore.case = TRUE) & grepl("with\\s*bg", origin, ignore.case = TRUE)) |
        (grepl("consensus",  origin, ignore.case = TRUE) & grepl("with\\s*bg", origin, ignore.case = TRUE))
    )

  build_go_networks(
    df = network_df_with_bg,
    label = "with_bg",
    size_by = size_by,
    min_shared_genes = min_shared_genes,
    sep_regex = sep_regex,
    min_gene_count = min_gene_count,
    analysis_cap_per_ontology = analysis_cap_per_ontology,
    plotting_cap_per_ontology = plotting_cap_per_ontology,
    edge_metric = edge_metric,
    min_jaccard = min_jaccard,
    use_sparse = use_sparse,
    outdir = outdir
  )

  # ---------- EXPLORATORY (with + no background) ----------
  message("🚀 Generating GO term networks (Exploratory — with BG + no BG)…")
  network_df_with_bg_and_nobg <- consensus_df %>%
    dplyr::mutate(
      ontology = .norm_ont(.data$ontology),
      sig_any  = .truthy(.data$significant_in_any),
      in_gsq   = if ("in_goseq" %in% names(.)) .truthy(.data$in_goseq) else FALSE,
      origin   = as.character(.data$origin %||% "")
    ) %>%
    dplyr::filter(
      .data$ontology %in% c("BP","MF","CC"),
      sig_any,
      .has_genes(.data$all_genes),
      in_gsq |
        grepl("g:?profiler", origin, ignore.case = TRUE) |
        grepl("consensus",  origin, ignore.case = TRUE)
    )

  build_go_networks(
    df = network_df_with_bg_and_nobg,
    label = "with_bg_and_nobg",
    size_by = size_by,
    min_shared_genes = min_shared_genes,
    sep_regex = sep_regex,
    min_gene_count = min_gene_count,
    analysis_cap_per_ontology = analysis_cap_per_ontology,
    plotting_cap_per_ontology = plotting_cap_per_ontology,
    edge_metric = edge_metric,
    min_jaccard = min_jaccard,
    use_sparse = use_sparse,
    outdir = outdir
  )

  # ----- Complexity comparison (read if present; don't create legacy) -----
  base_dir_canon  <- file.path(outdir, "networks")
  base_dir_legacy <- file.path(outdir, "Network_analysis")

  summary_bg_path_canon  <- file.path(base_dir_canon,  "summary_with_bg.csv")
  summary_exp_path_canon <- file.path(base_dir_canon,  "summary_with_bg_and_nobg.csv")

  if (!file.exists(summary_bg_path_canon) &&
      dir.exists(base_dir_legacy) &&
      file.exists(file.path(base_dir_legacy, "summary_with_bg.csv"))) {
    file.copy(file.path(base_dir_legacy, "summary_with_bg.csv"), summary_bg_path_canon, overwrite = TRUE)
  }
  if (!file.exists(summary_exp_path_canon) &&
      dir.exists(base_dir_legacy) &&
      file.exists(file.path(base_dir_legacy, "summary_with_bg_and_nobg.csv"))) {
    file.copy(file.path(base_dir_legacy, "summary_with_bg_and_nobg.csv"), summary_exp_path_canon, overwrite = TRUE)
  }

  comparison_outpath <- file.path(base_dir_canon, "network_complexity_comparison.csv")
  if (file.exists(summary_bg_path_canon) && file.exists(summary_exp_path_canon)) {
    comparison <- tryCatch({
      summary_bg  <- readr::read_csv(summary_bg_path_canon,  show_col_types = FALSE)
      summary_exp <- readr::read_csv(summary_exp_path_canon, show_col_types = FALSE)
      dplyr::full_join(summary_bg, summary_exp, by = "ontology", suffix = c("_bg", "_exploratory")) %>%
        dplyr::mutate(
          delta_terms      = total_terms_exploratory - total_terms_bg,
          delta_edges      = total_edges_exploratory - total_edges_bg,
          delta_avg_degree = avg_degree_exploratory - avg_degree_bg
        )
    }, error = function(e) { message("⚠️ Skipping network complexity comparison (read error): ", e$message); NULL })
    if (!is.null(comparison)) {
      readr::write_csv(comparison, comparison_outpath)
      print(comparison)
    }
  } else {
    message("ℹ Skipping network complexity comparison: summary files not found for both modes.")
  }

  # Mirror canonical tree → legacy ONLY if explicitly enabled
  if (legacy_on && dir.exists(canon_root)) .mirror_tree(canon_root, legacy_root)

  message("✅ Network generation complete.")
  invisible(NULL)
}


#' @name build_go_networks
#' @title Build Per-Ontology GO Term Networks (capped + scalable)
#' @description
#' Constructs GO-term overlap networks using a pre-capped set of terms per ontology,
#' ranked fairly: \emph{Consensus/g:Profiler} terms by the relevant consensus score, and
#' \emph{GOseq-only} terms by \code{−log10(p)} (auto-detected p column). Uses a
#' sparse \code{Matrix::tcrossprod} path by default for speed.
#'
#' @param df Subset of consensus table (columns: \code{term_id}, \code{term_name}, \code{ontology}, \code{all_genes},
#'   \code{origin}, \code{in_goseq}, \code{significant_in_any}, \code{consensus_score}, optionally \code{consensus_score_all} and a GOseq p-value column).
#' @param label Output subfolder under \code{outdir/networks}. (Expected: "with_bg" or "with_bg_and_nobg")
#' @param size_by \code{"gene_count"} (default) or \code{"consensus_score"} for node sizing.
#' @param min_shared_genes Edge threshold when \code{edge_metric="shared"}. Default \code{2}.
#' @param sep_regex Regex for splitting \code{all_genes}. Default \code{"[,;]"}.
#' @param min_gene_count Drop terms with < this many genes before ranking. Default \code{3}.
#' @param analysis_cap_per_ontology Hard cap per ontology before edges (for graph+stats). Default \code{500}.
#' @param plotting_cap_per_ontology Cap used only for plotting. Default \code{200}.
#' @param edge_metric \code{"shared"} or \code{"jaccard"}. Default \code{"shared"}.
#' @param min_jaccard Minimum Jaccard when \code{edge_metric="jaccard"}. Default \code{0.1}.
#' @param use_sparse Use sparse overlap path. Default \code{TRUE}.
#' @param outdir Base output directory.
#' @return Writes files; returns \code{NULL} (invisibly).
#' @export
build_go_networks <- function(df,
                              label,
                              size_by = "gene_count",
                              min_shared_genes = 2,
                              sep_regex = "[,;]",
                              min_gene_count = 3,
                              analysis_cap_per_ontology = 500,
                              plotting_cap_per_ontology = 200,
                              edge_metric = "shared",
                              min_jaccard = 0.1,
                              use_sparse = TRUE,
                              outdir = getwd()) {

  stopifnot(edge_metric %in% c("shared","jaccard"))

  # --- helpers (local scope) ---
  `%||%` <- function(a,b) if (!is.null(a)) a else b
  .mk <- function(...) { p <- file.path(...); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
  .count_nonempty <- function(x, sep_regex="[,;]") {
    toks <- trimws(unlist(strsplit(ifelse(is.na(x),"",x), sep_regex)))
    sum(nzchar(toks))
  }
  .find_p <- function(d) {
    # explicit candidates first
    cand <- c(
      "over_represented_pvalue","padj","p_adj","pvalue_adj","p_value_adj",
      "p_over","pval","p_value","p",
      # a few more common variants seen in the wild
      "p.over","p.adjust","p_adj_over","adj_p","adj.p"
    )
    hits <- intersect(cand, names(d))
    if (length(hits)) return(hits[1])

    # heuristic: any numeric column whose name looks like a p-value
    nm <- names(d)
    numeric_cols <- nm[vapply(d, is.numeric, logical(1))]
    regex_hits <- numeric_cols[grepl("^p$|p[._]?val|p[._]?value|over[._]?represented[._]?p", numeric_cols, ignore.case = TRUE)]
    if (length(regex_hits)) return(regex_hits[1])

    # nothing found -> return NULL (not NA!)
    NULL
  }
    .truthy <- function(x) {
    if (is.logical(x)) return(!is.na(x) & x)
    if (is.numeric(x)) return(!is.na(x) & x != 0)
    if (is.character(x)) return(tolower(trimws(x)) %in% c("true","1","yes","y"))
    rep(FALSE, length(x))
  }
  .norm_ont <- function(x) {
    x <- toupper(substr(as.character(x), 1, 2))
    x[x %in% c("BP","MF","CC")] <- x[x %in% c("BP","MF","CC")]
    x
  }
  .rank_attr <- function(d, label, sep_regex) {
    pcol <- .find_p(d)
    if (isTRUE(is.na(pcol)) || length(pcol) == 0) pcol <- NULL
        is_goseq <- (("in_goseq" %in% names(d)) & .truthy(d$in_goseq)) |
      grepl("goseq", d$origin %||% "", ignore.case = TRUE)
    cons_col <- if (identical(label, "with_bg_and_nobg") && "consensus_score_all" %in% names(d))
      "consensus_score_all" else "consensus_score"
    cons <- suppressWarnings(as.numeric(d[[cons_col]]))
    rank_val <- cons
    if (!is.null(pcol)) {
      pvec <- suppressWarnings(as.numeric(d[[pcol]]))
      use_p <- is_goseq & is.finite(pvec)
      rank_val[use_p] <- -log10(pmax(pvec[use_p], 1e-300))
    }
    bad <- !is.finite(rank_val)
    if (any(bad)) {
      gcnt <- if ("gene_count" %in% names(d)) d$gene_count else vapply(d$all_genes, function(x) .count_nonempty(x, sep_regex), integer(1))
      rank_val[bad] <- gcnt[bad]
    }
    rank_val
  }

  if (is.null(df) || !nrow(df)) { message("No rows for ", label, "; skipping."); return(invisible(NULL)) }
  if (!"consensus_score" %in% names(df)) df <- .echogo_add_consensus_scores(df)

  base_dir   <- .mk(outdir, "networks")
  output_dir <- .mk(base_dir, label)

  stats_rows <- list()
  set.seed(1L)

  for (ont in c("BP","MF","CC")) {
    message("🔵 [", label, "] ", ont, " — selecting terms")

    df_net <- df %>%
      dplyr::mutate(ontology = .norm_ont(.data$ontology)) %>%
      dplyr::filter(.data$ontology == ont, .truthy(.data$significant_in_any), !is.na(.data$term_id)) %>%
      dplyr::mutate(
        id    = .data$term_id,
        label = .data$term_name,
        gene_count = vapply(.data$all_genes, .count_nonempty, integer(1), sep_regex = sep_regex),
        consensus_score = suppressWarnings(as.numeric(.data$consensus_score))
      ) %>%
      dplyr::filter(.data$gene_count >= min_gene_count)

    if (nrow(df_net) < 2) {
      message("⏩ Skipping ", ont, ": too few terms after filter.")
      stats_rows[[length(stats_rows)+1]] <- tibble::tibble(ontology=ont,total_terms=0,total_edges=0,avg_degree=0)
      next
    }

    # ---- FAIR RANK (consensus-aware) + DETERMINISTIC TIEBREAKS ----
    df_net$rank_attr  <- .rank_attr(df_net, label, sep_regex)
    df_net$size_attr  <- ifelse(identical(size_by,"gene_count"), df_net$gene_count, df_net$consensus_score)

    if (is.finite(analysis_cap_per_ontology) && nrow(df_net) > analysis_cap_per_ontology) {
      avg_fold_sym <- if ("avg_fold_gprof_bg" %in% names(df_net)) rlang::sym("avg_fold_gprof_bg") else NULL
      df_net <- df_net %>%
        dplyr::arrange(
          dplyr::desc(.data$rank_attr),
          dplyr::desc(.data$gene_count),
          if (!is.null(avg_fold_sym)) dplyr::desc(!!avg_fold_sym) else dplyr::desc(.data$rank_attr),
          .data$term_id
        ) %>%
        dplyr::slice_head(n = analysis_cap_per_ontology)
      message("⚖️  ", ont, " (", label, "): analysis capped to ", nrow(df_net), " terms.")
    } else {
      df_net <- df_net %>%
        dplyr::arrange(
          dplyr::desc(.data$rank_attr),
          dplyr::desc(.data$gene_count),
          .data$term_id
        )
    }

    # ----- EDGE CONSTRUCTION -----
    message("🧮 ", ont, ": building overlaps (", if (use_sparse) "sparse" else "pairwise", "; metric=", edge_metric, ")")
    term_edges <- NULL

    tg <- df_net %>%
      dplyr::filter(!is.na(.data$all_genes)) %>%
      tidyr::separate_rows(all_genes, sep = sep_regex) %>%
      dplyr::mutate(all_genes = trimws(.data$all_genes)) %>%
      dplyr::filter(nzchar(.data$all_genes))

    ids <- unique(tg$id)
    if (length(ids) < 2) {
      message("⏩ Skipping ", ont, ": fewer than 2 terms with genes.")
      stats_rows[[length(stats_rows)+1]] <- tibble::tibble(ontology=ont,total_terms=dplyr::n_distinct(df_net$id),total_edges=0,avg_degree=0)
      next
    }

    if (use_sparse) {
      if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package needed for use_sparse=TRUE.")
      term_index <- match(tg$id, ids)
      gene_levels <- unique(tg$all_genes)
      gene_index <- match(tg$all_genes, gene_levels)

      M <- Matrix::sparseMatrix(i = term_index, j = gene_index, x = 1L,
                                dims = c(length(ids), length(gene_levels)),
                                dimnames = list(ids, gene_levels))

      S <- Matrix::tcrossprod(M)  # shared counts
      Matrix::diag(S) <- 0
      T <- Matrix::summary(S)     # i, j, x (non-zero)

      if (edge_metric == "shared") {
        if (nrow(T)) {
          T <- T[T$i < T$j & T$x >= min_shared_genes, , drop = FALSE]
          if (nrow(T)) {
            term_edges <- tibble::tibble(
              from   = rownames(S)[T$i],
              to     = colnames(S)[T$j],
              weight = as.numeric(T$x)
            )
          }
        }
      } else {
        if (nrow(T)) {
          deg <- Matrix::rowSums(M)
          T <- T[T$i < T$j & T$x >= 1, , drop = FALSE]
          if (nrow(T)) {
            numer <- as.numeric(T$x)
            denom <- deg[T$i] + deg[T$j] - numer
            jac   <- ifelse(denom > 0, numer/denom, 0)
            keep  <- jac >= min_jaccard
            if (any(keep)) {
              T <- T[keep, , drop = FALSE]
              term_edges <- tibble::tibble(
                from   = rownames(S)[T$i],
                to     = colnames(S)[T$j],
                weight = jac[keep]
              )
            }
          }
        }
      }

    } else {
      term_pairs <- utils::combn(ids, 2, simplify = FALSE)
      term_edges <- purrr::map_dfr(term_pairs, function(pair) {
        g1 <- tg$all_genes[tg$id == pair[1]]
        g2 <- tg$all_genes[tg$id == pair[2]]
        shared <- length(intersect(g1, g2))
        if (edge_metric == "shared") {
          if (shared >= min_shared_genes) tibble::tibble(from=pair[1], to=pair[2], weight=shared) else NULL
        } else {
          denom <- length(unique(c(g1, g2)))
          jac   <- if (denom > 0) shared/denom else 0
          if (jac >= min_jaccard) tibble::tibble(from=pair[1], to=pair[2], weight=jac) else NULL
        }
      })
    }

    if (is.null(term_edges) || !nrow(term_edges)) {
      message("⏩ Skipping ", ont, ": no edges with current thresholds.")
      stats_rows[[length(stats_rows)+1]] <- tibble::tibble(ontology=ont,total_terms=dplyr::n_distinct(df_net$id),total_edges=0,avg_degree=0)
      next
    }

    # ----- GRAPH + OUTPUTS -----
    nodes <- df_net %>%
      dplyr::transmute(
        id, label, genes_flat = ifelse(is.na(.data$all_genes),"",.data$all_genes),
        gene_count, consensus_score,
        size_attr = ifelse(identical(size_by,"gene_count"), gene_count, consensus_score)
      )

    g <- igraph::graph_from_data_frame(term_edges, vertices = nodes, directed = FALSE)
    igraph::V(g)$degree <- igraph::degree(g)
    comm <- igraph::cluster_louvain(g)$membership
    igraph::V(g)$community <- comm

    # community labels: top by degree within each community
    top_terms <- igraph::as_data_frame(g, "vertices") %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(community = comm, degree = igraph::V(g)$degree) %>%
      dplyr::group_by(community) %>%
      dplyr::arrange(dplyr::desc(degree)) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup() %>%
      dplyr::select(community, label) %>%
      dplyr::rename(top_term = label)
    top_map <- stats::setNames(top_terms$top_term, top_terms$community)
    igraph::V(g)$community_label <- paste0("Cluster: ", top_map[as.character(comm)])

    # write graphml
    graph_path <- file.path(output_dir, paste0("network_", ont, "_", size_by, ".graphml"))
    igraph::write_graph(g, graph_path, format = "graphml")

    # ---- Rmd-style layout + labels + ggsave (plots only) ----
    layout_tbl <- ggraph::create_layout(g, layout = "fr") %>%
      dplyr::mutate(
        community = igraph::V(g)$community,
        degree    = igraph::degree(g)
      )

    # join top term per community for legend labels (same as Rmd)
    top_terms_plot <- igraph::as_data_frame(g, what = "vertices") %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        community = igraph::V(g)$community,
        degree    = igraph::degree(g)
      ) %>%
      dplyr::group_by(community) %>%
      dplyr::arrange(dplyr::desc(degree)) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup() %>%
      dplyr::select(community, label) %>%
      dplyr::rename(top_term = label)

    layout_tbl <- layout_tbl %>%
      dplyr::left_join(top_terms_plot, by = "community") %>%
      dplyr::mutate(community_label = paste0("Cluster: ", top_term)) %>%
      dplyr::arrange(dplyr::desc(size_attr)) %>%
      dplyr::slice_head(n = plotting_cap_per_ontology)

    top_label_nodes <- layout_tbl %>%
      dplyr::group_by(community) %>%
      dplyr::arrange(dplyr::desc(size_attr)) %>%
      dplyr::slice_head(n = 3) %>%
      dplyr::ungroup() %>%
      dplyr::filter(
        !is.na(x), !is.na(y),
        is.finite(x), is.finite(y),
        !is.na(label), label != ""
      )

    p_net <- ggraph::ggraph(layout_tbl) +
      ggraph::geom_edge_link(alpha = 0.2) +
      ggraph::geom_node_point(
        ggplot2::aes(size = degree, color = community_label),
        alpha = 0.8
      ) +
      ggrepel::geom_label_repel(
        data = top_label_nodes,
        ggplot2::aes(x = x, y = y, label = label),
        box.padding = 0.4, max.overlaps = Inf, size = 4, fill = "white", alpha = 0.7
      ) +
      ggplot2::scale_size(range = c(3, 10)) +
      ggplot2::guides(color = ggplot2::guide_legend(title = "Top GO Term per Cluster")) +
      ggplot2::theme_void() +
      ggplot2::ggtitle(paste("GO", ont, "Network (", label, ", >", min_shared_genes, " shared genes)"))

    pdf_path <- file.path(output_dir, paste0("network_", ont, "_", size_by, "_filtered.pdf"))
    svg_path <- file.path(output_dir, paste0("network_", ont, "_", size_by, "_filtered.svg"))

    try(ggplot2::ggsave(filename = pdf_path, plot = p_net, width = 12, height = 10, limitsize = FALSE), silent = TRUE)
    if (requireNamespace("svglite", quietly = TRUE)) {
      try(ggplot2::ggsave(filename = svg_path, plot = p_net, width = 12, height = 10, device = svglite::svglite, limitsize = FALSE), silent = TRUE)
    }

    # Interactive visNetwork (labels for top 3 per community by size_attr)
    coords_vis <- igraph::layout_nicely(g)
    nodes_tbl  <- igraph::as_data_frame(g, what="vertices") %>% tibble::as_tibble()
    top_ids    <- nodes_tbl %>% dplyr::mutate(size_attr = ifelse(identical(size_by,"gene_count"), gene_count, consensus_score)) %>%
      dplyr::group_by(community) %>% dplyr::arrange(dplyr::desc(size_attr)) %>% dplyr::slice_head(n=3) %>% dplyr::pull(name)

    nodes_vis <- nodes_tbl %>%
      dplyr::mutate(
        id = name, group = as.character(community),
        hiddenLabel = ifelse(is.na(label)|label=="", id, label),
        label = ifelse(name %in% top_ids, hiddenLabel, ""),
        title = paste0("<b>Term:</b> ", hiddenLabel,
                       "<br><b>Genes:</b> ", ifelse(is.na(genes_flat),"",genes_flat),
                       "<br><b>Degree:</b> ", degree,
                       "<br><b>Gene count:</b> ", gene_count),
        x = coords_vis[,1]*100, y = coords_vis[,2]*100, value = degree
      )
    edges_vis <- igraph::as_data_frame(g, what="edges") %>% tibble::as_tibble()

    vis_obj <- visNetwork::visNetwork(nodes_vis, edges_vis) %>%
      visNetwork::visOptions(highlightNearest=TRUE, nodesIdSelection=TRUE) %>%
      visNetwork::visInteraction(dragNodes=TRUE, dragView=TRUE, zoomView=TRUE) %>%
      visNetwork::visIgraphLayout(layout="layout_nicely", physics=FALSE) %>%
      visNetwork::visEvents(click = "function(params) {
        var node = this.body.data.nodes.get(params.nodes[0]);
        if (node && node.hiddenLabel !== undefined) {
          var newLabel = (node.label === '') ? node.hiddenLabel : '';
          this.body.data.nodes.update({id: node.id, label: newLabel});
        }
      }")
    html_path <- file.path(output_dir, paste0("network_", ont, "_", size_by, "_filtered.html"))
    try(htmlwidgets::saveWidget(vis_obj, file = normalizePath(html_path, mustWork = FALSE), selfcontained = TRUE), silent = TRUE)

    # stats
    stats_rows[[length(stats_rows)+1]] <- tibble::tibble(
      ontology    = ont,
      total_terms = igraph::vcount(g),
      total_edges = igraph::gsize(g),
      avg_degree  = round(mean(igraph::degree(g)), 2)
    )
  }

  # summary
  summary_path <- file.path(base_dir, paste0("summary_", label, ".csv"))
  readr::write_csv(dplyr::bind_rows(stats_rows), summary_path)
  message("📝 Wrote summary: ", summary_path)

  invisible(NULL)
}

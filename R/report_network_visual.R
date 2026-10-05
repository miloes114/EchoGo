# Scoreless report network visualization -------------------------------------

.echogo_split_genes <- function(x, sep_regex = "[,;]") {
  vals <- unique(trimws(unlist(strsplit(ifelse(is.na(x), "", x), sep_regex))))
  vals[nzchar(vals)]
}

.echogo_report_network_data <- function(evidence, profiles, ontology,
                                        min_shared_genes = 2L,
                                        min_gene_count = 2L,
                                        max_terms = 120L) {
  empty_diag <- tibble::tibble(
    ontology = as.character(ontology),
    product = paste(profiles, collapse = "+"),
    input_terms = 0L,
    terms_with_gene_annotations = 0L,
    terms_after_gene_filter = 0L,
    candidate_pairs = 0L,
    pairs_sharing_at_least_one_gene = 0L,
    pairs_meeting_threshold = 0L,
    final_edges = 0L,
    min_shared_genes = as.integer(min_shared_genes),
    min_gene_count = as.integer(min_gene_count),
    status = "SKIPPED_NO_INPUT_TERMS"
  )
  if (is.null(evidence) || !is.data.frame(evidence) || !nrow(evidence)) {
    return(list(nodes = tibble::tibble(), edges = tibble::tibble(), diagnostics = empty_diag))
  }

  selected <- evidence |>
    dplyr::filter(
      .data$evidence_profile %in% profiles,
      .data$ontology %in% c(ontology, paste0("GO:", ontology))
    )

  input_terms <- nrow(selected)
  with_genes <- selected |>
    dplyr::filter(!is.na(.data$contributing_genes), nzchar(trimws(.data$contributing_genes)))

  nodes <- with_genes |>
    dplyr::mutate(
      gene_count = vapply(
        .data$contributing_genes,
        function(x) length(.echogo_split_genes(x)),
        integer(1)
      )
    ) |>
    dplyr::filter(.data$gene_count >= min_gene_count) |>
    dplyr::arrange(.data$display_order, dplyr::desc(.data$gene_count), .data$term_id) |>
    dplyr::slice_head(n = max_terms) |>
    dplyr::transmute(
      term_id = as.character(.data$term_id),
      term_name = .data$term_name,
      evidence_profile = .data$evidence_profile,
      alternative_context_support_n = .data$alternative_context_support_n,
      alternative_queried_context_n = .data$alternative_queried_context_n,
      alternative_context_support_fraction = .data$alternative_context_support_fraction,
      contributing_genes = .data$contributing_genes,
      gene_count = .data$gene_count
    )

  edges <- tibble::tibble(from = character(), to = character(), shared_gene_n = integer())
  candidate_pairs <- 0L
  pairs_ge_1 <- 0L
  pairs_ge_threshold <- 0L
  if (nrow(nodes) >= 2L) {
    pairs <- utils::combn(seq_len(nrow(nodes)), 2L, simplify = FALSE)
    candidate_pairs <- length(pairs)
    edge_rows <- lapply(pairs, function(pair) {
      left <- .echogo_split_genes(nodes$contributing_genes[[pair[[1]]]])
      right <- .echogo_split_genes(nodes$contributing_genes[[pair[[2]]]])
      shared <- length(intersect(left, right))
      if (shared >= 1L) pairs_ge_1 <<- pairs_ge_1 + 1L
      if (shared < min_shared_genes) return(NULL)
      pairs_ge_threshold <<- pairs_ge_threshold + 1L
      tibble::tibble(
        from = nodes$term_id[[pair[[1]]]],
        to = nodes$term_id[[pair[[2]]]],
        shared_gene_n = shared
      )
    })
    edges <- dplyr::bind_rows(edge_rows)
  }

  status <- dplyr::case_when(
    input_terms == 0L ~ "SKIPPED_NO_INPUT_TERMS",
    nrow(with_genes) == 0L ~ "SKIPPED_NO_GENE_ANNOTATIONS",
    nrow(nodes) < 2L ~ "SKIPPED_TOO_FEW_TERMS_AFTER_GENE_FILTER",
    nrow(edges) == 0L & pairs_ge_1 == 0L ~ "SKIPPED_NO_SHARED_GENES",
    nrow(edges) == 0L ~ "SKIPPED_BELOW_SHARED_GENE_THRESHOLD",
    TRUE ~ "GENERATED"
  )

  diagnostics <- tibble::tibble(
    ontology = as.character(ontology),
    product = paste(profiles, collapse = "+"),
    input_terms = as.integer(input_terms),
    terms_with_gene_annotations = as.integer(nrow(with_genes)),
    terms_after_gene_filter = as.integer(nrow(nodes)),
    candidate_pairs = as.integer(candidate_pairs),
    pairs_sharing_at_least_one_gene = as.integer(pairs_ge_1),
    pairs_meeting_threshold = as.integer(pairs_ge_threshold),
    final_edges = as.integer(nrow(edges)),
    min_shared_genes = as.integer(min_shared_genes),
    min_gene_count = as.integer(min_gene_count),
    status = status
  )

  list(nodes = nodes, edges = edges, diagnostics = diagnostics)
}

.echogo_network_skip_message <- function(diagnostics) {
  if (is.null(diagnostics) || !is.data.frame(diagnostics) || !nrow(diagnostics)) {
    return("No network diagnostic information was available.")
  }
  d <- diagnostics[1, , drop = FALSE]
  switch(
    d$status[[1]],
    SKIPPED_NO_INPUT_TERMS = "No qualifying GO terms were available for this network product.",
    SKIPPED_NO_GENE_ANNOTATIONS = paste0(
      d$input_terms[[1]],
      " qualifying terms were available, but none contained contributing-gene labels that could be used to construct a gene-overlap network."
    ),
    SKIPPED_TOO_FEW_TERMS_AFTER_GENE_FILTER = paste0(
      d$terms_with_gene_annotations[[1]],
      " terms contained contributing-gene labels, but fewer than two terms remained after requiring at least ",
      d$min_gene_count[[1]], " contributing genes per term."
    ),
    SKIPPED_NO_SHARED_GENES = paste0(
      d$terms_after_gene_filter[[1]],
      " terms passed the gene filter, but no pair shared a contributing gene."
    ),
    SKIPPED_BELOW_SHARED_GENE_THRESHOLD = paste0(
      d$terms_after_gene_filter[[1]], " terms passed the gene filter and ",
      d$pairs_sharing_at_least_one_gene[[1]],
      " term pairs shared at least one gene, but no pair shared the configured minimum of ",
      d$min_shared_genes[[1]], " genes."
    ),
    GENERATED = paste0(
      d$terms_after_gene_filter[[1]], " terms produced ", d$final_edges[[1]],
      " shared-gene edges."
    ),
    paste0("Network status: ", d$status[[1]], ".")
  )
}

#' Plot a gene-overlap network from scoreless EchoGO evidence
#' @keywords internal
plot_echogo_gene_overlap_network <- function(evidence, profiles, ontology,
                                             min_shared_genes = 2L,
                                             min_gene_count = 2L,
                                             max_terms = 120L,
                                             label_per_community = 3L) {
  dat <- .echogo_report_network_data(
    evidence = evidence,
    profiles = profiles,
    ontology = ontology,
    min_shared_genes = min_shared_genes,
    min_gene_count = min_gene_count,
    max_terms = max_terms
  )
  nodes <- dat$nodes
  edges <- dat$edges
  if (!nrow(nodes) || !nrow(edges)) {
    return(structure(NULL, echogo_network_diagnostics = dat$diagnostics))
  }

  g <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  igraph::V(g)$degree <- igraph::degree(g)
  membership <- if (igraph::ecount(g) > 0L && igraph::vcount(g) > 1L) {
    igraph::cluster_louvain(g)$membership
  } else {
    seq_len(igraph::vcount(g))
  }
  igraph::V(g)$community <- membership

  # igraph normalizes the vertex-key column to `name`. Restore the canonical
  # EchoGO term identifier explicitly before deterministic ordering. The prior
  # report renderer expected `term_id` here and failed silently for every
  # otherwise-valid graph.
  vertex_tbl <- igraph::as_data_frame(g, what = "vertices") |>
    tibble::as_tibble() |>
    dplyr::mutate(
      term_id = as.character(.data$name),
      degree = igraph::V(g)$degree,
      community = membership
    )

  top_terms <- vertex_tbl |>
    dplyr::group_by(.data$community) |>
    dplyr::arrange(dplyr::desc(.data$degree), dplyr::desc(.data$gene_count), .data$term_id) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup() |>
    dplyr::select(.data$community, top_term = .data$term_name)

  layout_tbl <- ggraph::create_layout(g, layout = "fr") |>
    dplyr::mutate(
      term_id = as.character(.data$name),
      community = membership,
      degree = igraph::V(g)$degree
    ) |>
    dplyr::left_join(top_terms, by = "community") |>
    dplyr::mutate(community_label = paste0("Module: ", .data$top_term))

  label_nodes <- layout_tbl |>
    dplyr::group_by(.data$community) |>
    dplyr::arrange(dplyr::desc(.data$degree), dplyr::desc(.data$gene_count), .data$term_id) |>
    dplyr::slice_head(n = label_per_community) |>
    dplyr::ungroup() |>
    dplyr::filter(
      is.finite(.data$x), is.finite(.data$y),
      !is.na(.data$term_name), nzchar(.data$term_name)
    )

  edge_max <- max(edges$shared_gene_n, na.rm = TRUE)
  edge_scale_max <- if (is.finite(edge_max) && edge_max > 0) edge_max else 1
  profile_shapes <- c(TARGET_ONLY = 21, TARGET_PLUS_CONTEXT = 22, ALTERNATIVE_CONTEXT = 21)

  p <- ggraph::ggraph(layout_tbl) +
    ggraph::geom_edge_link(
      ggplot2::aes(width = .data$shared_gene_n),
      alpha = .18, colour = "#80919c", show.legend = TRUE
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(
        size = .data$gene_count,
        fill = .data$community_label,
        shape = .data$evidence_profile
      ),
      colour = "#314854", stroke = .8, alpha = .9
    ) +
    ggplot2::scale_shape_manual(values = profile_shapes, drop = FALSE) +
    ggplot2::scale_size(range = c(3.5, 10), name = "Contributing genes") +
    ggraph::scale_edge_width(range = c(.3, 1.5), limits = c(1, edge_scale_max), name = "Shared genes") +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_void() +
    ggplot2::labs(
      title = paste0("GO ", ontology, " gene-overlap network"),
      subtitle = paste0("Edges join terms sharing at least ", min_shared_genes, " contributing genes.")
    )

  if (nrow(label_nodes)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_label_repel(
        data = label_nodes,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$term_name),
        box.padding = .4, max.overlaps = Inf, size = 3.4,
        fill = "white", alpha = .84
      )
    } else {
      p <- p + ggplot2::geom_label(
        data = label_nodes,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$term_name),
        size = 3.1, fill = "white", alpha = .84
      )
    }
  }

  attr(p, "echogo_network_nodes") <- nodes
  attr(p, "echogo_network_edges") <- edges
  attr(p, "echogo_network_diagnostics") <- dat$diagnostics
  p
}

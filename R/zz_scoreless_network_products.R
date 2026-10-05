# Canonical scoreless network products ---------------------------------------
# Loaded after network.R to replace the earlier machine-only primary network
# writer with two evidence-safe biological products and restored visuals.

.echogo_run_scoreless_primary_networks <- function(evidence, outdir,
                                                   min_shared_genes = 2,
                                                   min_gene_count = 3,
                                                   sep_regex = "[,;]") {
  products <- list(
    target_supported = c("TARGET_ONLY", "TARGET_PLUS_CONTEXT"),
    alternative_context_hypothesis = "ALTERNATIVE_CONTEXT"
  )
  summary_rows <- list()

  for (product in names(products)) {
    root <- file.path(outdir, "networks", product)
    dir.create(root, recursive = TRUE, showWarnings = FALSE)

    for (ont in c("BP", "MF", "CC")) {
      dat <- .echogo_report_network_data(
        evidence = evidence,
        profiles = products[[product]],
        ontology = ont,
        min_shared_genes = min_shared_genes,
        min_gene_count = min_gene_count,
        max_terms = 500L
      )
      nodes <- dat$nodes
      edges <- dat$edges

      readr::write_csv(nodes, file.path(root, paste0("network_nodes_", ont, ".csv")))
      readr::write_csv(edges, file.path(root, paste0("network_edges_", ont, ".csv")))

      graph <- NULL
      if (nrow(nodes)) {
        graph <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
        igraph::write_graph(
          graph,
          file.path(root, paste0("network_", ont, ".graphml")),
          format = "graphml"
        )
      }

      p <- tryCatch(
        plot_echogo_gene_overlap_network(
          evidence = evidence,
          profiles = products[[product]],
          ontology = ont,
          min_shared_genes = min_shared_genes,
          min_gene_count = min_gene_count,
          max_terms = 120L
        ),
        error = function(e) NULL
      )
      if (!is.null(p)) {
        .echogo_save_plot_variants(
          p, root, paste0("network_", ont, "_filtered"),
          width = 13, height = 10, dpi = 320
        )
      }

      if (!is.null(graph) && igraph::vcount(graph) &&
          requireNamespace("visNetwork", quietly = TRUE) &&
          requireNamespace("htmlwidgets", quietly = TRUE)) {
        vertex <- igraph::as_data_frame(graph, what = "vertices") |>
          tibble::as_tibble()
        if (!"name" %in% names(vertex) && "term_id" %in% names(vertex)) vertex$name <- vertex$term_id
        vertex$degree <- igraph::degree(graph)
        coords <- igraph::layout_nicely(graph)
        degree_cut <- if (length(vertex$degree)) {
          as.numeric(stats::quantile(vertex$degree, .75, na.rm = TRUE, names = FALSE))
        } else 0
        vertex <- vertex |>
          dplyr::mutate(
            id = .data$name,
            label = dplyr::if_else(
              .data$degree >= degree_cut,
              dplyr::coalesce(.data$term_name, .data$name),
              ""
            ),
            title = paste0(
              "<b>GO term:</b> ", dplyr::coalesce(.data$term_name, .data$name),
              "<br><b>Evidence:</b> ", .data$evidence_profile,
              "<br><b>Contributing genes:</b> ", .data$gene_count,
              "<br><b>Alternative contexts:</b> ", .data$alternative_context_support_n,
              "/", .data$alternative_queried_context_n,
              "<br><b>Degree:</b> ", .data$degree
            ),
            value = pmax(1, .data$gene_count),
            x = coords[,1] * 100,
            y = coords[,2] * 100
          )
        edge_tbl <- igraph::as_data_frame(graph, what = "edges") |>
          tibble::as_tibble()
        if (!"width" %in% names(edge_tbl)) {
          edge_tbl$width <- if ("shared_gene_n" %in% names(edge_tbl)) edge_tbl$shared_gene_n else 1
        }
        vis <- visNetwork::visNetwork(vertex, edge_tbl) |>
          visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) |>
          visNetwork::visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE) |>
          visNetwork::visIgraphLayout(layout = "layout_nicely", physics = FALSE)
        html <- file.path(root, paste0("network_", ont, "_filtered.html"))
        ok <- tryCatch({
          htmlwidgets::saveWidget(vis, html, selfcontained = TRUE)
          TRUE
        }, error = function(e) FALSE)
        if (!ok) {
          try(htmlwidgets::saveWidget(vis, html, selfcontained = FALSE), silent = TRUE)
        }
      }

      summary_rows[[length(summary_rows) + 1L]] <- tibble::tibble(
        product = product,
        ontology = ont,
        nodes = nrow(nodes),
        edges = nrow(edges),
        min_shared_genes = min_shared_genes,
        min_gene_count = min_gene_count
      )
    }
  }

  summary <- dplyr::bind_rows(summary_rows)
  summary_path <- file.path(outdir, "networks", "scoreless_network_summary.csv")
  readr::write_csv(summary, summary_path)

  # Preserve the established machine-readable aggregate path for downstream
  # consumers and older validation contracts.  The biological report uses the
  # two partitioned products above; this compatibility export is only a view
  # of the same non-NO_PRIMARY_SUPPORT rows and does not create a third
  # interpretive product.
  compatibility_root <- file.path(outdir, "networks", "primary_custom_evidence")
  dir.create(compatibility_root, recursive = TRUE, showWarnings = FALSE)
  for (ont in c("BP", "MF", "CC")) {
    dat <- .echogo_report_network_data(
      evidence = evidence,
      profiles = c("TARGET_ONLY", "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT"),
      ontology = ont,
      min_shared_genes = min_shared_genes,
      min_gene_count = min_gene_count,
      max_terms = 500L
    )
    readr::write_csv(dat$nodes,
                     file.path(compatibility_root, paste0("primary_network_nodes_", ont, ".csv")))
    readr::write_csv(dat$edges,
                     file.path(compatibility_root, paste0("primary_network_edges_", ont, ".csv")))
  }
  readr::write_csv(summary,
                   file.path(compatibility_root, "primary_network_summary.csv"))
  invisible(file.path(outdir, "networks"))
}

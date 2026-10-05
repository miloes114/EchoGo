# Canonical scoreless network products ---------------------------------------
#
# Networks are pipeline products. The report should consume these assets rather
# than silently rebuilding an analysis during knitting.

.echogo_network_product_definitions <- function() {
  list(
    target_supported = c("TARGET_ONLY", "TARGET_PLUS_CONTEXT"),
    alternative_context_hypothesis = "ALTERNATIVE_CONTEXT"
  )
}

.echogo_write_network_graphml <- function(nodes, edges, path) {
  if (!requireNamespace("igraph", quietly = TRUE) || !nrow(nodes) || !nrow(edges)) return(NA_character_)
  graph <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  igraph::write_graph(graph, path, format = "graphml")
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

.echogo_write_network_html <- function(nodes, edges, path) {
  if (!requireNamespace("visNetwork", quietly = TRUE) ||
      !requireNamespace("htmlwidgets", quietly = TRUE) ||
      !nrow(nodes) || !nrow(edges)) {
    return(NA_character_)
  }

  vis_nodes <- nodes |>
    dplyr::transmute(
      id = .data$term_id,
      label = dplyr::if_else(
        is.na(.data$term_name) | !nzchar(.data$term_name),
        .data$term_id,
        .data$term_name
      ),
      value = .data$gene_count,
      group = .data$evidence_profile,
      title = paste0(
        "<b>", .data$term_id, "</b><br/>",
        dplyr::coalesce(.data$term_name, ""), "<br/>",
        "Evidence: ", .data$evidence_profile, "<br/>",
        "Contributing genes: ", .data$gene_count, "<br/>",
        "Other references recovering term: ",
        .data$alternative_context_support_n, "/", .data$alternative_queried_context_n
      )
    )
  vis_edges <- edges |>
    dplyr::transmute(
      from = .data$from,
      to = .data$to,
      value = .data$shared_gene_n,
      title = paste0("Shared contributing genes: ", .data$shared_gene_n)
    )

  widget <- visNetwork::visNetwork(vis_nodes, vis_edges, width = "100%", height = "760px") |>
    visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) |>
    visNetwork::visPhysics(stabilization = TRUE)

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  htmlwidgets::saveWidget(widget, path, selfcontained = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

#' Write the two canonical scoreless gene-overlap network products
#'
#' @param evidence Canonical exact-term EchoGO evidence.
#' @param outdir EchoGO run root.
#' @param min_shared_genes Minimum shared contributing genes required for an edge.
#' @param min_gene_count Minimum contributing genes required for a term to enter the graph.
#' @param max_terms Maximum terms per product and ontology used for graph construction.
#' @param label_per_community Maximum labels per graph community in the static plot.
#' @return Invisible path to the canonical network directory.
#' @keywords internal
.echogo_write_scoreless_network_products <- function(
    evidence,
    outdir,
    min_shared_genes = 2L,
    min_gene_count = 2L,
    max_terms = 120L,
    label_per_community = 3L
) {
  root <- file.path(outdir, "networks")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  status_rows <- list()
  product_definitions <- .echogo_network_product_definitions()

  for (product in names(product_definitions)) {
    profiles <- product_definitions[[product]]
    for (ont in c("BP", "MF", "CC")) {
      product_dir <- file.path(root, product, ont)
      dir.create(product_dir, recursive = TRUE, showWarnings = FALSE)

      dat <- .echogo_report_network_data(
        evidence = evidence,
        profiles = profiles,
        ontology = ont,
        min_shared_genes = min_shared_genes,
        min_gene_count = min_gene_count,
        max_terms = max_terms
      )
      nodes <- dat$nodes
      edges <- dat$edges
      diagnostics <- dat$diagnostics
      # Avoid dplyr data-mask collision with the diagnostics `product` column.
      diagnostics$semantic_product <- product
      diagnostics$status_message <- .echogo_network_skip_message(dat$diagnostics)

      node_path <- file.path(product_dir, "nodes.csv")
      edge_path <- file.path(product_dir, "edges.csv")
      readr::write_csv(nodes, node_path)
      readr::write_csv(edges, edge_path)

      graphml_path <- .echogo_write_network_graphml(
        nodes,
        edges,
        file.path(product_dir, "network.graphml")
      )

      png_path <- pdf_path <- svg_path <- html_path <- NA_character_
      render_error <- NA_character_
      if (identical(diagnostics$status[[1]], "GENERATED")) {
        p <- tryCatch(
          plot_echogo_gene_overlap_network(
            evidence = evidence,
            profiles = profiles,
            ontology = ont,
            min_shared_genes = min_shared_genes,
            min_gene_count = min_gene_count,
            max_terms = max_terms,
            label_per_community = label_per_community
          ),
          error = function(e) {
            render_error <<- conditionMessage(e)
            NULL
          }
        )
        if (!is.null(p)) {
          variants <- tryCatch(
            .echogo_save_plot_variants(
              p,
              directory = product_dir,
              stem = "network",
              width = 13,
              height = 10
            ),
            error = function(e) {
              render_error <<- conditionMessage(e)
              list(png = NA_character_, pdf = NA_character_, svg = NA_character_)
            }
          )
          png_path <- variants$png
          pdf_path <- variants$pdf
          svg_path <- variants$svg
          html_path <- tryCatch(
            .echogo_write_network_html(nodes, edges, file.path(product_dir, "network.html")),
            error = function(e) {
              render_error <<- conditionMessage(e)
              NA_character_
            }
          )
        }
      }

      if (!is.na(render_error) && nzchar(render_error)) {
        diagnostics$status[[1]] <- "ERROR_RENDERING"
        diagnostics$status_message[[1]] <- paste0("Network data were generated, but rendering failed: ", render_error)
      }

      status_row <- diagnostics |>
        dplyr::mutate(
          nodes_file = normalizePath(node_path, winslash = "/", mustWork = FALSE),
          edges_file = normalizePath(edge_path, winslash = "/", mustWork = FALSE),
          graphml_file = graphml_path,
          png_file = if (!is.na(png_path)) normalizePath(png_path, winslash = "/", mustWork = FALSE) else NA_character_,
          pdf_file = if (!is.na(pdf_path)) normalizePath(pdf_path, winslash = "/", mustWork = FALSE) else NA_character_,
          svg_file = if (!is.na(svg_path)) normalizePath(svg_path, winslash = "/", mustWork = FALSE) else NA_character_,
          html_file = html_path
        )
      readr::write_csv(status_row, file.path(product_dir, "status.csv"))
      status_rows[[paste(product, ont, sep = "_")]] <- status_row
    }
  }

  readr::write_csv(dplyr::bind_rows(status_rows), file.path(root, "network_product_status.csv"))
  invisible(root)
}

.echogo_read_network_status <- function(network_dir) {
  path <- file.path(network_dir, "network_product_status.csv")
  if (!file.exists(path)) return(tibble::tibble())
  tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) tibble::tibble())
}

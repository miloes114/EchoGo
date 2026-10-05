test_that("RRvGO similarity guard handles malformed small-ontology outputs", {
  invalid_objects <- list(
    NULL,
    numeric(),
    NA_real_,
    matrix(NA_real_, nrow = 2, ncol = 2),
    matrix(1, nrow = 1, ncol = 2)
  )

  for (simMatrix in invalid_objects) {
    expect_true(.echogo_invalid_similarity_matrix(simMatrix))
  }

  expect_false(.echogo_invalid_similarity_matrix(diag(2)))
})

test_that("RRvGO attaches prepared primary provenance without caller internals", {
  skip_if_not_installed("rrvgo")
  skip_if_not_installed("GO.db")
  skip_if_not_installed("org.Dr.eg.db")

  # This is a normal scoreless evidence input: callers do not know about, and
  # must not supply, the RRvGO-private provenance field.
  input <- tibble::tibble(
    term_id = c("GO:0006355", "GO:0006412", "GO:0006950"),
    term_name = c("regulation of transcription", "translation", "response to stress"),
    ontology = "BP",
    evidence_profile = c("TARGET_ONLY", "ALTERNATIVE_CONTEXT", "TARGET_ONLY"),
    primary_evidence = TRUE,
    representative_order = 1:3,
    display_order = 1:3
  )
  expect_false(".echogo_rrvgo_origin" %in% names(input))

  prepared <- .echogo_prepare_rrvgo_primary_terms(input)
  expect_true(".echogo_rrvgo_origin" %in% names(prepared))
  expect_identical(prepared$.echogo_rrvgo_origin, input$evidence_profile)

  output <- tempfile("rrvgo-prepared-provenance-")
  run_rrvgo_consensus_analysis(
    df_input = input,
    label = "provenance_regression",
    output_base = output,
    semantic_reference_orgdb = "org.Dr.eg.db",
    semantic_reference_role = "target_reference",
    ontologies = "BP"
  )
  clusters <- list.files(output, pattern = "rrvgo_BP_clusters\\.csv$", recursive = TRUE, full.names = TRUE)
  expect_length(clusters, 1L)
  reduced <- readr::read_csv(clusters, show_col_types = FALSE)
  expect_true("origin" %in% names(reduced))
  expect_true(all(reduced$origin %in% unique(input$evidence_profile)))
  expect_identical(unique(reduced$semantic_method), "Rel")
  expect_identical(unique(reduced$semantic_reference_orgdb), "org.Dr.eg.db")
  expect_identical(unique(reduced$semantic_reference_role), "target_reference")
  expect_identical(unique(reduced$semantic_product), "target_supported")
  expect_true(all(reduced$go %in% input$term_id))
})

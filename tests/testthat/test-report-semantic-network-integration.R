test_that("scoreless network renderer restores term_id after igraph normalization", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")

  evidence <- tibble::tibble(
    term_id = c("GO:0000001", "GO:0000002", "GO:0000003"),
    term_name = c("term one", "term two", "term three"),
    evidence_profile = c("TARGET_ONLY", "TARGET_PLUS_CONTEXT", "TARGET_PLUS_CONTEXT"),
    ontology = "BP",
    alternative_context_support_n = c(0L, 2L, 1L),
    alternative_queried_context_n = 3L,
    alternative_context_support_fraction = c(0, 2/3, 1/3),
    display_order = 1:3,
    contributing_genes = c("A;B;C", "B;C;D", "B;C;E")
  )

  p <- plot_echogo_gene_overlap_network(
    evidence,
    profiles = c("TARGET_ONLY", "TARGET_PLUS_CONTEXT"),
    ontology = "BP",
    min_shared_genes = 2L,
    min_gene_count = 2L
  )

  expect_s3_class(p, "ggplot")
  diag <- attr(p, "echogo_network_diagnostics")
  expect_equal(diag$status[[1]], "GENERATED")
  expect_equal(diag$final_edges[[1]], 3L)
  expect_true(all(attr(p, "echogo_network_nodes")$term_id %in% evidence$term_id))
})

test_that("network diagnostics distinguish legitimate sparse skips", {
  evidence <- tibble::tibble(
    term_id = c("GO:0000001", "GO:0000002"),
    term_name = c("term one", "term two"),
    evidence_profile = c("ALTERNATIVE_CONTEXT", "ALTERNATIVE_CONTEXT"),
    ontology = "MF",
    alternative_context_support_n = c(1L, 1L),
    alternative_queried_context_n = 3L,
    alternative_context_support_fraction = c(1/3, 1/3),
    display_order = 1:2,
    contributing_genes = c("A;B", "C;D")
  )

  dat <- EchoGO:::.echogo_report_network_data(
    evidence,
    profiles = "ALTERNATIVE_CONTEXT",
    ontology = "MF",
    min_shared_genes = 2L,
    min_gene_count = 2L
  )
  expect_equal(dat$diagnostics$status[[1]], "SKIPPED_NO_SHARED_GENES")
  expect_equal(dat$diagnostics$pairs_sharing_at_least_one_gene[[1]], 0L)
  expect_match(EchoGO:::.echogo_network_skip_message(dat$diagnostics), "no pair shared")
})

test_that("matching RRvGO sidecars can be discovered without missing compatibility metadata", {
  tmp <- tempfile("echogo-rrvgo-")
  dir.create(file.path(tmp, "rrvgo_target_supported", "OrgDb=org.Test.eg.db"), recursive = TRUE)
  path <- file.path(tmp, "rrvgo_target_supported", "OrgDb=org.Test.eg.db", "rrvgo_BP_clusters.csv")
  readr::write_csv(
    tibble::tibble(
      go = c("GO:0000001", "GO:0000002"),
      term = c("term one", "term two"),
      semantic_product = "target_supported",
      semantic_reference_orgdb = "org.Test.eg.db",
      semantic_reference_role = "proxy",
      semantic_method = "Rel"
    ),
    path
  )

  evidence <- tibble::tibble(
    term_id = c("GO:0000001", "GO:0000002", "GO:0000003"),
    evidence_profile = c("TARGET_ONLY", "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT")
  )

  found <- EchoGO:::.echogo_discover_rrvgo_products(
    tmp,
    evidence,
    configuration = list(run_rrvgo = FALSE)
  )

  expect_true(found$enabled)
  expect_equal(found$state, "DISCOVERED_MATCHING_SIDECARS")
  expect_equal(found$semantic_reference_orgdb, "org.Test.eg.db")
  expect_equal(found$semantic_reference_role, "proxy")
  expect_equal(found$semantic_method, "Rel")
  expect_equal(nrow(found$valid_products), 1L)
})

test_that("RRvGO sidecars from a different evidence set are rejected", {
  tmp <- tempfile("echogo-rrvgo-mismatch-")
  dir.create(file.path(tmp, "rrvgo_target_supported", "OrgDb=org.Test.eg.db"), recursive = TRUE)
  path <- file.path(tmp, "rrvgo_target_supported", "OrgDb=org.Test.eg.db", "rrvgo_BP_clusters.csv")
  readr::write_csv(
    tibble::tibble(
      go = c("GO:9999999"),
      term = "stale term",
      semantic_product = "target_supported",
      semantic_reference_orgdb = "org.Test.eg.db",
      semantic_reference_role = "proxy",
      semantic_method = "Rel"
    ),
    path
  )

  evidence <- tibble::tibble(
    term_id = c("GO:0000001"),
    evidence_profile = "TARGET_ONLY"
  )

  found <- EchoGO:::.echogo_discover_rrvgo_products(
    tmp,
    evidence,
    configuration = list(run_rrvgo = FALSE)
  )

  expect_false(found$enabled)
  expect_equal(nrow(found$valid_products), 0L)
  expect_equal(found$products$status[[1]], "REJECTED_EVIDENCE_MISMATCH")
})

test_that("RRvGO status manifest round-trips explicit generated and skipped states", {
  tmp <- tempfile("echogo-rrvgo-status-")
  rows <- list(
    EchoGO:::.echogo_rrvgo_status_row(
      "target_supported", "BP", "GENERATED",
      input_term_n = 12L, valid_term_n = 12L,
      semantic_reference_orgdb = "org.Test.eg.db",
      semantic_reference_role = "target_reference"
    ),
    EchoGO:::.echogo_rrvgo_status_row(
      "alternative_context_hypothesis", "MF", "SKIPPED_TOO_FEW_TERMS",
      reason = "Only one valid GO term was available.",
      input_term_n = 1L, valid_term_n = 1L,
      semantic_reference_orgdb = "org.Test.eg.db",
      semantic_reference_role = "target_reference"
    )
  )
  path <- EchoGO:::.echogo_write_rrvgo_status(rows, tmp)
  expect_true(file.exists(path))
  tab <- EchoGO:::.echogo_read_rrvgo_status(tmp)
  expect_equal(nrow(tab), 2L)
  expect_equal(tab$status, c("GENERATED", "SKIPPED_TOO_FEW_TERMS"))
  expect_match(EchoGO:::.echogo_rrvgo_status_message(tab[2, ]), "insufficient")
})

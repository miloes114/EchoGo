test_that("ontology normalization is explicit and stable", {
  expect_equal(
    .echogo_normalize_ontology_code(c(
      "BP", "GO:BP", "biological process", "Biological_Process",
      "MF", "GO:MF", "molecular function",
      "CC", "GO:CC", "cellular-component"
    )),
    c("BP", "BP", "BP", "BP", "MF", "MF", "MF", "CC", "CC", "CC")
  )
  expect_true(is.na(.echogo_normalize_ontology_code("unknown")))
})

.report_scope_fixture <- function() {
  tibble::tibble(
    term_id = c("GO:BP1", "GO:BP2", "GO:MF1", "GO:MF2", "GO:CC1", "GO:CC2"),
    term_name = c(
      "biological process one", "biological process two",
      "molecular function one", "molecular function two",
      "cellular component one", "cellular component two"
    ),
    ontology = c("BP", "GO:BP", "MF", "GO:MF", "CC", "GO:CC"),
    evidence_profile = c(
      "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT",
      "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT",
      "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT"
    ),
    target_goseq_supported = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
    target_goseq_adjusted_p = c(.001, NA, .002, NA, .003, NA),
    alternative_context_support_n = c(2L, 1L, 2L, 1L, 2L, 1L),
    alternative_queried_context_n = rep(3L, 6),
    alternative_context_support_fraction = c(2/3, 1/3, 2/3, 1/3, 2/3, 1/3),
    contributing_genes = c(
      "g1;g2;g3", "g1;g2",
      "m1;m2;m3", "m1;m2",
      "c1;c2;c3", "c1;c2"
    ),
    display_order = seq_len(6)
  )
}

test_that("Key Findings input is strictly ontology-specific", {
  ev <- .report_scope_fixture()

  bp <- .echogo_key_findings_data(ev, ontology = "BP", max_target_recovered = 10, max_hypotheses = 10)
  mf <- .echogo_key_findings_data(ev, ontology = "MF", max_target_recovered = 10, max_hypotheses = 10)
  cc <- .echogo_key_findings_data(ev, ontology = "CC", max_target_recovered = 10, max_hypotheses = 10)

  expect_setequal(bp$term_id, c("GO:BP1", "GO:BP2"))
  expect_setequal(mf$term_id, c("GO:MF1", "GO:MF2"))
  expect_setequal(cc$term_id, c("GO:CC1", "GO:CC2"))
  expect_length(intersect(bp$term_id, mf$term_id), 0)
  expect_length(intersect(bp$term_id, cc$term_id), 0)
  expect_length(intersect(mf$term_id, cc$term_id), 0)
})

test_that("Evidence Landscape term selection is strictly ontology-specific", {
  ev <- .report_scope_fixture()

  bp <- .echogo_landscape_terms(ev, "BP", max_terms_per_profile = 20)
  mf <- .echogo_landscape_terms(ev, "MF", max_terms_per_profile = 20)
  cc <- .echogo_landscape_terms(ev, "CC", max_terms_per_profile = 20)

  expect_setequal(bp$term_id, c("GO:BP1", "GO:BP2"))
  expect_setequal(mf$term_id, c("GO:MF1", "GO:MF2"))
  expect_setequal(cc$term_id, c("GO:CC1", "GO:CC2"))
  expect_true(all(.echogo_normalize_ontology_code(bp$ontology) == "BP"))
  expect_true(all(.echogo_normalize_ontology_code(mf$ontology) == "MF"))
  expect_true(all(.echogo_normalize_ontology_code(cc$ontology) == "CC"))
})

test_that("network products cannot leak terms across ontologies", {
  ev <- .report_scope_fixture()

  bp <- .echogo_report_network_data(
    ev,
    profiles = c("TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT"),
    ontology = "BP",
    min_shared_genes = 2,
    min_gene_count = 2,
    max_terms = 20
  )
  mf <- .echogo_report_network_data(
    ev,
    profiles = c("TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT"),
    ontology = "MF",
    min_shared_genes = 2,
    min_gene_count = 2,
    max_terms = 20
  )
  cc <- .echogo_report_network_data(
    ev,
    profiles = c("TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT"),
    ontology = "CC",
    min_shared_genes = 2,
    min_gene_count = 2,
    max_terms = 20
  )

  expect_setequal(bp$nodes$term_id, c("GO:BP1", "GO:BP2"))
  expect_setequal(mf$nodes$term_id, c("GO:MF1", "GO:MF2"))
  expect_setequal(cc$nodes$term_id, c("GO:CC1", "GO:CC2"))
  expect_equal(bp$diagnostics$input_terms, 2L)
  expect_equal(mf$diagnostics$input_terms, 2L)
  expect_equal(cc$diagnostics$input_terms, 2L)
  expect_equal(bp$diagnostics$final_edges, 1L)
  expect_equal(mf$diagnostics$final_edges, 1L)
  expect_equal(cc$diagnostics$final_edges, 1L)
})

test_that("unsupported ontology labels fail rather than widening the selection", {
  ev <- .report_scope_fixture()
  expect_error(.echogo_subset_evidence_ontology(ev, "ALL"), "Expected BP, MF, or CC")
})

test_that("report profile labels remain plain-language views of scoreless evidence", {
  evidence <- tibble::tibble(
    term_id = c("GO:0000001", "GO:0000002", "GO:0000003"),
    term_name = c("alpha process", "beta process", "gamma process"),
    ontology = "BP",
    evidence_profile = c("TARGET_ONLY", "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT"),
    target_goseq_supported = c(TRUE, TRUE, FALSE),
    alternative_context_support_n = c(0, 1, 1),
    alternative_queried_context_n = c(2, 2, 2),
    alternative_context_support_fraction = c(0, .5, .5),
    display_order = 1:3,
    contributing_genes = c("a;b", "b;c", "c;d")
  )

  p <- EchoGO:::plot_echogo_profile_counts(evidence)
  expect_s3_class(p, "ggplot")
  expect_false(any(grepl("consensus_score|EQI", c(p$labels$title, p$labels$x, p$labels$y), ignore.case = TRUE)))
})

test_that("Evidence Landscape is a direct view of source qualification", {
  contexts <- tibble::tibble(
    context_code = c("target", "ref2", "ref3"),
    context_label = c("target_reference", "reference_two", "reference_three"),
    context_role = c("TARGET", "ALTERNATIVE", "ALTERNATIVE")
  )
  evidence <- tibble::tibble(
    term_id = c("GO:0000001", "GO:0000002", "GO:0000003"),
    term_name = c("alpha process", "beta process", "gamma process"),
    ontology = "BP",
    evidence_profile = c("TARGET_ONLY", "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT"),
    target_goseq_supported = c(TRUE, TRUE, FALSE),
    alternative_context_support_n = c(0, 1, 1),
    alternative_queried_context_n = c(2, 2, 2),
    alternative_context_support_fraction = c(0, .5, .5),
    display_order = 1:3
  )
  sources <- tibble::tibble(
    term_id = c("GO:0000002", "GO:0000003", "GO:0000003"),
    source_type = "gprofiler",
    background_mode = "custom_experimental_background",
    context_label = c("reference_two", "reference_two", "reference_three"),
    context_role = "ALTERNATIVE",
    source_qualifies = c(TRUE, TRUE, FALSE)
  )

  dat <- EchoGO:::.echogo_landscape_long(evidence, sources, contexts, "BP", 12)
  expect_equal(nrow(dat$matrix), nrow(evidence) * (nrow(contexts) + 1L))

  beta_ref2 <- dat$matrix[
    dat$matrix$term_id == "GO:0000002" & dat$matrix$source_column == "reference_two",
    , drop = FALSE
  ]
  gamma_ref3 <- dat$matrix[
    dat$matrix$term_id == "GO:0000003" & dat$matrix$source_column == "reference_three",
    , drop = FALSE
  ]
  expect_true(beta_ref2$qualifies)
  expect_false(gamma_ref3$qualifies)
  expect_equal(
    unique(dat$matrix$alternative_context_support_n[dat$matrix$term_id == "GO:0000002"]),
    1
  )
})

test_that("gProfiler reader is not shadowed by an empty canonical plot directory", {
  td <- tempfile("echogo-gp-")
  dir.create(td, recursive = TRUE)
  canonical <- file.path(td, "custom_experimental_background")
  legacy <- file.path(td, "with_custom_background")
  dir.create(canonical)
  dir.create(legacy)

  file.create(file.path(canonical, "gprofiler_reference_two_custom_experimental_background_GO_BP_lollipop.pdf"))

  readr::write_csv(
    tibble::tibble(
      term_id = "GO:0000002",
      term_name = "beta process",
      source = "GO:BP",
      p_value = 0.01,
      fold_enrichment = 2,
      intersection = "a,b",
      intersection_size = 2
    ),
    file.path(legacy, "gprofiler_reference_two_with_bg.csv")
  )

  species_map <- c(ref2 = "reference_two")
  out <- EchoGO:::.echogo_read_gprofiler_long(
    td,
    species_map,
    mode = "custom_experimental_background",
    target_context = NULL
  )

  expect_equal(nrow(out), 1L)
  expect_equal(out$term_id, "GO:0000002")
  expect_match(out$provenance_file, "with_custom_background", fixed = TRUE)
})

test_that("report network products keep target and hypothesis profiles separate", {
  evidence <- tibble::tibble(
    term_id = c("GO:0000001", "GO:0000002", "GO:0000003"),
    term_name = c("alpha", "beta", "gamma"),
    ontology = "BP",
    evidence_profile = c("TARGET_ONLY", "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT"),
    alternative_context_support_n = c(0, 1, 1),
    alternative_queried_context_n = c(2, 2, 2),
    alternative_context_support_fraction = c(0, .5, .5),
    display_order = 1:3,
    contributing_genes = c("a;b;c", "b;c;d", "b;c;e")
  )

  target <- EchoGO:::.echogo_report_network_data(
    evidence, c("TARGET_ONLY", "TARGET_PLUS_CONTEXT"), "BP", min_shared_genes = 1
  )
  hyp <- EchoGO:::.echogo_report_network_data(
    evidence, "ALTERNATIVE_CONTEXT", "BP", min_shared_genes = 1
  )

  expect_true(all(target$nodes$evidence_profile %in% c("TARGET_ONLY", "TARGET_PLUS_CONTEXT")))
  expect_true(all(hyp$nodes$evidence_profile == "ALTERNATIVE_CONTEXT"))
  expect_length(intersect(target$nodes$term_id, hyp$nodes$term_id), 0L)
})
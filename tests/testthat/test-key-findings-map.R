test_that("key findings map uses transparent profile-specific display ordering", {
  evidence <- tibble::tibble(
    term_id = c("GO:0000001", "GO:0000002", "GO:0000003", "GO:0000004", "GO:0000005", "GO:0000006"),
    term_name = c("target broad", "target narrow", "target only strong", "target only weak", "hyp broad", "hyp narrow"),
    ontology = "BP",
    evidence_profile = c(
      "TARGET_PLUS_CONTEXT", "TARGET_PLUS_CONTEXT",
      "TARGET_ONLY", "TARGET_ONLY",
      "ALTERNATIVE_CONTEXT", "ALTERNATIVE_CONTEXT"
    ),
    target_goseq_adjusted_p = c(0.03, 0.001, 0.002, 0.04, NA, NA),
    alternative_context_support_n = c(5L, 2L, 0L, 0L, 6L, 2L),
    alternative_queried_context_n = rep(8L, 6),
    alternative_context_support_fraction = c(5, 2, 0, 0, 6, 2) / 8,
    contributing_genes = c("a;b;c", "d;e", "f", "g;h", "i;j;k", "l"),
    display_order = 1:6
  )

  out <- .echogo_key_findings_data(
    evidence,
    ontology = "BP",
    max_target_recovered = 2,
    max_target_only = 2,
    max_hypotheses = 2
  )

  expect_equal(
    out$term_id,
    c("GO:0000002", "GO:0000001", "GO:0000003", "GO:0000004", "GO:0000005", "GO:0000006")
  )
  expect_equal(out$contributing_gene_n, c(2L, 3L, 1L, 2L, 3L, 1L))
  expect_false(any(grepl("score", names(out), ignore.case = TRUE)))
})

test_that("alternative hypotheses are not ranked by gProfiler p-values", {
  evidence <- tibble::tibble(
    term_id = c("GO:1000001", "GO:1000002"),
    term_name = c("more references", "fewer references"),
    ontology = "MF",
    evidence_profile = "ALTERNATIVE_CONTEXT",
    target_goseq_adjusted_p = NA_real_,
    alternative_context_support_n = c(5L, 2L),
    alternative_queried_context_n = c(8L, 8L),
    alternative_context_support_fraction = c(5 / 8, 2 / 8),
    contributing_genes = c("a;b", "c"),
    display_order = c(2L, 1L)
  )

  out <- .echogo_key_findings_data(evidence, ontology = "MF", max_hypotheses = 2)
  expect_equal(out$term_id, c("GO:1000001", "GO:1000002"))
  expect_equal(out$support_label, c("5/8", "2/8"))
})

test_that("NO_PRIMARY_SUPPORT never enters the key findings map", {
  evidence <- tibble::tibble(
    term_id = c("GO:2000001", "GO:2000002"),
    term_name = c("kept", "excluded"),
    ontology = "CC",
    evidence_profile = c("TARGET_ONLY", "NO_PRIMARY_SUPPORT"),
    target_goseq_adjusted_p = c(0.01, NA),
    alternative_context_support_n = c(0L, 0L),
    alternative_queried_context_n = c(3L, 3L),
    alternative_context_support_fraction = c(0, 0),
    contributing_genes = c("a", "b"),
    display_order = 1:2
  )

  out <- .echogo_key_findings_data(evidence, ontology = "CC")
  expect_equal(out$term_id, "GO:2000001")
})

test_that("key findings map renders with literal recurrence labels", {
  evidence <- tibble::tibble(
    term_id = c("GO:3000001", "GO:3000002", "GO:3000003"),
    term_name = c("target recovered", "target only", "new hypothesis"),
    ontology = "BP",
    evidence_profile = c("TARGET_PLUS_CONTEXT", "TARGET_ONLY", "ALTERNATIVE_CONTEXT"),
    target_goseq_adjusted_p = c(0.01, 0.02, NA),
    alternative_context_support_n = c(3L, 0L, 2L),
    alternative_queried_context_n = c(4L, 4L, 4L),
    alternative_context_support_fraction = c(.75, 0, .5),
    contributing_genes = c("a;b", "c", "d;e;f"),
    display_order = 1:3
  )

  p <- plot_echogo_key_findings_map(evidence, ontology = "BP")
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  labels <- unlist(lapply(built$data, function(x) if ("label" %in% names(x)) as.character(x$label) else character()))
  expect_true(all(c("3/4", "0/4", "2/4") %in% labels))
})

test_that("theme summary reports reference contribution without cluster scoring", {
  reduced <- tibble::tibble(
    go = c("GO:4000001", "GO:4000002", "GO:4000003"),
    term = c("theme one", "member one", "theme two"),
    cluster = c("GO:4000001", "GO:4000001", "GO:4000003"),
    score = c(-1, -2, -3)
  )
  evidence <- tibble::tibble(
    term_id = reduced$go,
    alternative_queried_context_n = 4L,
    contributing_genes = c("a;b", "b;c", "d"),
    evidence_profile = c("TARGET_PLUS_CONTEXT", "TARGET_PLUS_CONTEXT", "TARGET_ONLY")
  )
  sources <- tibble::tibble(
    term_id = c("GO:4000001", "GO:4000002", "GO:4000002"),
    source_type = "gprofiler",
    background_mode = "custom_experimental_background",
    context_role = "ALTERNATIVE",
    context_label = c("ref1", "ref2", "ref3"),
    source_qualifies = TRUE
  )

  out <- .echogo_theme_summary_data(reduced, evidence, sources)
  expect_equal(out$reference_support_n[out$cluster == "GO:4000001"], 3L)
  expect_equal(out$exact_term_n[out$cluster == "GO:4000001"], 2L)
  expect_equal(out$representative_term[out$cluster == "GO:4000001"], "theme one")
  expect_false(any(grepl("consensus|integrated|evidence_score", names(out), ignore.case = TRUE)))
})

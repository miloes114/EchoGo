test_that("the report exposes the target-anchored interpretation model", {
  report <- system.file("reports", "echogo_report.Rmd", package = "EchoGO")
  if (!nzchar(report)) {
    root_source <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/")
    report <- file.path(root_source, "inst", "reports", "echogo_report.Rmd")
  }
  expect_true(file.exists(report))
  report <- paste(readLines(report, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(report, "Experiment and representability")
  expect_match(report, "What the target analysis already says")
  expect_match(report, "New alternative-context hypotheses")
  expect_match(report, "Target-supported semantic themes")
  expect_match(report, "Alternative-context hypothesis semantic themes")
  expect_match(report, "Hypothesis drill-down")
  expect_match(report, "Default-domain exploratory hypotheses")
})

test_that("the offline demo uses explicit inputs without retired score-era results", {
  demo <- system.file("extdata", "echogo_demo", package = "EchoGO")
  cache <- file.path(demo, "cache")
  manifest <- jsonlite::read_json(file.path(cache, "run_manifest.json"), simplifyVector = TRUE)

  expect_true(dir.exists(cache))
  expect_identical(manifest$echogo_version, "0.1.4")
  expect_false(dir.exists(file.path(demo, "cached_gprofiler_v0.1.3")))
  expect_true(file.exists(file.path(echogo_demo_results_path(), "README.txt")))
  expect_false(nzchar(system.file("extdata", "echogo_demo_results_v0.1.3", package = "EchoGO")))
})

test_that("report template identifies default-domain evidence without changing primary profiles", {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) skip("rmarkdown is not installed")
  report_template <- system.file("reports", "echogo_report.Rmd", package = "EchoGO")
  if (!nzchar(report_template)) {
    root_source <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/")
    report_template <- file.path(root_source, "inst", "reports", "echogo_report.Rmd")
  }
  root <- tempfile("echogo-report-contract-")
  dirs <- list(base = root, consensus = file.path(root, "consensus"), rrvgo = file.path(root, "rrvgo"))
  dir.create(dirs$consensus, recursive = TRUE)
  dir.create(file.path(root, "gprofiler", "submitted_vectors"), recursive = TRUE)
  dir.create(file.path(root, "diagnostics"), recursive = TRUE)
  readr::write_csv(data.frame(
    term_id = "GO:0000001", term_name = "example term", ontology = "BP",
    evidence_profile = "ALTERNATIVE_CONTEXT", target_goseq_supported = FALSE,
    alternative_context_support_n = 1, alternative_context_support_fraction = 1,
    recovered_alternative_contexts = "mouse", contributing_genes = "GENE1", display_order = 1
  ), file.path(dirs$consensus, "term_evidence_exact.csv"))
  readr::write_csv(data.frame(
    term_id = "GO:0000001", term_name = "example term", ontology = "BP",
    source_type = "gprofiler", context_label = "mouse", context_role = "ALTERNATIVE",
    context_rationale = "example", background_mode = "default_domain_exploratory",
    domain_scope = "default_domain", source_adjusted_p_value = 0.01,
    fold_enrichment = 2, effective_query_n = 2, effective_background_or_domain_n = 100
  ), file.path(dirs$consensus, "term_source_provenance_long.csv"))
  readr::write_csv(data.frame(tested = TRUE, significant = TRUE, resolved_name = "GENE1", included_background = TRUE, included_foreground = TRUE, annotation_match = TRUE), file.path(root, "gprofiler", "submitted_vectors", "mapping_table.csv"))
  jsonlite::write_json(list(run_rrvgo = FALSE, run_exploratory_default_domain = TRUE), file.path(root, "diagnostics", "run_configuration.json"), auto_unbox = TRUE)
  out <- rmarkdown::render(report_template, output_dir = root, params = list(dirs = dirs), quiet = TRUE)
  html <- paste(readLines(out, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(html, "annotated statistical domain")
  expect_match(html, "ALTERNATIVE_CONTEXT")
})

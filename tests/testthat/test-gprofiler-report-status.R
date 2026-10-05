test_that("gProfiler report status distinguishes no source rows from no qualifying terms", {
  src <- tibble::tibble(
    source_type = c("gprofiler", "gprofiler"),
    background_mode = c("custom_experimental_background", "custom_experimental_background"),
    context_label = c("refA", "refA"),
    ontology = c("BP", "BP"),
    term_id = c("GO:0000001", "GO:0000002"),
    source_qualifies = c(FALSE, FALSE)
  )

  no_terms <- EchoGO:::.echogo_gprofiler_report_status(src, "refA", "BP")
  expect_identical(no_terms$status[[1]], "NO_QUALIFYING_TERMS")
  expect_equal(no_terms$qualifying_term_n[[1]], 0L)

  missing <- EchoGO:::.echogo_gprofiler_report_status(src, "refB", "BP")
  expect_identical(missing$status[[1]], "NO_SOURCE_ROWS")
})

test_that("gProfiler report status requires a figure when qualifying terms exist", {
  src <- tibble::tibble(
    source_type = "gprofiler",
    background_mode = "custom_experimental_background",
    context_label = "refA",
    ontology = "MF",
    term_id = "GO:0000001",
    source_qualifies = TRUE
  )
  status <- EchoGO:::.echogo_gprofiler_report_status(src, "refA", "MF")
  expect_identical(status$status[[1]], "FIGURE_EXPECTED")
  expect_equal(status$qualifying_term_n[[1]], 1L)
})

test_that("default-domain results use exploratory rather than qualifying status language", {
  src <- tibble::tibble(
    source_type = "gprofiler",
    background_mode = "default_domain_exploratory",
    context_label = "refA",
    ontology = "GO:BP",
    term_id = c("GO:0000001", "GO:0000002"),
    source_qualifies = c(TRUE, FALSE)
  )
  status <- EchoGO:::.echogo_gprofiler_report_status(src, "refA", "BP", exploratory = TRUE)
  expect_identical(status$status[[1]], "EXPLORATORY_PLOT_EXPECTED")
  expect_equal(status$exploratory_term_n[[1]], 2L)
  expect_true(is.na(status$qualifying_term_n[[1]]))
  expect_match(status$message[[1]], "default-domain exploratory")
  expect_false(grepl("qualifying|significant|matched-background", status$message[[1]], ignore.case = TRUE))

  empty <- src[FALSE, ]
  empty_status <- EchoGO:::.echogo_gprofiler_report_status(empty, "refA", "BP", exploratory = TRUE)
  expect_identical(empty_status$status[[1]], "NO_SOURCE_ROWS")
  expect_false(grepl("qualifying|significant|matched-background", empty_status$message[[1]], ignore.case = TRUE))
})

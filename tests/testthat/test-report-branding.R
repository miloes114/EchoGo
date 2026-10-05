test_that("first-screen report cards summarize only frozen evidence profiles", {
  evidence <- tibble::tibble(
    evidence_profile = c(
      "TARGET_ONLY", "TARGET_PLUS_CONTEXT", "TARGET_PLUS_CONTEXT",
      "ALTERNATIVE_CONTEXT", "NO_PRIMARY_SUPPORT"
    ),
    alternative_queried_context_n = 3L
  )
  contexts <- tibble::tibble(
    context_label = c("target", "ref1", "ref2", "ref3"),
    context_role = c("TARGET", "ALTERNATIVE", "ALTERNATIVE", "ALTERNATIVE")
  )

  out <- .echogo_primary_metric_cards(evidence, contexts)
  expect_equal(out$value, c(3L, 2L, 1L, 1L, 3L))
  expect_equal(
    out$label,
    c(
      "GO terms supported by Target GOseq",
      "Target + context",
      "Target only",
      "Context-derived hypotheses",
      "Selected annotation contexts"
    )
  )
  expect_false(any(grepl("score|significance|confidence", out$label, ignore.case = TRUE)))
})

test_that("representability statement is descriptive and derived from mapping audit", {
  mapping <- tibble::tibble(
    tested = c(TRUE, TRUE, TRUE, TRUE),
    significant = c(TRUE, TRUE, FALSE, FALSE),
    included_background = c(TRUE, TRUE, FALSE, FALSE),
    included_foreground = c(TRUE, FALSE, FALSE, FALSE),
    resolved_name = c("a", "b", NA, "b"),
    exclusion_reason = c(NA, NA, "no portable canonical name", "duplicate_resolved_name")
  )

  txt <- .echogo_representability_statement(mapping)
  expect_match(txt, "2 of 4 tested target entities")
  expect_match(txt, "50.0%")
  expect_match(txt, "matched-background comparative analysis")
})

test_that("missing logo is a graceful presentation-only state", {
  expect_true(is.na(.echogo_logo_source_path(tempfile(fileext = ".png"))))
  expect_true(is.na(.echogo_stage_report_logo(tempfile("echogo-report-"), tempfile(fileext = ".png"))))
  expect_equal(.echogo_report_logo_html(NA_character_), "")
})

test_that("provided logo can be staged without touching scientific evidence", {
  source <- tempfile(fileext = ".png")
  writeBin(as.raw(c(137, 80, 78, 71, 13, 10, 26, 10)), source)
  report_dir <- tempfile("echogo-report-")
  staged <- .echogo_stage_report_logo(report_dir, source)

  expect_true(file.exists(staged))
  expect_equal(basename(staged), "echogo-logo.png")
  html <- .echogo_report_logo_html("echogo-logo.png")
  expect_match(html, "hero-echogo-logo")
  expect_match(html, "EchoGO logo")
})

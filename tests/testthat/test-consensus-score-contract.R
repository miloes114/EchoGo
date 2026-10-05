legacy_evidence_fixture <- function() {
  tibble::tibble(
    term_id = c("GO:0000001", "GO:0000002", "GO:0000003", "GO:0000004"),
    term_name = paste("fixture", 1:4),
    in_goseq = c(TRUE, TRUE, FALSE, TRUE),
    in_hsapiens = c(FALSE, TRUE, TRUE, FALSE),
    in_mmusculus = c(FALSE, FALSE, FALSE, FALSE),
    min_pval_goseq = c(0.01, 0.01, NA, 0.20)
  )
}

test_that("the deprecated score adapter emits scoreless legacy evidence", {
  expect_warning(out <- score_consensus_terms(legacy_evidence_fixture()), "deprecated")
  expect_false(any(c("consensus_score", "consensus_score_all") %in% names(out)))
  expect_identical(
    out$evidence_profile,
    c("TARGET_ONLY", "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT", "NO_PRIMARY_SUPPORT")
  )
})

test_that("legacy GOseq row presence is not target support without qualifying FDR", {
  x <- legacy_evidence_fixture()
  x$in_hsapiens <- FALSE
  x$min_pval_goseq <- c(0.049, 0.051, NA, NA)
  expect_warning(out <- score_consensus_terms(x), "deprecated")
  expect_identical(out$target_goseq_supported, c(TRUE, FALSE, FALSE, FALSE))
})

test_that("legacy adapter never reconstructs a rejected heuristic score", {
  x <- legacy_evidence_fixture()
  x$min_pval_gprof_bg <- c(1e-100, 0.01, 0.02, 0.03)
  x$avg_fold_gprof_bg <- c(100, 4, 2, 6)
  x$depth <- c(2L, 5L, 11L, 20L)
  expect_warning(out <- score_consensus_terms(x), "deprecated")
  expect_false(any(grepl("score|comp_|fold|depth", names(out))))
})

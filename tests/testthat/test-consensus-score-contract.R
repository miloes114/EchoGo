score_fixture <- function(n = 1L) {
  tibble::tibble(
    term_id = sprintf("GO:%07d", seq_len(n)), term_name = paste("fixture", seq_len(n)),
    in_goseq = rep(FALSE, n), in_hsapiens = rep(FALSE, n), in_mmusculus = rep(FALSE, n),
    in_hsapiens_nobg = rep(FALSE, n), in_mmusculus_nobg = rep(FALSE, n),
    min_pval_gprof_bg = rep(NA_real_, n), min_pval_gprof_nobg = rep(NA_real_, n),
    min_pval_goseq = rep(NA_real_, n), avg_fold_gprof_bg = rep(NA_real_, n),
    avg_fold_gprof_nobg = rep(NA_real_, n), fold_enrichment_goseq = rep(NA_real_, n),
    depth = rep(12L, n)
  )
}

test_that("GOseq TRUE FALSE NA are handled element-wise", {
  x <- score_fixture(3); x$in_goseq <- c(TRUE, FALSE, NA)
  out <- score_consensus_terms(x)
  expect_identical(out$comp_goseq, c(1L, 0L, 0L))
  expect_identical(out$sources_count, c(1L, 0L, 0L))
  expect_identical(out$source_origin, c("GOseq", "", ""))
  expect_equal(out$consensus_score, c(1, 0, 0), tolerance = 1e-12)
  expect_equal(out$consensus_score_all, c(1, 0, 0), tolerance = 1e-12)
})

test_that("complete score formulas match a hand calculation", {
  x <- score_fixture(); x$in_goseq <- TRUE; x$in_hsapiens <- TRUE
  x$in_hsapiens_nobg <- TRUE; x$in_mmusculus_nobg <- TRUE
  x$min_pval_gprof_bg <- 1e-3; x$min_pval_gprof_nobg <- 1e-6; x$min_pval_goseq <- 1e-2
  x$avg_fold_gprof_bg <- 3; x$avg_fold_gprof_nobg <- 7; x$fold_enrichment_goseq <- 1
  out <- score_consensus_terms(x)
  conservative <- 1 + 0.8 * 0.5 + 0.4 * 2 + 0.4 * 1 + 0.4 * 0.5
  exploratory <- 1 + 0.7 * 0.5 + 0.3 + 0.35 * 2 + 0.25 * 3 + 0.35 + 0.35
  expect_equal(out$consensus_score, conservative, tolerance = 1e-12)
  expect_equal(out$consensus_score_all, exploratory, tolerance = 1e-12)
})

test_that("score transformations are monotone and capped", {
  weak <- score_fixture(); weak$in_goseq <- TRUE; weak$in_hsapiens <- TRUE
  weak$min_pval_gprof_bg <- 1e-2; weak$avg_fold_gprof_bg <- 2; weak$fold_enrichment_goseq <- 2
  strong <- weak; strong$min_pval_gprof_bg <- 1e-8; strong$avg_fold_gprof_bg <- 20; strong$fold_enrichment_goseq <- 20
  capped <- strong; capped$min_pval_gprof_bg <- 1e-100; capped$avg_fold_gprof_bg <- 200; capped$fold_enrichment_goseq <- 200
  shallow <- strong; shallow$depth <- 3L
  missing_depth <- strong; missing_depth$depth <- NA_integer_
  weak_s <- score_consensus_terms(weak); strong_s <- score_consensus_terms(strong)
  capped_s <- score_consensus_terms(capped); shallow_s <- score_consensus_terms(shallow)
  missing_s <- score_consensus_terms(missing_depth)
  expect_gte(strong_s$consensus_score, weak_s$consensus_score)
  expect_equal(capped_s$comp_p_strict, strong_s$comp_p_strict)
  expect_equal(capped_s$comp_fold_bg, strong_s$comp_fold_bg)
  expect_equal(capped_s$comp_fold_gs, strong_s$comp_fold_gs)
  expect_lt(shallow_s$consensus_score, strong_s$consensus_score)
  expect_equal(missing_s$consensus_score, strong_s$consensus_score)
})

test_that("no-background evidence cannot change the conservative score", {
  base <- score_fixture(); changed <- base
  changed$in_hsapiens_nobg <- TRUE; changed$in_mmusculus_nobg <- TRUE
  changed$min_pval_gprof_nobg <- 1e-8; changed$avg_fold_gprof_nobg <- 8
  base_s <- score_consensus_terms(base); changed_s <- score_consensus_terms(changed)
  expect_equal(changed_s$consensus_score, base_s$consensus_score)
  expect_gt(changed_s$consensus_score_all, base_s$consensus_score_all)
})

test_that("source provenance distinguishes all support combinations", {
  x <- score_fixture(6)
  x$in_goseq <- c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
  x$in_hsapiens <- c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE)
  x$min_pval_gprof_bg <- c(NA, .01, .01, NA, .01, NA)
  x$in_hsapiens_nobg <- c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE)
  x$min_pval_gprof_nobg <- c(NA, NA, NA, .01, .01, NA)
  x$min_pval_goseq <- ifelse(x$in_goseq, .01, NA)
  out <- score_consensus_terms(x)
  expect_identical(out$source_origin, c("GOseq", "g:Profiler_BG", "GOseq+g:Profiler_BG",
    "g:Profiler_noBG", "GOseq+g:Profiler_BG+g:Profiler_noBG", ""))
  expect_identical(out$origin, c("GO terms - GOseq only", "GO terms - g:Profiler only (with BG)",
    "GO terms - Consensus (with BG)", "GO terms - g:Profiler only (no BG)",
    "GO terms - Consensus (with BG)", "Other"))
})

test_that("ambiguous support columns are rejected", {
  x <- score_fixture(); x$in_goseq <- "TRUE"
  expect_error(score_consensus_terms(x), "logical support columns")
  x <- score_fixture(); x$in_hsapiens <- 1L
  expect_error(score_consensus_terms(x), "logical support columns")
})

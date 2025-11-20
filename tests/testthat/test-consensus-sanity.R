test_that("consensus scores match canonical formulas on a real run", {
  skip_if_not_installed("openxlsx")
  run_dir  <- "E:/Gigascience submission/EchoGo_Gigascience_submission/EchoGO/my_first_run/results"
  cons_dir <- file.path(run_dir, "consensus_enrichment")
  xlsx_in  <- file.path(cons_dir, "consensus_enrichment_exploration_clean.xlsx")
  skip_if_not(file.exists(xlsx_in), "No exploration table found")

  cons <- openxlsx::read.xlsx(xlsx_in)

  nz   <- function(x, v=0) { x[is.na(x)] <- v; x }
  zfun <- function(x) log2(1 + pmax(x, 0))
  mlog <- function(p) -log10(pmax(p, .Machine$double.xmin))

  strict_ref <- with(cons, {
    1.5 * as.integer(nz(in_goseq, FALSE)) +
      0.5 * nz(num_species_gprof_bg, 0) +
      zfun(nz(avg_fold_gprof_bg, 0)) +
      zfun(nz(fold_enrichment_goseq, 0)) +
      mlog(nz(min_pval_gprof_bg, 1))
  })

  min_p_any <- pmin(
    nz(cons$min_pval_gprof_bg,   1),
    nz(cons$min_pval_gprof_nobg, 1),
    nz(cons$min_pval_goseq,      1),
    na.rm = TRUE
  )

  all_ref <- with(cons, {
    1.5 * as.integer(nz(in_goseq, FALSE)) +
      0.5 * nz(num_species_gprof_bg,   0) +
      0.5 * nz(num_species_gprof_nobg, 0) +
      zfun(nz(avg_fold_gprof_bg,    0)) +
      zfun(nz(avg_fold_gprof_nobg,  0)) +
      zfun(nz(fold_enrichment_goseq,0)) +
      mlog(min_p_any)
  })

  # Allow tiny FP noise only
  expect_lt(max(abs(nz(cons$consensus_score, 0)     - strict_ref), na.rm = TRUE), 1e-9)
  expect_lt(max(abs(nz(cons$consensus_score_all, 0) - all_ref),    na.rm = TRUE), 1e-9)

  # Also ensure the package scorer agrees on a sample
  if (exists("score_consensus_terms")) {
    rescored <- score_consensus_terms(cons)
    expect_lt(max(abs(nz(rescored$consensus_score,     0) - strict_ref), na.rm = TRUE), 1e-9)
    expect_lt(max(abs(nz(rescored$consensus_score_all, 0) - all_ref),    na.rm = TRUE), 1e-9)
  }
})

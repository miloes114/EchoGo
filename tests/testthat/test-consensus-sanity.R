test_that("deterministic score-contract fixtures replace the obsolete external test", {
  fixture <- testthat::test_path("test-consensus-score-contract.R")
  expect_true(file.exists(fixture))
})

test_that("RRvGO similarity guard handles malformed small-ontology outputs", {
  invalid_objects <- list(
    NULL,
    numeric(),
    NA_real_,
    matrix(NA_real_, nrow = 2, ncol = 2),
    matrix(1, nrow = 1, ncol = 2)
  )

  for (simMatrix in invalid_objects) {
    expect_true(.echogo_invalid_similarity_matrix(simMatrix))
  }

  expect_false(.echogo_invalid_similarity_matrix(diag(2)))
})

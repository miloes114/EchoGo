test_that("exploratory g:Profiler plots retain source rows without inferential filtering", {
  root <- tempfile("echogo-exploratory-plot-")
  custom <- file.path(root, "custom_experimental_background")
  exploratory <- file.path(root, "default_domain_exploratory")
  dir.create(custom, recursive = TRUE)
  dir.create(exploratory, recursive = TRUE)
  fixture <- data.frame(
    term_id = "GO:0000001", term_name = "fixture exploratory term",
    source = "GO:BP", fold_enrichment = 2, p_value = 0.8,
    intersection_size = 2L
  )
  readr::write_csv(fixture, file.path(exploratory, "gprofiler_human_nobg.csv"))
  readr::write_csv(fixture, file.path(custom, "gprofiler_human_with_bg.csv"))

  EchoGO:::gprofiler_make_lollipops(root, c(hsapiens = "human"))

  expect_true(file.exists(file.path(
    exploratory, "gprofiler_human_default_domain_exploratory_GO_BP_lollipop.pdf"
  )))
  expect_true(file.exists(file.path(custom, "gprofiler_human_NO_SIGNIFICANT_RESULTS.txt")))
  expect_false(file.exists(file.path(exploratory, "gprofiler_human_NO_EXPLORATORY_RESULTS.txt")))
})

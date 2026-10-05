test_that("plot-only canonical gProfiler directories cannot shadow cached vectors", {
  td <- tempfile("echogo-gp-cache-")
  dir.create(td, recursive = TRUE)
  canonical <- file.path(td, "custom_experimental_background")
  legacy <- file.path(td, "with_custom_background")
  dir.create(canonical)
  dir.create(legacy)

  file.create(file.path(canonical, "example_plot.pdf"))
  writeLines(c("A", "B"), file.path(legacy, "gprofiler_reference_query.txt"))
  writeLines(c("A", "B", "C"), file.path(legacy, "gprofiler_reference_background.txt"))

  paths <- EchoGO:::.echogo_gprofiler_mode_paths(td, "custom_experimental_background")
  expect_equal(normalizePath(paths[[1]], winslash = "/"), normalizePath(legacy, winslash = "/"))
})

test_that("exact gProfiler result candidates remain canonical-first", {
  td <- tempfile("echogo-gp-results-")
  dir.create(td, recursive = TRUE)
  canonical <- file.path(td, "custom_experimental_background")
  legacy <- file.path(td, "with_custom_background")
  dir.create(canonical)
  dir.create(legacy)

  # Even if the legacy directory owns cached vectors, an actual canonical
  # scientific CSV has precedence for that exact reference species.
  writeLines("A", file.path(legacy, "x_query.txt"))
  writeLines("A", file.path(legacy, "x_background.txt"))
  readr::write_csv(
    tibble::tibble(term_id = "GO:0000001", p_value = 0.01),
    file.path(canonical, "gprofiler_reference_with_bg.csv")
  )
  readr::write_csv(
    tibble::tibble(term_id = "GO:0000002", p_value = 0.02),
    file.path(legacy, "gprofiler_reference_with_bg.csv")
  )

  hit <- EchoGO:::.echogo_find_gprofiler_result_file(
    td, "custom_experimental_background", "reference"
  )
  expect_match(hit, "custom_experimental_background", fixed = TRUE)
})

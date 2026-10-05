# This file exercises the bundled cache only. Keep every case offline so a
# package test never reaches the live g:Profiler species endpoint.
options(
  EchoGO.taxonomy_online = FALSE,
  EchoGO.species_autoupdate = FALSE
)

local_species_cache <- function() {
  cache <- tempfile("echogo-species-cache-")
  dir.create(cache, recursive = TRUE)
  cache
}

test_that("species table has taxonomy + tags offline", {
  cache <- local_species_cache()
  testthat::local_mocked_bindings(.echogo_cache_dir = function() cache, .package = "EchoGO")
  tbl <- echogo_species_table(refresh = FALSE)
  expect_true(all(c("organism","name","ncbi","superkingdom","kingdom","phylum",
                    "class","order","family","genus","tags") %in% names(tbl)))
  expect_true(any(!is.na(tbl$order)))
})

test_that("echogo_resolve returns valid organism IDs only", {
  cache <- local_species_cache()
  testthat::local_mocked_bindings(.echogo_cache_dir = function() cache, .package = "EchoGO")
  ids <- echogo_resolve("tag:AnimalModels OR order:Perciformes")
  tbl <- echogo_species_table(refresh = FALSE)
  expect_type(ids, "character")
  expect_true(length(ids) >= 1)
  expect_true(all(!is.na(ids) & nzchar(ids)))
  expect_true(all(ids %in% tbl$organism))
})

test_that("tag-only selection returns the tagged subset", {
  cache <- local_species_cache()
  testthat::local_mocked_bindings(.echogo_cache_dir = function() cache, .package = "EchoGO")
  tbl <- echogo_species_table(refresh = FALSE)
  expected <- tbl$organism[
    vapply(tbl$tags, function(x) "AnimalModels" %in% x, logical(1))
  ]

  selected <- echogo_select_species(tags = "AnimalModels", refresh = FALSE)

  expect_setequal(selected, expected)
  expect_lt(length(selected), nrow(tbl))
})

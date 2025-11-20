test_that("species table has taxonomy + tags offline", {
  old <- getOption("EchoGO.taxonomy_online")
  options(EchoGO.taxonomy_online = FALSE)
  on.exit(options(EchoGO.taxonomy_online = old), add = TRUE)

  tbl <- echogo_species_table(refresh = FALSE)
  expect_true(all(c("organism","name","ncbi","superkingdom","kingdom","phylum",
                    "class","order","family","genus","tags") %in% names(tbl)))
  expect_true(any(!is.na(tbl$order)))
})

test_that("echogo_resolve returns valid organism IDs only", {
  ids <- echogo_resolve("tag:AnimalModels OR order:Perciformes")
  tbl <- echogo_species_table(refresh = FALSE)
  expect_type(ids, "character")
  expect_true(length(ids) >= 1)
  expect_true(all(!is.na(ids) & nzchar(ids)))
  expect_true(all(ids %in% tbl$organism))
})

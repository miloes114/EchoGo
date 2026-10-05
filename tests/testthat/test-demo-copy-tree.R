test_that("packaged demo copy handles mixed top-level files and nested cache directories", {
  root <- tempfile("echogo-demo-copy-")
  src <- file.path(root, "source")
  dst <- file.path(root, "destination")
  dir.create(src, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  writeLines("root input", file.path(src, "DE_results_demo.tsv"))
  dir.create(file.path(src, "cache", "bg", "drerio"), recursive = TRUE)
  writeLines(
    "cached response",
    file.path(src, "cache", "bg", "drerio", "GO_BP.csv")
  )
  dir.create(file.path(src, "empty_nested_directory"), recursive = TRUE)

  copied <- .echogo_copy_tree(src, dst)

  expect_true(file.exists(file.path(dst, "DE_results_demo.tsv")))
  expect_true(file.exists(file.path(
    dst,
    "cache",
    "bg",
    "drerio",
    "GO_BP.csv"
  )))
  expect_true(dir.exists(file.path(dst, "empty_nested_directory")))
  expect_equal(readLines(file.path(dst, "DE_results_demo.tsv")), "root input")
  expect_equal(
    readLines(file.path(
      dst,
      "cache",
      "bg",
      "drerio",
      "GO_BP.csv"
    )),
    "cached response"
  )
  expect_length(copied, 2L)
})

test_that("packaged demo copy has matched source and destination vectors on overwrite", {
  root <- tempfile("echogo-demo-copy-overwrite-")
  src <- file.path(root, "source")
  dst <- file.path(root, "destination")
  dir.create(file.path(src, "cache", "nested"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  writeLines("v1", file.path(src, "input.tsv"))
  writeLines("cache-v1", file.path(src, "cache", "nested", "response.tsv"))
  .echogo_copy_tree(src, dst)

  writeLines("v2", file.path(src, "input.tsv"))
  writeLines("cache-v2", file.path(src, "cache", "nested", "response.tsv"))
  expect_silent(.echogo_copy_tree(src, dst, overwrite = TRUE))

  expect_equal(readLines(file.path(dst, "input.tsv")), "v2")
  expect_equal(readLines(file.path(dst, "cache", "nested", "response.tsv")), "cache-v2")
})

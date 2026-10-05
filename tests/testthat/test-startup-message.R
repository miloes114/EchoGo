test_that("startup remains readable without ANSI colours", {
  withr::local_options(cli.num_colors = 1L)
  startup <- testthat::capture_messages(
    EchoGO:::.onAttach("", "EchoGO")
  )
  text <- paste(startup, collapse = "\n")

  expect_match(text, paste0("EchoGO v", utils::packageVersion("EchoGO")), fixed = TRUE)
  expect_false(grepl("\033", text, fixed = TRUE))
  expect_true(all(nchar(strsplit(text, "\n", fixed = TRUE)[[1]]) <= 80L))
})

test_that("colour styling preserves content and startup suppression", {
  skip_if_not_installed("cli")
  withr::local_options(cli.num_colors = 1L)
  plain <- testthat::capture_messages(EchoGO:::.onAttach("", "EchoGO"))
  options(cli.num_colors = 256L)
  colour <- testthat::capture_messages(EchoGO:::.onAttach("", "EchoGO"))
  expect_identical(cli::ansi_strip(colour), plain)
  expect_true(any(grepl("\033", colour, fixed = TRUE)))
  expect_length(testthat::capture_messages(
    suppressPackageStartupMessages(EchoGO:::.onAttach("", "EchoGO"))
  ), 0L)
})

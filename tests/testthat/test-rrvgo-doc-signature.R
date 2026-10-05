test_that("RRvGO compatibility wrapper and help usage remain variadic", {
  expect_identical(names(formals(EchoGO::run_rrvgo_consensus_analysis)), "...")

  rd_file <- file.path(
    testthat::test_path("..", "..", "man"),
    "run_rrvgo_consensus_analysis.Rd"
  )
  if (file.exists(rd_file)) {
    rd <- tools::parse_Rd(rd_file)
  } else {
    topic <- utils::help("run_rrvgo_consensus_analysis", package = "EchoGO")
    expect_true(length(topic) == 1L)
    rd <- utils:::.getHelpFile(topic)
  }
  help_text <- paste(capture.output(tools::Rd2txt(rd)), collapse = "\n")
  expect_match(help_text, "run_rrvgo_consensus_analysis\\(\\.\\.\\.")
  expect_match(help_text, "Arguments forwarded")
})

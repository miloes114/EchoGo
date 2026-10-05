test_that("own-data quickstart help declares the target-anchored configuration", {
  lines <- capture.output(echogo_help())
  begin <- match("Use your own data:", lines)
  end <- grep("^Reference-based RNA-seq", lines)[1]
  expect_false(is.na(begin))
  expect_false(is.na(end))
  examples <- parse(text = lines[seq.int(begin + 1L, end - 1L)])
  calls <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("echogo_run")), examples)
  expect_length(calls, 1L)
  call <- calls[[1]]
  expect_true(all(c("input_dir", "outdir", "species", "target_context",
                    "run_exploratory_default_domain", "run_rrvgo",
                    "semantic_reference_orgdb", "semantic_reference_role") %in% names(call)))
  expect_identical(call$run_exploratory_default_domain, FALSE)
  expect_identical(call$semantic_reference_role, "target_reference")
  # The wrapper's ... must forward a valid current run_full_echogo signature.
  call[[1]] <- as.name("run_full_echogo")
  expect_silent(match.call(run_full_echogo, call))
  expect_true(any(grepl("REQUIRED precomputed target GOseq", lines, fixed = TRUE)))
  expect_true(any(grepl("Target GOseq stays required", lines, fixed = TRUE)))
})

test_that("installed workflow help requires GOseq independently of target context", {
  source_rd <- testthat::test_path("..", "..", "man", "run_full_echogo.Rd")
  rd <- if (file.exists(source_rd)) tools::parse_Rd(source_rd) else
    tools::Rd_db("EchoGO")[["run_full_echogo.Rd"]]
  text <- gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = " "))
  expect_match(text, "Target GOseq provides the primary target evidence", fixed = TRUE)
  expect_match(text, "itself is required for the standard workflow", fixed = TRUE)
  expect_match(text, "target GOseq remains the experimental anchor", fixed = TRUE)
  expect_false(grepl("\\(optional\\) a GOseq|GOseq enriched categories; used when present", text))
})

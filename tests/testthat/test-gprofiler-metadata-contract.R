test_that("gProfiler runs preserve exact vectors and complete mode metadata", {
  fake_gost <- function(query, organism, custom_bg, ...) {
    domain <- if (organism == "hsapiens") 8L else 7L
    query_n <- if (organism == "hsapiens") 2L else 1L
    list(
      result = data.frame(
        term_id = "GO:0000001",
        term_name = "fixture term",
        source = "GO:BP",
        p_value = 0.01,
        intersection = "GENEA,GENEC",
        intersection_size = 2L,
        query_size = query_n,
        term_size = 3L,
        effective_domain_size = domain,
        stringsAsFactors = FALSE
      ),
      meta = list(
        source = "mocked",
        database_version = paste0("mock-", organism),
        effective_query_size = query_n,
        effective_domain_size = domain
      )
    )
  }
  testthat::local_mocked_bindings(.echogo_gost = fake_gost)
  out <- tempfile("echogo-gprofiler-")
  options(EchoGO.force_rerun_gprofiler = TRUE)
  withr::defer(options(EchoGO.force_rerun_gprofiler = NULL))

  result <- run_gprofiler_cross_species(
    de_genes = c("GENEA", "GENEC"),
    bg_genes = c("GENEA", "GENEB", "GENEC", "GENED"),
    species = c(hsapiens = "human", mmusculus = "mouse"),
    outdir = out,
    do_no_bg = TRUE,
    sleep_sec = 0,
    verbose = FALSE,
    significance_rule = list(type = "threshold", padj_threshold = 0.05),
    resolver_definition = list(
      route = "shared_portable_canonical_organism_context",
      explicit_species_specific_ortholog_mapping = FALSE
    )
  )

  manifest_path <- file.path(out, "run_manifest.json")
  expect_true(file.exists(manifest_path))
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
  expect_length(manifest$runs, 4L)
  expect_identical(manifest$vector_contract, "shared_portable_canonical_organism_context")
  expect_false(manifest$explicit_species_specific_ortholog_mapping)

  for (entry in manifest$runs) {
    expect_true(file.exists(file.path(out, entry$query_file)))
    expect_true(file.exists(file.path(out, entry$metadata_file)))
    expect_identical(readLines(file.path(out, entry$query_file)), c("GENEA", "GENEC"))
    expect_identical(unname(tools::md5sum(file.path(out, entry$query_file))), entry$query_hash)
    if (entry$background_mode == "custom_background") {
      expect_true(file.exists(file.path(out, entry$background_file)))
      expect_identical(readLines(file.path(out, entry$background_file)), c("GENEA", "GENEB", "GENEC", "GENED"))
      expect_identical(unname(tools::md5sum(file.path(out, entry$background_file))), entry$background_hash)
    } else {
      expect_length(entry$background_file, 0L)
    }
  }

  human <- Filter(function(x) x$organism_code == "hsapiens" && x$background_mode == "custom_background", manifest$runs)[[1]]
  mouse <- Filter(function(x) x$organism_code == "mmusculus" && x$background_mode == "custom_background", manifest$runs)[[1]]
  expect_identical(human$submitted_foreground_count, mouse$submitted_foreground_count)
  expect_false(identical(human$effective_query_size, mouse$effective_query_size))
  expect_true(all(file.exists(result$paths$written)))
})

test_that("gProfiler custom mode rejects missing or degenerate backgrounds", {
  expect_error(
    run_gprofiler_cross_species(c("A"), bg_genes = NULL, species = "hsapiens", verbose = FALSE),
    "custom background"
  )
  expect_error(
    run_gprofiler_cross_species(c("A", "B"), bg_genes = c("A", "B"), species = "hsapiens", verbose = FALSE),
    "same vector"
  )
  expect_error(
    run_gprofiler_cross_species(c("A", "C"), bg_genes = c("A", "B"), species = "hsapiens", verbose = FALSE),
    "subset"
  )
})

test_that("gProfiler metadata distinguishes empty responses from request failures", {
  run_with_mock <- function(mock) {
    testthat::with_mocked_bindings(
      {
        out <- tempfile("echogo-gprofiler-status-")
        withr::local_options(list(EchoGO.force_rerun_gprofiler = TRUE))
        run_gprofiler_cross_species(
          de_genes = "A",
          bg_genes = c("A", "B"),
          species = "hsapiens",
          outdir = out,
          do_no_bg = FALSE,
          sleep_sec = 0,
          verbose = FALSE
        )
        jsonlite::read_json(file.path(out, "run_manifest.json"), simplifyVector = FALSE)
      },
      .echogo_gost = mock,
      .package = "EchoGO"
    )
  }

  empty <- run_with_mock(function(...) NULL)
  expect_identical(empty$runs[[1]]$status, "no_results")

  expect_warning(
    failed <- run_with_mock(function(...) stop("simulated transport failure")),
    "simulated transport failure"
  )
  expect_identical(failed$runs[[1]]$status, "request_failed")
})

demo_source_dir <- function() {
  system.file("extdata", "echogo_demo", package = "EchoGO")
}

test_that("packaged demo contains the canonical denominator-aware TSV only", {
  demo_dir <- demo_source_dir()
  tsv <- file.path(demo_dir, "GOseq_enrichment_demo.tsv")
  csv <- file.path(demo_dir, "GOseq_enrichment_demo.csv")

  expect_true(file.exists(tsv))
  expect_false(file.exists(csv))

  goseq <- utils::read.delim(
    tsv,
    sep = "\t",
    quote = "\"",
    check.names = FALSE,
    nrows = 5
  )
  expect_identical(
    names(goseq),
    c(
      "category",
      "term",
      "ontology",
      "numDEInCat",
      "numInCat",
      "total_significant_genes",
      "total_tested_genes",
      "over_represented_FDR",
      "gene_ids"
    )
  )
})

test_that("quoted CSV gene lists remain in one gene_ids column", {
  csv <- tempfile(fileext = ".csv")
  writeLines(
    c(
      "category,term,ontology,numDEInCat,numInCat,over_represented_FDR,gene_ids",
      "GO:0000001,example,BP,2,5,0.01,\"gene_a,gene_b,gene_c\""
    ),
    csv
  )

  parsed <- EchoGO:::.echogo_read_delim_robust(
    csv,
    expected = c(
      "category", "term", "ontology", "numdeincat", "numincat",
      "over_represented_fdr", "gene_ids"
    )
  )

  expect_identical(ncol(parsed), 7L)
  expect_identical(parsed$gene_ids, "gene_a,gene_b,gene_c")
  expect_silent(EchoGO:::.echogo_check_goseq_parse(parsed))
})

test_that("resolver honors pattern rank and prefers TSV over CSV", {
  input_dir <- tempfile("echogo-resolver-")
  dir.create(input_dir)
  tsv <- file.path(input_dir, "sample.GOseq.enriched.tsv")
  csv <- file.path(input_dir, "sample.GOseq.enriched.csv")
  file.create(tsv, csv)

  found <- EchoGO:::.echogo_find_one(
    file.path(input_dir, c("*.GOseq.enriched.tsv", "*.GOseq.enriched.csv")),
    label = "GOseq"
  )

  expect_identical(normalizePath(found), normalizePath(tsv))
})

test_that("resolver rejects ambiguity within the highest matching rank", {
  input_dir <- tempfile("echogo-ambiguous-")
  dir.create(input_dir)
  file.create(
    file.path(input_dir, "a.GOseq.enriched.tsv"),
    file.path(input_dir, "b.GOseq.enriched.tsv")
  )

  expect_error(
    EchoGO:::.echogo_find_one(
      file.path(input_dir, "*.GOseq.enriched.tsv"),
      label = "GOseq"
    ),
    "Multiple candidate GOseq files matched"
  )
})

test_that("missing RRvGO reports the dependency instead of matrix sparsity", {
  expect_error(
    EchoGO:::.echogo_require_rrvgo(FALSE),
    "package 'rrvgo'.*not installed",
    fixed = FALSE
  )
})

test_that("OrgDb installer defaults to the active first library", {
  default_lib <- formals(echogo_install_orgdb)$lib
  expect_identical(deparse(default_lib), ".libPaths()[1]")
  expect_identical(eval(default_lib), .libPaths()[1])
})

test_that("quickstart passes explicit demo files and cleans stale results", {
  captured <- NULL

  testthat::local_mocked_bindings(
    echogo_require_orgdb = function(...) invisible(TRUE),
    run_full_echogo = function(...) {
      captured <<- list(...)
      list(files = list())
    },
    .validate_quickstart_result = function(...) invisible(TRUE),
    .package = "EchoGO"
  )

  outdir <- tempfile("echogo-quickstart-")
  copied <- echogo_quickstart(run_demo = FALSE, outdir = outdir)
  stale <- file.path(copied, "results", "stale.txt")
  dir.create(dirname(stale), recursive = TRUE)
  file.create(stale)

  echogo_quickstart(
    run_demo = TRUE,
    outdir = outdir,
    make_report = FALSE,
    clean = TRUE
  )

  expect_false(file.exists(stale))
  expect_identical(basename(captured$goseq_file), "GOseq_enrichment_demo.tsv")
  expect_identical(basename(captured$trinotate_file), "Trinotate_demo.tsv")
  expect_identical(basename(captured$de_file), "DE_results_demo.tsv")
  expect_identical(basename(captured$count_matrix_file), "counts_demo.tsv")
  expect_identical(captured$target_context, "drerio")
  expect_identical(captured$semantic_reference_orgdb, "org.Dr.eg.db")
  expect_identical(captured$semantic_reference_role, "target_reference")
  expect_false(captured$run_exploratory_default_domain)
})

test_that("packaged demo snapshot identifies the v0.1.4 offline product", {
  frozen <- echogo_demo_results_path()
  readme <- file.path(frozen, "README.txt")
  expect_true(file.exists(readme))
  text <- paste(readLines(readme, warn = FALSE), collapse = "\n")
  expect_match(text, "v0.1.4")
  expect_false(grepl("consensus_score|true_consensus|no-background", text, ignore.case = TRUE))
})

test_that("offline basic quickstart produces consensus without an optional OrgDb", {
  orgdb_checked <- FALSE
  cache <- tempfile("echogo-species-cache-")
  dir.create(cache)
  testthat::local_mocked_bindings(
    echogo_require_orgdb = function(...) {
      orgdb_checked <<- TRUE
      stop("basic quickstart must not enter semantic OrgDb preparation")
    },
    .echogo_cache_dir = function() cache,
    .package = "EchoGO"
  )
  res <- echogo_quickstart(
    run_demo = TRUE,
    outdir = tempfile("echogo-live-"),
    make_report = FALSE,
    clean = TRUE,
    live_gprofiler = FALSE
  )

  expect_true(file.exists(res$files$exact_term_evidence))
  expect_gt(nrow(readr::read_csv(res$files$exact_term_evidence, show_col_types = FALSE)), 0L)
  expect_false(orgdb_checked)
})

test_that("offline quickstart never calls the species network fetcher with an empty cache", {
  cache <- tempfile("echogo-species-empty-")
  dir.create(cache, recursive = TRUE)
  network_called <- FALSE
  testthat::local_mocked_bindings(
    .echogo_cache_dir = function() cache,
    .echogo_fetch_species_online = function(...) {
      network_called <<- TRUE
      stop("NETWORK SHOULD NOT BE CALLED")
    },
    .package = "EchoGO"
  )
  old <- options(EchoGO.species_autoupdate = TRUE, EchoGO.taxonomy_online = TRUE)
  on.exit(options(old), add = TRUE)
  res <- echogo_quickstart(
    run_demo = TRUE, outdir = tempfile("echogo-empty-cache-"),
    make_report = FALSE, full = FALSE, live_gprofiler = FALSE
  )
  expect_true(file.exists(res$files$consensus_xlsx))
  expect_false(network_called)
})

test_that("offline quickstart never refreshes a stale species cache", {
  cache <- tempfile("echogo-species-stale-")
  dir.create(cache, recursive = TRUE)
  fallback <- EchoGO:::.echogo_species_fallback()
  readr::write_csv(fallback, file.path(cache, "gprofiler_species.csv"))
  jsonlite::write_json(
    list(source = "api", timestamp = 0, ok = TRUE, version = "0.1.4"),
    file.path(cache, "gprofiler_species_meta.json"), auto_unbox = TRUE
  )
  network_called <- FALSE
  testthat::local_mocked_bindings(
    .echogo_cache_dir = function() cache,
    .echogo_fetch_species_online = function(...) {
      network_called <<- TRUE
      stop("NETWORK SHOULD NOT BE CALLED")
    },
    .package = "EchoGO"
  )
  old <- options(EchoGO.species_autoupdate = TRUE, EchoGO.taxonomy_online = TRUE)
  on.exit(options(old), add = TRUE)
  res <- echogo_quickstart(
    run_demo = TRUE, outdir = tempfile("echogo-stale-cache-"),
    make_report = FALSE, full = FALSE, live_gprofiler = FALSE
  )
  expect_true(file.exists(res$files$consensus_xlsx))
  expect_false(network_called)
})

test_that("offline quickstart restores species and taxonomy options on success and error", {
  cache <- tempfile("echogo-species-restore-")
  dir.create(cache, recursive = TRUE)
  testthat::local_mocked_bindings(
    .echogo_cache_dir = function() cache,
    .echogo_fetch_species_online = function(...) stop("NETWORK SHOULD NOT BE CALLED"),
    .package = "EchoGO"
  )
  old <- options(EchoGO.species_autoupdate = "before-auto", EchoGO.taxonomy_online = "before-tax")
  on.exit(options(old), add = TRUE)
  expect_no_error(echogo_quickstart(
    run_demo = TRUE, outdir = tempfile("echogo-restore-ok-"),
    make_report = FALSE, full = FALSE, live_gprofiler = FALSE
  ))
  expect_identical(getOption("EchoGO.species_autoupdate"), "before-auto")
  expect_identical(getOption("EchoGO.taxonomy_online"), "before-tax")

  testthat::local_mocked_bindings(
    run_full_echogo = function(...) stop("intentional downstream failure"),
    .validate_quickstart_result = function(...) invisible(TRUE),
    .package = "EchoGO"
  )
  expect_error(echogo_quickstart(
    run_demo = TRUE, outdir = tempfile("echogo-restore-error-"),
    make_report = FALSE, full = FALSE, live_gprofiler = FALSE
  ), "intentional downstream failure")
  expect_identical(getOption("EchoGO.species_autoupdate"), "before-auto")
  expect_identical(getOption("EchoGO.taxonomy_online"), "before-tax")
})

test_that("live quickstart mode does not apply the offline species guard", {
  captured <- NULL
  testthat::local_mocked_bindings(
    run_full_echogo = function(...) {
      captured <<- list(
        species_autoupdate = getOption("EchoGO.species_autoupdate"),
        taxonomy_online = getOption("EchoGO.taxonomy_online")
      )
      list(files = list())
    },
    .validate_quickstart_result = function(...) invisible(TRUE),
    .package = "EchoGO"
  )
  old <- options(EchoGO.species_autoupdate = TRUE, EchoGO.taxonomy_online = TRUE)
  on.exit(options(old), add = TRUE)
  echogo_quickstart(
    run_demo = TRUE, outdir = tempfile("echogo-live-mode-"),
    make_report = FALSE, full = FALSE, live_gprofiler = TRUE
  )
  expect_true(captured$species_autoupdate)
  expect_true(captured$taxonomy_online)
})

test_that("full quickstart passes an explicit demo semantic reference without installing", {
  run_args <- NULL
  testthat::local_mocked_bindings(
    echogo_require_orgdb = function(...) stop("quickstart must not install dependencies"),
    run_full_echogo = function(...) {
      run_args <<- list(...)
      list(files = list())
    },
    .validate_quickstart_result = function(...) invisible(TRUE),
    .package = "EchoGO"
  )

  echogo_quickstart(
    run_demo = TRUE,
    outdir = tempfile("echogo-full-"),
    make_report = FALSE,
    clean = TRUE,
    full = TRUE,
    live_gprofiler = FALSE
  )

  expect_true(run_args$run_rrvgo)
  expect_true(run_args$run_evaluation)
  expect_true(run_args$run_exploratory_default_domain)
  expect_identical(run_args$semantic_reference_orgdb, "org.Dr.eg.db")
  expect_identical(run_args$semantic_reference_role, "target_reference")
  expect_identical(run_args$target_context, "drerio")
})

test_that("GO depth preserves mixed-ID row alignment", {
  input <- tibble::tibble(
    term_id = c(
      "GO:0008150", "prefix GO:0003674 suffix", "KEGG:04010",
      "not_a_go_identifier", NA_character_, ""
    )
  )

  expect_warning(
    result <- suppressMessages(EchoGO:::.echogo_add_go_depth(input)),
    NA
  )
  expect_equal(nrow(result), nrow(input))
  expect_length(result$depth, nrow(input))
  expect_type(result$depth, "integer")
  expect_true(all(is.na(result$depth[c(3, 4, 5, 6)])))
  if (requireNamespace("GO.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    expect_false(anyNA(result$depth[c(1, 2)]))
  }
})

test_that("echogo_run accepts overrides and uses the configured report title", {
  captured <- NULL

  testthat::local_mocked_bindings(
    echogo_require_orgdb = function(...) invisible(TRUE),
    run_full_echogo = function(...) {
      captured <<- list(...)
      invisible(captured)
    },
    .package = "EchoGO"
  )

  input_dir <- tempfile("echogo-run-")
  dir.create(input_dir)
  writeLines(
    c(
      "species: [hsapiens, mmusculus]",
      "target_context: NO_TARGET",
      "run_rrvgo: false",
      "report_title: Configured report"
    ),
    file.path(input_dir, "config.yml")
  )

  echogo_run(
    input_dir = input_dir,
    outdir = file.path(input_dir, "results"),
    make_report = FALSE,
    run_rrvgo = FALSE
  )

  expect_false(captured$make_report)
  expect_false(captured$run_rrvgo)
  expect_identical(captured$report_title, "Configured report")
  expect_identical(unname(captured$species), c("hsapiens", "mmusculus"))
})

test_that("reference resolver rejects ambiguous contrast folders", {
  input_dir <- tempfile("echogo-reference-ambiguous-")
  dir.create(input_dir)
  file.create(
    file.path(input_dir, "a.GOseq.enriched.tsv"),
    file.path(input_dir, "b.GOseq.enriched.tsv"),
    file.path(input_dir, "Trinotate_for_EchoGO.tsv")
  )

  expect_error(
    echogo_resolve_reference_inputs(input_dir),
    "Multiple GOseq enriched files"
  )
})

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
})

test_that("frozen demo outputs remain non-empty and include an HTML report", {
  frozen <- echogo_demo_results_path()
  consensus_file <- file.path(
    frozen,
    "consensus",
    "consensus_enrichment_results_with_and_without_bg.xlsx"
  )
  reports <- list.files(
    file.path(frozen, "report"),
    pattern = "\\.html$",
    full.names = TRUE
  )

  expect_true(file.exists(consensus_file))
  expect_gt(nrow(openxlsx::read.xlsx(consensus_file)), 0L)
  expect_true(length(reports) > 0L)
  expect_true(all(file.info(reports)$size > 0))

  report_html <- paste(
    readLines(reports[[1]], warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )
  asset_tags <- regmatches(
    report_html,
    gregexpr(
      "(?:src|href|data)=[\"']assets/[^\"']+[\"']",
      report_html,
      perl = TRUE,
      ignore.case = TRUE
    )
  )[[1]]
  asset_refs <- sub(
    "^(?:src|href|data)=[\"']",
    "",
    asset_tags,
    perl = TRUE,
    ignore.case = TRUE
  )
  asset_refs <- unique(utils::URLdecode(sub("[\"']$", "", asset_refs)))

  expect_gt(length(asset_refs), 0L)
  asset_paths <- file.path(dirname(reports[[1]]), asset_refs)
  missing_assets <- asset_refs[!file.exists(asset_paths)]
  expect_length(missing_assets, 0L)
})

test_that("offline quickstart produces consensus without network access", {
  res <- echogo_quickstart(
    run_demo = TRUE,
    outdir = tempfile("echogo-live-"),
    make_report = FALSE,
    clean = TRUE,
    live_gprofiler = FALSE
  )

  expect_true(file.exists(res$files$consensus_xlsx))
  expect_gt(nrow(openxlsx::read.xlsx(res$files$consensus_xlsx)), 0L)
})

test_that("GO depth preserves rows with mixed identifier types", {
  input <- tibble::tibble(
    term_id = c("GO:0008150", "KEGG:04010", NA_character_, "prefix GO:0003674 suffix")
  )

  expect_warning(
    result <- suppressMessages(EchoGO:::.echogo_add_go_depth(input)),
    NA
  )
  expect_equal(nrow(result), nrow(input))
  expect_length(result$depth, nrow(input))
  expect_true(is.na(result$depth[[2]]))
  expect_true(is.na(result$depth[[3]]))
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
      "orgdb: org.Mm.eg.db",
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

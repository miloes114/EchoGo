test_that("packaged demo is a valid matched-universe scientific fixture", {
  demo <- system.file("extdata", "echogo_demo", package = "EchoGO")
  expect_true(nzchar(demo))
  de <- readr::read_tsv(file.path(demo, "DE_results_demo.tsv"), show_col_types = FALSE)
  counts <- readr::read_tsv(file.path(demo, "counts_demo.tsv"), show_col_types = FALSE)
  annotation <- readr::read_tsv(file.path(demo, "Trinotate_demo.tsv"), show_col_types = FALSE)
  goseq <- readr::read_tsv(file.path(demo, "GOseq_enrichment_demo.tsv"), show_col_types = FALSE)

  expect_equal(nrow(de), 120L)
  expect_equal(nrow(counts), 120L)
  expect_equal(nrow(annotation), 118L)
  expect_false(anyNA(de$gene_id))
  expect_equal(anyDuplicated(de$gene_id), 0L)
  expect_false(anyNA(counts$gene_id))
  expect_equal(anyDuplicated(counts$gene_id), 0L)
  expect_true(any(de$significant))
  expect_true(any(!de$significant))
  expect_setequal(de$gene_id, counts$gene_id)

  sets <- prepare_gprofiler_gene_sets(
    de,
    annotation = annotation,
    count_matrix = counts,
    use_trinotate_universe = TRUE
  )
  expect_length(sets$foreground_original, 24L)
  expect_length(sets$background_original, 120L)
  expect_length(sets$foreground_unmapped, 0L)
  expect_length(sets$background_unmapped, 2L)
  expect_equal(sets$duplicates_collapsed$foreground, 1L)
  expect_equal(sets$duplicates_collapsed$background, 1L)
  expect_length(setdiff(sets$foreground_resolved, sets$background_resolved), 0L)
  expect_lt(length(sets$foreground_resolved), length(sets$background_resolved))

  folded <- .echogo_compute_goseq_fold(goseq, 24L, 120L)
  expect_equal(folded$foldEnrichment[[1]], 5, tolerance = 1e-12)
  expect_true(all(goseq$total_significant_genes == 24L))
  expect_true(all(goseq$total_tested_genes == 120L))
})

test_that("packaged demo cache is explicit, complete, and offline", {
  demo <- system.file("extdata", "echogo_demo", package = "EchoGO")
  expect_true(nzchar(demo))
  cache <- file.path(demo, "cached_gprofiler_v0.1.3")
  manifest <- jsonlite::read_json(file.path(cache, "run_manifest.json"), simplifyVector = FALSE)
  expect_identical(manifest$execution, "cached_demo_fixture")
  expect_identical(manifest$vector_contract, "shared_portable_canonical_organism_context")
  expect_false(manifest$explicit_species_specific_ortholog_mapping)
  expect_length(manifest$runs, 6L)
  for (entry in manifest$runs) {
    expect_true(file.exists(file.path(cache, entry$result_file)))
    expect_true(file.exists(file.path(cache, entry$metadata_file)))
    expect_true(file.exists(file.path(cache, entry$query_file)))
    if (entry$background_mode == "custom_background") {
      expect_true(file.exists(file.path(cache, entry$background_file)))
    }
  }
})

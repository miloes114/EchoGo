test_that("GOseq fold uses significant and tested denominators", {
  terms <- data.frame(
    numDEInCat = c(24, 2),
    numInCat = c(24, 20)
  )
  out <- .echogo_compute_goseq_fold(
    terms,
    total_significant_genes = 24,
    total_tested_genes = 120
  )
  expect_equal(out$foldEnrichment, c(5, 0.5), tolerance = 1e-12)
  expect_identical(out$total_significant_genes, c(24L, 24L))
  expect_identical(out$total_tested_genes, c(120L, 120L))
})

test_that("GOseq denominator validation rejects incompatible universes", {
  expect_error(
    .echogo_compute_goseq_fold(data.frame(numDEInCat = 25, numInCat = 25), 24, 120),
    "numDEInCat exceeds"
  )
  expect_error(
    .echogo_compute_goseq_fold(data.frame(numDEInCat = 1, numInCat = 121), 24, 120),
    "numInCat exceeds"
  )
  expect_error(
    .echogo_compute_goseq_fold(data.frame(numDEInCat = 2, numInCat = 1), 24, 120),
    "numDEInCat exceeds numInCat"
  )
  expect_error(
    .echogo_compute_goseq_fold(data.frame(numDEInCat = 1, numInCat = 1), 0, 120),
    "positive"
  )
  expect_error(
    .echogo_compute_goseq_fold(data.frame(numDEInCat = 1, numInCat = 1), 24, NA),
    "positive"
  )
})

test_that("declared significant-only DE tables use an explicit larger universe", {
  de <- data.frame(gene_id = paste0("g", 1:24))
  counts <- data.frame(gene_id = paste0("g", 1:120), sample = 1L)
  annotation <- data.frame(
    transcript_id = paste0("g", 1:120),
    EggNM.Preferred_name = paste0("G", 1:120)
  )
  sets <- prepare_gprofiler_gene_sets(
    de,
    annotation = annotation,
    count_matrix = counts,
    de_table_significant_only = TRUE,
    use_trinotate_universe = TRUE
  )
  expect_length(sets$foreground_original, 24L)
  expect_length(sets$background_original, 120L)
  out <- .echogo_compute_goseq_fold(data.frame(numDEInCat = 24, numInCat = 24), 24, 120)
  expect_equal(out$foldEnrichment, 5)
})

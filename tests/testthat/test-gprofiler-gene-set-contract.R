gene_set_fixture <- function() {
  list(
    de = data.frame(
      gene_id = c("tx1", "tx2", "tx3", "tx4", "tx5", "tx6"),
      log2FoldChange = c(2, 0.2, -2, 3, 0.1, 0),
      padj = c(0.001, 0.001, 0.02, NA, 0.8, 0.9),
      stringsAsFactors = FALSE
    ),
    annotation = data.frame(
      transcript_id = c("tx1", "tx2", "tx3", "tx4", "tx5"),
      sprot_Top_BLASTX_hit = c("GENEA^Metazoa", "GENEB^Metazoa", NA, "GENED^Metazoa", "GENEB^Metazoa"),
      EggNM.Preferred_name = c(NA, NA, "GENEC", NA, NA),
      stringsAsFactors = FALSE
    )
  )
}

test_that("foreground contains only genes passing the explicit DE rule", {
  f <- gene_set_fixture()
  sets <- prepare_gprofiler_gene_sets(f$de, annotation = f$annotation)
  expect_identical(sets$foreground_original, c("tx1", "tx3"))
  expect_identical(sets$foreground_resolved, c("GENEA", "GENEC"))
  expect_identical(sets$significance_rule$padj_threshold, 0.05)
  expect_identical(sets$significance_rule$log2fc_threshold, 1)
})

test_that("tested background records duplicate and unmapped genes", {
  f <- gene_set_fixture()
  sets <- prepare_gprofiler_gene_sets(
    f$de,
    tested_gene_ids = c("tx1", "tx2", "tx3", "tx5", "tx6"),
    annotation = f$annotation,
    use_trinotate_universe = TRUE
  )
  expect_identical(sets$background_resolved, c("GENEA", "GENEB", "GENEC"))
  expect_identical(sets$background_unmapped, "tx6")
  expect_equal(sets$duplicates_collapsed$background, 1L)
  expect_true(all(c(
    "original_gene_id", "portable_name", "portable_name_source", "native_symbol",
    "swissprot_accession", "swissprot_gene_symbol", "eggnog_preferred_name",
    "annotation_taxonomy", "included_in_gprofiler", "exclusion_reason",
    "original_id", "tested", "significant", "resolved_name", "resolution_source",
    "included_foreground", "included_background", "duplicate_group"
  ) %in% names(sets$mapping_table)))
})

test_that("foreground and background use the same resolver", {
  f <- gene_set_fixture()
  sets <- prepare_gprofiler_gene_sets(f$de, annotation = f$annotation)
  row <- sets$mapping_table[sets$mapping_table$original_id == "tx1", ]
  expect_true(row$included_foreground)
  expect_true(row$included_background)
  expect_identical(row$resolved_name, "GENEA")
})

test_that("foreground outside the explicit tested universe is rejected", {
  f <- gene_set_fixture()
  expect_error(
    prepare_gprofiler_gene_sets(
      f$de,
      tested_gene_ids = c("tx2", "tx3", "tx5", "tx6"),
      annotation = f$annotation
    ),
    "significant gene.*not present in the tested universe"
  )
})

test_that("missing and degenerate backgrounds are rejected", {
  f <- gene_set_fixture()
  expect_error(
    prepare_gprofiler_gene_sets(
      f$de[f$de$padj <= 0.05 & !is.na(f$de$padj), ],
      annotation = f$annotation,
      de_table_significant_only = TRUE
    ),
    "tested-gene universe"
  )
  expect_error(
    prepare_gprofiler_gene_sets(
      f$de,
      tested_gene_ids = c("tx1", "tx3"),
      annotation = f$annotation
    ),
    "foreground and background resolve to the same vector"
  )
})

test_that("canonical-name collisions collapse deterministically", {
  f <- gene_set_fixture()
  sets <- prepare_gprofiler_gene_sets(f$de, annotation = f$annotation)
  expect_equal(sum(sets$mapping_table$resolved_name == "GENEB", na.rm = TRUE), 2L)
  expect_equal(sum(sets$mapping_table$included_background & sets$mapping_table$resolved_name == "GENEB"), 1L)
  expect_equal(sets$duplicates_collapsed$background, 1L)
  geneb <- which(!is.na(sets$mapping_table$resolved_name) & sets$mapping_table$resolved_name == "GENEB")
  expect_true(all(sets$mapping_table$duplicate_group[geneb] == "GENEB"))
})

test_that("indexed resolver preserves column priority across duplicate annotation rows", {
  de <- data.frame(
    gene_id = c("tx1", "tx2"),
    significant = c(TRUE, FALSE)
  )
  annotation <- data.frame(
    transcript_id = c("tx1", "tx1", "tx2"),
    sprot_Top_BLASTX_hit = c(NA, "PRIMARY^details", NA),
    EggNM.Preferred_name = c("SECONDARY", "LATER", "BACKGROUND")
  )
  sets <- prepare_gprofiler_gene_sets(de, annotation = annotation)
  tx1 <- sets$mapping_table[sets$mapping_table$original_id == "tx1", ]
  expect_identical(tx1$resolved_name, "PRIMARY")
  expect_identical(tx1$resolution_source, "swissprot_gene_symbol")
})

test_that("raw target IDs never become gProfiler names", {
  de <- data.frame(
    gene_id = c("TRINITY_DN1_c0_g1", "TRINITY_DN2_c0_g1", "nativeA"),
    significant = c(TRUE, FALSE, FALSE)
  )
  annotation <- data.frame(
    transcript_id = de$gene_id,
    sprot_Top_BLASTX_hit = c(NA, "7227.FBpp0070001", NA),
    EggNM.Preferred_name = c("HSP90AA1", NA, NA)
  )
  sets <- prepare_gprofiler_gene_sets(de, annotation = annotation)
  expect_identical(sets$foreground_resolved, "HSP90AA1")
  expect_false(any(grepl("^TRINITY", sets$background_resolved)))
  excluded <- sets$mapping_table[sets$mapping_table$original_id == "TRINITY_DN2_c0_g1", ]
  expect_false(excluded$included_in_gprofiler)
  expect_identical(excluded$exclusion_reason, "no portable canonical name")
  expect_identical(sets$resolver_definition$route, "shared_portable_canonical_organism_context")
})

test_that("mislabeled seed ortholog falls through to eggNOG preferred name", {
  de <- data.frame(gene_id = c("egln3", "actb1"), significant = c(TRUE, FALSE))
  annotation <- data.frame(
    gene_id = de$gene_id,
    sprot_Top_BLASTX_hit = c("7955.ENSDARP00000124991", "7955.ENSDARP00000000001"),
    EggNM.Preferred_name = c("EGLN3", "ACTB")
  )
  sets <- prepare_gprofiler_gene_sets(
    de, annotation = annotation, annotation_id_column = "gene_id",
    use_trinotate_universe = TRUE
  )
  expect_identical(sets$foreground_resolved, "EGLN3")
  expect_identical(sets$background_resolved, c("EGLN3", "ACTB"))
  row <- sets$mapping_table[sets$mapping_table$original_id == "egln3", ]
  expect_true(is.na(row$swissprot_gene_symbol))
  expect_identical(row$portable_name_source, "eggnog_preferred_name")
})

test_that("SwissProt accessions and explicit TrEMBL hits are not submitted as symbols", {
  de <- data.frame(
    gene_id = c("tx1", "tx2", "tx3"),
    significant = c(TRUE, FALSE, FALSE)
  )
  annotation <- data.frame(
    transcript_id = de$gene_id,
    sprot_Top_BLASTX_hit = c(
      "sp|P04637|P53_HUMAN^Cellular tumor antigen p53",
      "P12345^accession without a gene symbol",
      "tr|A0A024RBG1|A0A024RBG1_HUMAN^unreviewed entry"
    ),
    EggNM.Preferred_name = c(NA, "GENEB", "GENEC")
  )
  sets <- prepare_gprofiler_gene_sets(de, annotation = annotation)
  expect_identical(sets$foreground_resolved, "P53")
  expect_identical(sets$background_resolved, c("P53", "GENEB", "GENEC"))
  expect_identical(
    sets$mapping_table$portable_name_source,
    c("swissprot_gene_symbol", "eggnog_preferred_name", "eggnog_preferred_name")
  )
})

test_that("deprecated universe flag does not change submitted vectors", {
  f <- gene_set_fixture()
  without_flag <- prepare_gprofiler_gene_sets(
    f$de, annotation = f$annotation, use_trinotate_universe = FALSE
  )
  with_flag <- prepare_gprofiler_gene_sets(
    f$de, annotation = f$annotation, use_trinotate_universe = TRUE
  )
  expect_identical(with_flag$foreground_resolved, without_flag$foreground_resolved)
  expect_identical(with_flag$background_resolved, without_flag$background_resolved)
  expect_identical(with_flag$mapping_table, without_flag$mapping_table)
})

test_that("accession-like gene symbols remain portable outside accession fields", {
  de <- data.frame(
    gene_id = c("N6AMT2", "zgc195001", "background_gene"),
    significant = c(TRUE, FALSE, FALSE)
  )
  annotation <- data.frame(
    gene_id = de$gene_id,
    sprot_Top_BLASTX_hit = NA_character_,
    EggNM.Preferred_name = c("N6AMT2", "zgc195001", "BACKGROUND")
  )
  sets <- prepare_gprofiler_gene_sets(
    de, annotation = annotation, annotation_id_column = "gene_id"
  )
  expect_identical(
    sets$background_resolved,
    c("N6AMT2", "zgc195001", "BACKGROUND")
  )
})

test_that("portable native symbols support reference-based routes", {
  de <- data.frame(gene_id = c("egln3", "actb1"), significant = c(TRUE, FALSE))
  annotation <- data.frame(
    gene_id = de$gene_id,
    sprot_Top_BLASTX_hit = NA_character_,
    EggNM.Preferred_name = NA_character_
  )
  sets <- prepare_gprofiler_gene_sets(
    de, annotation = annotation, annotation_id_column = "gene_id",
    use_trinotate_universe = TRUE
  )
  expect_identical(sets$foreground_resolved, "egln3")
  expect_true(all(sets$mapping_table$portable_name_source == "native_symbol"))
})

test_that("available eggNOG fields are retained as provenance, not resolver scores", {
  de <- data.frame(gene_id = c("tx1", "tx2"), significant = c(TRUE, FALSE))
  annotation <- data.frame(
    transcript_id = de$gene_id,
    EggNM.Preferred_name = c("GENEA", "GENEB"),
    EggNM.seed_ortholog = c("seedA", "seedB"),
    EggNM.seed_evalue = c("1e-20", "2e-10"),
    EggNM.seed_score = c("200", "100"),
    EggNM.OGs = c("COG0001", "COG0002"),
    EggNM.GOs = c("GO:0000001", "GO:0000002")
  )
  sets <- prepare_gprofiler_gene_sets(de, annotation = annotation)
  row <- sets$mapping_table[sets$mapping_table$original_id == "tx1", ]
  expect_identical(row$eggnog_seed_ortholog, "seedA")
  expect_identical(row$eggnog_seed_evalue, "1e-20")
  expect_identical(row$eggnog_ogs, "COG0001")
  expect_false(any(grepl("score|confidence", names(sets$mapping_table), ignore.case = TRUE) &
                   names(sets$mapping_table) %in% c("consensus_score", "confidence_score")))
})

test_that("Gammarus-like sparse annotation preserves the experimental bottleneck", {
  de <- data.frame(
    gene_id = c("TRINITY_DN1_c0_g1", "TRINITY_DN2_c0_g1", "TRINITY_DN3_c0_g1", "TRINITY_DN4_c0_g1"),
    significant = c(TRUE, FALSE, FALSE, FALSE)
  )
  annotation <- data.frame(
    transcript_id = c("TRINITY_DN1_c0_g1", "TRINITY_DN2_c0_g1"),
    EggNM.Preferred_name = c("GENEA", "GENEB")
  )
  sets <- prepare_gprofiler_gene_sets(de, annotation = annotation)
  expect_identical(sets$foreground_original, "TRINITY_DN1_c0_g1")
  expect_identical(sets$background_original, de$gene_id)
  expect_identical(sets$foreground_resolved, "GENEA")
  expect_identical(sets$background_resolved, c("GENEA", "GENEB"))
  expect_identical(sets$background_unmapped, c("TRINITY_DN3_c0_g1", "TRINITY_DN4_c0_g1"))
})

test_that("logical significance columns reject ambiguous values", {
  f <- gene_set_fixture()
  f$de$significant <- c("yes", "no", "maybe", "no", "no", "no")
  expect_error(
    prepare_gprofiler_gene_sets(f$de, annotation = f$annotation),
    "ambiguous value"
  )
})

make_reference_fixture <- function(root, n = 120L, n_de = 24L) {
  genes <- sprintf("gene%03d", seq_len(n))
  proteins <- sprintf("prot%03d", seq_len(n))

  counts <- data.frame(
    gene_id = genes,
    control_1 = 100L + seq_len(n),
    control_2 = 105L + seq_len(n),
    treated_1 = 110L + seq_len(n),
    check.names = FALSE
  )
  readr::write_tsv(counts, file.path(root, "allcounts_table.txt"))

  dge <- data.frame(
    ID = genes,
    log2FoldChange = rep(2, n),
    padj = rep(0.001, n),
    significant = c(rep(TRUE, n_de), rep(FALSE, n - n_de))
  )
  readr::write_csv(dge, file.path(root, "dge_Treated_vs_Control.csv"))

  starts <- seq(1L, by = 1000L, length.out = n)
  gff <- unlist(lapply(seq_len(n), function(i) {
    end <- starts[[i]] + 299L + i
    c(
      paste("chr1", ".", "gene", starts[[i]], end, ".", "+", ".",
            paste0("ID=gene-", genes[[i]], ";Name=", genes[[i]]), sep = "\t"),
      paste("chr1", ".", "mRNA", starts[[i]], end, ".", "+", ".",
            paste0("ID=rna-", genes[[i]], ";Parent=gene-", genes[[i]]), sep = "\t"),
      paste("chr1", ".", "CDS", starts[[i]], end, ".", "+", "0",
            paste0("ID=cds-", proteins[[i]], ";Parent=rna-", genes[[i]],
                   ";protein_id=", proteins[[i]]), sep = "\t")
    )
  }))
  writeLines(c("##gff-version 3", gff), file.path(root, "genome.gff3"))

  fasta <- unlist(Map(
    function(protein, i) c(
      paste0(">", protein, " synthetic protein"),
      paste0("M", paste(rep("A", 99L + (i %% 10L)), collapse = ""))
    ),
    proteins,
    seq_along(proteins)
  ))
  writeLines(fasta, file.path(root, "proteins.faa"))

  invisible(list(genes = genes, proteins = proteins, n_de = n_de))
}

write_emapper_fixture <- function(path, fixture) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  header <- paste(
    c("#query", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
      "max_annot_lvl", "COG_category", "Description", "Preferred_name", "GOs"),
    collapse = "\t"
  )
  rows <- vapply(seq_along(fixture$proteins), function(i) {
    go <- if (i <= fixture$n_de) "GO:0006355" else "GO:0003674"
    paste(
      fixture$proteins[[i]], "ortholog", "1e-50", "200", "OG", "33208|Metazoa",
      "K", "synthetic", fixture$genes[[i]], go,
      sep = "\t"
    )
  }, character(1))
  writeLines(c("## synthetic eggNOG-mapper output", header, rows), path)
}

test_that("reference preparation is resumable and creates a loadable bundle", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("goseq")

  root <- tempfile("echogo-reference-")
  dir.create(root)
  fixture <- make_reference_fixture(root)

  first <- echogo_prepare_reference_inputs(
    root = root,
    gff_file = "genome.gff3",
    protein_fasta = "proteins.faa",
    reference_label = "Synthetic",
    verbose = FALSE
  )
  expect_identical(first$status, "awaiting_emapper")
  expect_true(file.exists(first$background_fasta))
  expect_true(file.exists(first$emapper_command))
  expect_match(readLines(first$emapper_command), "--itype proteins")

  write_emapper_fixture(first$emapper_file, fixture)
  second <- echogo_prepare_reference_inputs(
    root = root,
    gff_file = "genome.gff3",
    protein_fasta = "proteins.faa",
    reference_label = "Synthetic",
    verbose = FALSE
  )

  expect_identical(second$status, "complete")
  expect_true(all(second$summary$enriched_terms > 0L))
  expected <- c(
    "allcounts_table.txt",
    "dge_Treated_vs_Control.csv",
    "dge_Treated_vs_Control.GOseq.enriched.tsv",
    "Trinotate_for_EchoGO.tsv",
    "Synthetic_eggNOG_for_EchoGO.tsv"
  )
  expect_true(all(file.exists(file.path(second$output_dir, expected))))

  prepared_annotation <- readr::read_tsv(
    second$trinotate_file, show_col_types = FALSE
  )
  expect_true(all(c(
    "gene_id", "native_symbol", "sprot_Top_BLASTX_hit",
    "swissprot_accession", "swissprot_gene_symbol",
    "EggNM.Preferred_name", "EggNM.seed_ortholog",
    "EggNM.max_annot_lvl", "annotation_taxonomy"
  ) %in% names(prepared_annotation)))
  expect_true(all(is.na(prepared_annotation$sprot_Top_BLASTX_hit)))
  expect_true(all(prepared_annotation$EggNM.Preferred_name == prepared_annotation$gene_id))
  expect_true(all(prepared_annotation$EggNM.seed_ortholog == "ortholog"))
  expect_true(all(grepl("Metazoa", prepared_annotation$annotation_taxonomy)))

  annotated <- load_and_annotate_goseq(
    goseq_file = file.path(
      second$output_dir,
      "dge_Treated_vs_Control.GOseq.enriched.tsv"
    ),
    trinotate_file = second$trinotate_file,
    de_file = file.path(second$output_dir, "dge_Treated_vs_Control.csv"),
    count_matrix_file = file.path(second$output_dir, "allcounts_table.txt"),
    output_dir = file.path(root, "annotated")
  )
  expect_gt(nrow(annotated), 0L)
  expect_true(all(c("category", "term", "ontology") %in% names(annotated)))
  expect_true(all(c("total_significant_genes", "total_tested_genes") %in% names(annotated)))
  fold_fixture <- annotated[
    annotated$numDEInCat == 24 & annotated$numInCat == 24,
    "foldEnrichment",
    drop = TRUE
  ]
  if (length(fold_fixture)) expect_equal(fold_fixture[[1]], 5, tolerance = 1e-12)
})

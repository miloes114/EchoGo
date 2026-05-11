test_that("RRvGO mode tables keep strict and conservative sets semantically honest", {
  df <- tibble::tibble(
    term_id = paste0("GO:", sprintf("%07d", 1:6)),
    term_name = paste("term", 1:6),
    ontology = c("BP", "BP", "BP", "MF", "CC", "KEGG"),
    significant_in_any = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
    in_goseq = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
    origin = c(
      "GO terms - Consensus (with BG)",
      "GO terms - GOseq only",
      "GO terms - g:Profiler only (with BG)",
      "GO terms - Consensus (no BG)",
      "GO terms - Consensus (with BG)",
      "KEGG terms - Consensus (with BG)"
    )
  )

  rr_modes <- EchoGO:::.echogo_rrvgo_mode_tables(df)

  expect_setequal(rr_modes$true_consensus_with_bg$term_id, "GO:0000001")
  expect_setequal(
    rr_modes$conservative_bg_supported$term_id,
    c("GO:0000001", "GO:0000002", "GO:0000003")
  )
  expect_setequal(
    rr_modes$exploratory_all_significant$term_id,
    c("GO:0000001", "GO:0000002", "GO:0000003", "GO:0000004")
  )

  expect_false("GO:0000004" %in% rr_modes$true_consensus_with_bg$term_id)
  expect_false("GO:0000004" %in% rr_modes$conservative_bg_supported$term_id)
  expect_false(any(rr_modes$exploratory_all_significant$ontology == "KEGG"))
})

scoreless_output_fixture <- function() {
  tibble::tibble(
    term_id = c("GO:0000001", "GO:0000002", "GO:0000003"),
    term_name = c("target", "alternative", "default only"),
    ontology = c("BP", "BP", "BP"),
    evidence_profile = c("TARGET_ONLY", "ALTERNATIVE_CONTEXT", "NO_PRIMARY_SUPPORT"),
    primary_evidence = c(TRUE, TRUE, FALSE),
    target_goseq_supported = c(TRUE, FALSE, FALSE),
    alternative_context_support_n = c(0L, 1L, 0L),
    alternative_context_support_fraction = c(0, 0.5, 0),
    alternative_queried_context_n = c(2L, 2L, 2L),
    gprof_custom_context_support_n = c(0L, 1L, 0L),
    queried_context_n = c(2L, 2L, 2L),
    default_domain_exploratory_context_support_n = c(0L, 0L, 1L),
    display_order = 1:3,
    representative_order = 1:3,
    contributing_genes = c("A,B,C", "A,B,D", "A,B,E")
  )
}

test_that("RRvGO candidate preparation excludes default-domain-only rows", {
  candidates <- .echogo_prepare_rrvgo_primary_terms(scoreless_output_fixture())
  expect_identical(candidates$term_id, c("GO:0000001", "GO:0000002"))
  expect_true(all(candidates$representative_order_source == "non_inferential_representative_order"))
  expect_false(any(grepl("consensus|pval", names(candidates))))
})

test_that("primary networks exclude default-domain-only rows", {
  output <- tempfile("scoreless-networks-")
  run_all_networks(scoreless_output_fixture(), outdir = output, min_gene_count = 1, min_shared_genes = 1)
  nodes <- readr::read_csv(
    file.path(output, "networks", "primary_custom_evidence", "primary_network_nodes_BP.csv"),
    show_col_types = FALSE
  )
  expect_identical(nodes$term_id, c("GO:0000001", "GO:0000002"))
  expect_false(any(grepl("score", names(nodes))))
})

test_that("scoreless evaluation writes primary descriptive diagnostics only", {
  output <- tempfile("scoreless-evaluation-")
  consensus_dir <- file.path(output, "consensus")
  dir.create(consensus_dir, recursive = TRUE)
  evidence_path <- file.path(consensus_dir, "term_evidence_exact.xlsx")
  openxlsx::write.xlsx(scoreless_output_fixture(), evidence_path)
  readr::write_csv(
    tibble::tibble(
      term_id = "GO:0000002", context_code = "hsapiens", context_label = "human",
      context_role = "ALTERNATIVE", background_mode = "custom_experimental_background",
      source_qualifies = TRUE, submitted_foreground_n = 4, submitted_background_n = 10,
      effective_query_n = 3, effective_background_or_domain_n = 9
    ),
    file.path(consensus_dir, "term_source_provenance_long.csv")
  )
  goseq_path <- file.path(output, "goseq.csv")
  readr::write_csv(tibble::tibble(category = "GO:0000001", over_represented_FDR = 0.01), goseq_path)
  evaluate_consensus_vs_goseq(evidence_path, goseq_path, file.path(output, "evaluation"))
  counts <- readr::read_csv(file.path(output, "evaluation", "primary_evidence_profile_counts.csv"), show_col_types = FALSE)
  expect_identical(counts$evidence_profile, c("ALTERNATIVE_CONTEXT", "TARGET_ONLY"))
  expect_true(file.exists(file.path(output, "evaluation", "README_scoreless_evaluation.txt")))
  expect_false(dir.exists(file.path(output, "evaluation", "exploratory_no_bg")))
})

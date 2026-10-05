evidence_source_row <- function(term_id, source_type, adjusted_p,
                                context_label = NA_character_,
                                context_role = if (identical(source_type, "goseq")) "TARGET_GOSEQ" else "ALTERNATIVE",
                                background_mode = if (identical(source_type, "goseq")) "target_annotation_goseq" else "custom_experimental_background") {
  tibble::tibble(
    term_id = term_id, term_name = paste("term", term_id), ontology = "GO:BP",
    source_type = source_type, context_code = context_label,
    context_label = context_label, context_role = context_role,
    background_mode = background_mode,
    domain_scope = if (identical(background_mode, "custom_experimental_background")) "custom" else "annotated",
    source_adjusted_p = adjusted_p, fold_enrichment = 2,
    contributing_genes = "GENE_A,GENE_B", intersection_size = 2,
    query_size = 4, effective_domain_size = 10,
    submitted_foreground_n = 4, submitted_background_n = 10,
    effective_query_n = 3, effective_background_or_domain_n = 9,
    provenance_file = "fixture.csv"
  )
}

evidence_fixture <- function(target_context = "drerio", include_exploratory = FALSE) {
  contexts <- .echogo_context_configuration(
    c(drerio = "zebrafish", hsapiens = "human", mmusculus = "mouse"), target_context
  )
  goseq <- dplyr::bind_rows(
    evidence_source_row("GO:0000001", "goseq", 0.01),
    evidence_source_row("GO:0000002", "goseq", 0.01),
    evidence_source_row("GO:0000003", "goseq", 0.20)
  )
  custom <- dplyr::bind_rows(
    evidence_source_row("GO:0000001", "gprofiler", 0.01, "zebrafish", "TARGET"),
    evidence_source_row("GO:0000002", "gprofiler", 0.01, "human", "ALTERNATIVE"),
    evidence_source_row("GO:0000003", "gprofiler", 0.02, "human", "ALTERNATIVE")
  )
  exploratory <- if (include_exploratory) dplyr::bind_rows(
    evidence_source_row("GO:0000001", "gprofiler", 1e-9, "mouse", "ALTERNATIVE", "default_domain_exploratory"),
    evidence_source_row("GO:0000004", "gprofiler", 1e-9, "mouse", "ALTERNATIVE", "default_domain_exploratory")
  ) else tibble::tibble()
  assemble_echogo_evidence(goseq, custom, exploratory, contexts)
}

test_that("GOseq support is threshold based and primary profiles are exact", {
  out <- evidence_fixture()
  e <- out$exact_terms
  expect_identical(e$evidence_profile, c("TARGET_ONLY", "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT"))
  expect_identical(e$target_goseq_supported, c(TRUE, TRUE, FALSE))
  expect_identical(e$alternative_context_support_n, c(0L, 1L, 1L))
  expect_false(any(c("consensus_score", "consensus_score_all") %in% names(e)))
  expect_false(any(grepl("min_pval|max_pval|best_p", names(e))))
})

test_that("explicit target contexts are excluded from alternative recurrence", {
  target <- evidence_fixture("drerio")$exact_terms
  all_alternative <- evidence_fixture(NULL)$exact_terms
  expect_equal(target$alternative_context_support_n[target$term_id == "GO:0000001"], 0L)
  expect_equal(all_alternative$alternative_context_support_n[all_alternative$term_id == "GO:0000001"], 1L)
  expect_equal(all_alternative$alternative_queried_context_n, rep(3L, 3L))
})

test_that("researcher-selected contexts are retained and invalid target errors", {
  contexts <- .echogo_context_configuration(c(drerio = "zebrafish", hsapiens = "human"), "drerio")
  expect_identical(contexts$context_code, c("drerio", "hsapiens"))
  expect_identical(contexts$context_role, c("TARGET", "ALTERNATIVE"))
  expect_error(.echogo_context_configuration(c(hsapiens = "human"), "drerio"), "not among")
})

test_that("default-domain rows remain separate from primary profiles and recurrence", {
  primary <- evidence_fixture(include_exploratory = FALSE)$exact_terms
  with_default <- evidence_fixture(include_exploratory = TRUE)$exact_terms
  cols <- c("term_id", "evidence_profile", "alternative_context_support_n", "alternative_context_support_fraction")
  expect_equal(with_default[match(primary$term_id, with_default$term_id), cols], primary[, cols])
  extra <- with_default[with_default$term_id == "GO:0000004", , drop = FALSE]
  expect_identical(extra$evidence_profile, "NO_PRIMARY_SUPPORT")
  expect_identical(extra$default_domain_exploratory_context_support_n, 1L)
})

test_that("context-local p changes do not create an integrated p statistic", {
  first <- evidence_fixture()
  modified <- first$source_provenance
  modified$source_adjusted_p[modified$term_id == "GO:0000003" & modified$context_label == "human"] <- 0.001
  second <- assemble_echogo_evidence(
    modified[modified$source_type == "goseq", ],
    modified[modified$background_mode == "custom_experimental_background", ],
    tibble::tibble(),
    .echogo_context_configuration(c(drerio = "zebrafish", hsapiens = "human", mmusculus = "mouse"), "drerio")
  )
  cols <- c("term_id", "evidence_profile", "alternative_context_support_n", "display_order")
  expect_equal(first$exact_terms[, cols], second$exact_terms[, cols])
  expect_false(any(grepl("pval|adjusted_p.*context", names(second$exact_terms))))
})

test_that("alternative-context drill-down retains source and optional provenance", {
  out <- evidence_fixture()
  provenance <- tibble::tibble(resolved_name = "GENE_A", native_symbol = "nativeA", annotation_match = TRUE)
  drill <- echogo_hypothesis_drilldown(out$exact_terms, out$source_provenance, "GO:0000003", provenance)
  expect_identical(drill$term_evidence$evidence_profile, "ALTERNATIVE_CONTEXT")
  expect_true("human" %in% drill$source_provenance$context_label)
  expect_true("GENE_A" %in% drill$contributing_genes)
  expect_equal(nrow(drill$annotation_provenance), 1L)
})

test_that("unavailable annotation fields are represented without failure", {
  out <- evidence_fixture()
  path <- tempfile("echogo-provenance-")
  .echogo_write_annotation_provenance(tibble::tibble(original_id = "x", resolved_name = "GENE_A"), out$exact_terms, path)
  expect_true(file.exists(file.path(path, "annotation_provenance_gene_level.csv")))
  expect_true(file.exists(file.path(path, "annotation_provenance_summary.csv")))
})

test_that("canonical files produce exact and long-form scoreless outputs", {
  root <- tempfile("scoreless-build-")
  custom <- file.path(root, "gprofiler", "custom_experimental_background")
  dir.create(custom, recursive = TRUE)
  readr::write_tsv(
    tibble::tibble(category = c("GO:0000001", "GO:0000002"),
                   term = c("target", "added"),
                   over_represented_FDR = c(0.01, 0.20)),
    file.path(root, "goseq.tsv")
  )
  readr::write_csv(
    tibble::tibble(term_id = c("GO:0000001", "GO:0000002"), term_name = c("target", "added"),
                   source = "GO:BP", p_value = c(0.01, 0.01), fold_enrichment = 2,
                   intersection = "GENE_A", query_size = 2, effective_domain_size = 10),
    file.path(custom, "gprofiler_human_with_bg.csv")
  )
  jsonlite::write_json(list(submitted_foreground_count = 2, submitted_background_count = 10,
                            effective_query_size = 1, effective_domain_size = 8),
                       file.path(custom, "gprofiler_human_with_bg_metadata.json"), auto_unbox = TRUE)
  out <- build_consensus_table(
    file.path(root, "goseq.tsv"), file.path(root, "gprofiler"), c(hsapiens = "human"),
    file.path(root, "consensus"), annotation_provenance = tibble::tibble(original_id = "x", resolved_name = "GENE_A")
  )
  expect_identical(out$evidence_profile, c("TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT"))
  expect_true(file.exists(file.path(root, "consensus", "term_evidence_exact.csv")))
  expect_true(file.exists(file.path(root, "consensus", "term_source_provenance_long.csv")))
  expect_true(file.exists(file.path(root, "consensus", "annotation_provenance_gene_level.csv")))
  expect_false(any(c("consensus_score", "consensus_score_all") %in% names(out)))
})

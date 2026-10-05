context_rationale_fixture <- function(context_metadata = NULL) {
  species <- c(drerio = "zebrafish", hsapiens = "human", mmusculus = "mouse")
  contexts <- .echogo_context_configuration(species, "drerio", context_metadata)
  source_row <- function(term_id, source_type, adjusted_p, context_label = NA_character_) {
    index <- match(context_label, contexts$context_label)
    tibble::tibble(
      term_id = term_id,
      term_name = term_id,
      ontology = "BP",
      source_type = source_type,
      context_code = if (is.na(index)) NA_character_ else contexts$context_code[[index]],
      context_label = if (is.na(index)) "target_goseq" else contexts$context_label[[index]],
      context_role = if (is.na(index)) "TARGET_GOSEQ" else contexts$context_role[[index]],
      context_rationale = if (is.na(index)) NA_character_ else contexts$context_rationale[[index]],
      selection_basis = if (is.na(index)) NA_character_ else contexts$selection_basis[[index]],
      background_mode = if (identical(source_type, "goseq")) {
        "target_annotation_goseq"
      } else {
        "custom_experimental_background"
      },
      domain_scope = if (identical(source_type, "goseq")) {
        "target_goseq_tested_universe"
      } else {
        "custom"
      },
      source_adjusted_p = adjusted_p,
      fold_enrichment = 2,
      contributing_genes = "GENE_A,GENE_B",
      intersection_size = 2,
      query_size = 4,
      effective_domain_size = 10,
      submitted_foreground_n = 4,
      submitted_background_n = 10,
      effective_query_n = 3,
      effective_background_or_domain_n = 9,
      provenance_file = "fixture.csv"
    )
  }
  goseq <- dplyr::bind_rows(
    source_row("GO:0000001", "goseq", 0.01),
    source_row("GO:0000002", "goseq", 0.20),
    source_row("GO:0000003", "goseq", 0.01)
  )
  custom <- dplyr::bind_rows(
    source_row("GO:0000001", "gprofiler", 0.01, "human"),
    source_row("GO:0000002", "gprofiler", 0.01, "human"),
    source_row("GO:0000003", "gprofiler", 0.01, "zebrafish")
  )
  assemble_echogo_evidence(goseq, custom, context_configuration = contexts)
}

write_context_rationale_builder_fixture <- function(root) {
  custom <- file.path(root, "gprofiler", "custom_experimental_background")
  dir.create(custom, recursive = TRUE)
  readr::write_tsv(
    tibble::tibble(
      category = c("GO:0000001", "GO:0000002"),
      term = c("target", "alternative"),
      over_represented_FDR = c(0.01, 0.20)
    ),
    file.path(root, "goseq.tsv")
  )
  readr::write_csv(
    tibble::tibble(
      term_id = c("GO:0000001", "GO:0000002"),
      term_name = c("target", "alternative"),
      source = "GO:BP",
      p_value = c(0.01, 0.01),
      fold_enrichment = 2,
      intersection = "GENE_A",
      query_size = 2,
      effective_domain_size = 10
    ),
    file.path(custom, "gprofiler_human_with_bg.csv")
  )
  jsonlite::write_json(
    list(
      submitted_foreground_count = 2,
      submitted_background_count = 10,
      effective_query_size = 1,
      effective_domain_size = 8
    ),
    file.path(custom, "gprofiler_human_with_bg_metadata.json"),
    auto_unbox = TRUE
  )
}

test_that("context rationale is optional and preserves existing evidence semantics", {
  species <- c(drerio = "zebrafish", hsapiens = "human", mmusculus = "mouse")
  without_metadata <- .echogo_context_configuration(species, "drerio")
  empty_metadata <- .echogo_context_configuration(
    species,
    "drerio",
    data.frame(context_code = character(), context_rationale = character())
  )
  expect_equal(without_metadata, empty_metadata)
  expect_true(all(is.na(without_metadata$context_rationale)))
  expect_true(all(is.na(without_metadata$selection_basis)))

  bare <- context_rationale_fixture()
  annotated <- context_rationale_fixture(data.frame(
    context_code = c("drerio", "hsapiens"),
    context_rationale = c("experimental model", "ecologically relevant reference"),
    selection_basis = c("experimental_relevance", "ecological_relevance")
  ))
  exact_fields <- c(
    "term_id", "evidence_profile", "target_goseq_supported",
    "alternative_context_support_n", "alternative_context_support_fraction",
    "gprof_custom_context_support_n", "display_order", "representative_order"
  )
  expect_equal(bare$exact_terms[, exact_fields], annotated$exact_terms[, exact_fields])
  expect_equal(bare$source_provenance$source_qualifies, annotated$source_provenance$source_qualifies)

  bare_rrvgo <- .echogo_prepare_rrvgo_primary_terms(bare$exact_terms)
  annotated_rrvgo <- .echogo_prepare_rrvgo_primary_terms(annotated$exact_terms)
  expect_equal(bare_rrvgo, annotated_rrvgo)
})

test_that("context rationale round-trips through canonical configuration and long provenance", {
  root <- tempfile("context-rationale-build-")
  write_context_rationale_builder_fixture(root)
  metadata <- data.frame(
    context_code = "hsapiens",
    context_rationale = "Ecologically relevant\t reference fish",
    selection_basis = "ecological relevance",
    stringsAsFactors = FALSE
  )
  out <- build_consensus_table(
    goseq_file = file.path(root, "goseq.tsv"),
    gprofiler_dir = file.path(root, "gprofiler"),
    species_map = c(hsapiens = "human"),
    output_dir = file.path(root, "consensus"),
    context_metadata = metadata
  )
  configuration <- readr::read_csv(
    file.path(root, "consensus", "context_configuration.csv"), show_col_types = FALSE
  )
  provenance <- readr::read_csv(
    file.path(root, "consensus", "term_source_provenance_long.csv"), show_col_types = FALSE
  )
  expect_identical(configuration$context_rationale, "Ecologically relevant reference fish")
  expect_identical(configuration$selection_basis, "ecological relevance")
  human_rows <- provenance[!is.na(provenance$context_code) & provenance$context_code == "hsapiens", , drop = FALSE]
  expect_true(all(human_rows$context_rationale == "Ecologically relevant reference fish"))
  expect_true(all(human_rows$selection_basis == "ecological relevance"))

  drill <- echogo_hypothesis_drilldown(out, attr(out, "source_provenance"), "GO:0000002")
  expect_identical(drill$context_selection$context_label, "human")
  expect_identical(drill$context_selection$context_rationale, "Ecologically relevant reference fish")
})

test_that("target roles, network outputs, and evaluation outputs ignore context rationale", {
  bare <- context_rationale_fixture()
  annotated <- context_rationale_fixture(data.frame(
    context_label = "human",
    context_rationale = "comparative ecological hypothesis",
    stringsAsFactors = FALSE
  ))
  bare_target <- bare$exact_terms[bare$exact_terms$term_id == "GO:0000003", , drop = FALSE]
  annotated_target <- annotated$exact_terms[annotated$exact_terms$term_id == "GO:0000003", , drop = FALSE]
  expect_equal(bare_target$alternative_context_support_n, 0L)
  expect_equal(bare_target$alternative_context_support_n, annotated_target$alternative_context_support_n)

  network_files <- function(evidence) {
    root <- tempfile("context-rationale-network-")
    network_evidence <- dplyr::mutate(
      evidence,
      primary_evidence = .data$evidence_profile != "NO_PRIMARY_SUPPORT"
    )
    run_all_networks(network_evidence, outdir = root, min_gene_count = 1, min_shared_genes = 1)
    readr::read_csv(
      file.path(root, "networks", "primary_custom_evidence", "primary_network_nodes_BP.csv"),
      show_col_types = FALSE
    )
  }
  expect_equal(network_files(bare$exact_terms), network_files(annotated$exact_terms))

  evaluation_files <- function(evidence, provenance) {
    root <- tempfile("context-rationale-evaluation-")
    consensus_dir <- file.path(root, "consensus")
    dir.create(consensus_dir, recursive = TRUE)
    consensus_file <- file.path(consensus_dir, "term_evidence_exact.xlsx")
    openxlsx::write.xlsx(evidence, consensus_file)
    readr::write_csv(provenance, file.path(consensus_dir, "term_source_provenance_long.csv"))
    goseq_file <- file.path(root, "goseq.csv")
    readr::write_csv(
      tibble::tibble(category = "GO:0000001", over_represented_FDR = 0.01),
      goseq_file
    )
    evaluation_dir <- file.path(root, "evaluation")
    evaluate_consensus_vs_goseq(consensus_file, goseq_file, evaluation_dir)
    list(
      profiles = readr::read_csv(file.path(evaluation_dir, "primary_evidence_profile_counts.csv"), show_col_types = FALSE),
      recurrence = readr::read_csv(file.path(evaluation_dir, "primary_annotation_context_recurrence.csv"), show_col_types = FALSE),
      coverage = readr::read_csv(file.path(evaluation_dir, "custom_background_context_recognition_coverage.csv"), show_col_types = FALSE)
    )
  }
  expect_equal(
    evaluation_files(bare$exact_terms, bare$source_provenance),
    evaluation_files(annotated$exact_terms, annotated$source_provenance)
  )
})

test_that("missing or invalid context rationale metadata has an explicit contract", {
  contexts <- .echogo_context_configuration(
    c(hsapiens = "human", mmusculus = "mouse"),
    context_metadata = data.frame(context_label = "human", stringsAsFactors = FALSE)
  )
  expect_true(is.na(contexts$context_rationale[contexts$context_label == "human"]))
  expect_error(
    .echogo_context_configuration(
      c(hsapiens = "human"),
      context_metadata = data.frame(context_code = "rnorvegicus", context_rationale = "not queried")
    ),
    "unqueried context"
  )
  expect_error(
    .echogo_context_configuration(
      c(hsapiens = "human", mmusculus = "mouse"),
      context_metadata = data.frame(context_code = "hsapiens", context_label = "mouse")
    ),
    "same queried context"
  )
})

test_that("Decision 0008 partitions only the intended semantic products", {
  x <- tibble::tibble(
    term_id = sprintf("GO:%07d", 1:7),
    ontology = c("BP", "MF", "CC", "BP", "BP", "BP", "BP"),
    evidence_profile = c(
      "TARGET_ONLY", "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT",
      "NO_PRIMARY_SUPPORT", "TARGET_ONLY", "ALTERNATIVE_CONTEXT", "TARGET_ONLY"
    ),
    primary_evidence = c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE),
    target_goseq_supported = c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
    alternative_context_support_n = c(0L, 1L, 3L, 0L, 0L, 2L, 0L),
    alternative_context_recurrence = c(0L, 1L, 3L, 0L, 0L, 2L, 0L),
    target_goseq_adjusted_p = c(.001, .2, .3, .4, .005, .6, .7),
    source_specific_p_value = c(.01, .02, .03, .04, .05, .06, .07),
    source_p_value_target = c(.01, .02, .03, .04, .05, .06, .07),
    source_p_adjusted_target = c(.1, .2, .3, .4, .5, .6, .7),
    submitted_foreground_vector = c("A", "B", "C", "D", "E", "F", "G"),
    submitted_background_vector = rep("BG", 7),
    context_role = c("TARGET", "TARGET", "ALTERNATIVE", "NONE", "TARGET", "ALTERNATIVE", "TARGET"),
    context_rationale = paste("rationale", 1:7),
    representative_order = 1:7
  )

  parts <- .echogo_rrvgo_partition_inputs(x)
  expect_named(parts, c("target_supported", "alternative_context_hypothesis"))
  expect_identical(parts$target_supported$term_id, sprintf("GO:%07d", c(1, 2, 5, 7)))
  expect_identical(parts$alternative_context_hypothesis$term_id, "GO:0000003")
  expect_false(any(parts$target_supported$evidence_profile == "ALTERNATIVE_CONTEXT"))
  expect_false(any(parts$alternative_context_hypothesis$evidence_profile != "ALTERNATIVE_CONTEXT"))

  expected <- x[x$primary_evidence &
    x$evidence_profile %in% c("TARGET_ONLY", "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT") &
    x$ontology %in% c("BP", "MF", "CC"), , drop = FALSE]
  recombined <- dplyr::bind_rows(parts) |>
    dplyr::arrange(.data$term_id)
  expected <- expected |>
    dplyr::arrange(.data$term_id)

  invariant_columns <- c(
    "term_id", "target_goseq_supported", "evidence_profile",
    "alternative_context_support_n", "alternative_context_recurrence",
    "target_goseq_adjusted_p", "source_specific_p_value",
    "source_p_value_target", "source_p_adjusted_target",
    "submitted_foreground_vector", "submitted_background_vector",
    "context_role", "context_rationale"
  )
  expect_equal(
    dplyr::select(recombined, dplyr::all_of(invariant_columns)),
    dplyr::select(expected, dplyr::all_of(invariant_columns))
  )
  expect_true(all(recombined$term_id %in% expected$term_id))
})

test_that("Decision 0008 partition rejects malformed inputs", {
  expect_error(.echogo_rrvgo_partition_inputs(tibble::tibble(term_id = "GO:0000001")),
               "requires columns")
  expect_error(.echogo_rrvgo_partition_inputs(list()), "must be a data frame")
})

test_that("own-data context selection has no arbitrary package default", {
  withr::local_options(EchoGO.default_species = NULL)

  expect_error(
    .echogo_resolve_enrichment_contexts(
      NULL, argument_missing = TRUE, caller = "test"
    ),
    "requires researcher-selected enrichment contexts"
  )

  withr::local_options(EchoGO.default_species = c("drerio", "strutta"))
  expect_warning(
    contexts <- .echogo_resolve_enrichment_contexts(
      NULL, argument_missing = TRUE, caller = "test"
    ),
    "compatibility fallback"
  )
  expect_identical(contexts, c("drerio", "strutta"))
  expect_null(getOption("EchoGO.default_orgdb"))
})

test_that("target context distinguishes target, explicit no-target and compatibility", {
  target <- .echogo_resolve_target_declaration("drerio", caller = "test")
  expect_identical(target$value, "drerio")
  expect_identical(target$state, "TARGET_DECLARED")

  no_target <- .echogo_resolve_target_declaration(NA_character_, caller = "test")
  expect_null(no_target$value)
  expect_true(is.na(no_target$configured_value))
  expect_identical(no_target$state, "NO_TARGET_DECLARED")

  yaml_no_target <- .echogo_resolve_target_declaration("NO_TARGET", caller = "test")
  expect_identical(yaml_no_target$state, "NO_TARGET_DECLARED")

  expect_warning(
    compatibility <- .echogo_resolve_target_declaration(
      NULL, argument_missing = TRUE, caller = "test"
    ),
    "v0.1.3 compatibility"
  )
  expect_identical(compatibility$state, "COMPATIBILITY_NO_TARGET")
})

test_that("semantic reference preferred and compatibility aliases resolve deterministically", {
  preferred <- .echogo_resolve_semantic_reference(
    semantic_reference_orgdb = "org.Dr.eg.db",
    semantic_reference_role = "target_reference",
    run_rrvgo = TRUE,
    caller = "test"
  )
  expect_identical(preferred$orgdb, "org.Dr.eg.db")
  expect_identical(preferred$role, "target_reference")
  expect_identical(preferred$source, "semantic_reference_orgdb")

  expect_warning(
    legacy <- .echogo_resolve_semantic_reference(
      orgdb = "org.Dr.eg.db",
      semantic_reference_role = "target_reference",
      run_rrvgo = TRUE,
      caller = "test"
    ),
    "compatibility alias"
  )
  expect_identical(legacy$orgdb, preferred$orgdb)
  expect_identical(legacy$role, preferred$role)

  both <- .echogo_resolve_semantic_reference(
    semantic_reference_orgdb = "org.Dr.eg.db",
    orgdb = "org.Dr.eg.db",
    semantic_reference_role = "target_reference",
    run_rrvgo = TRUE,
    caller = "test"
  )
  expect_identical(both$orgdb, "org.Dr.eg.db")

  expect_error(
    .echogo_resolve_semantic_reference(
      semantic_reference_orgdb = "org.Dr.eg.db",
      orgdb = "org.Dm.eg.db",
      semantic_reference_role = "proxy",
      run_rrvgo = TRUE,
      caller = "test"
    ),
    "Conflicting RRvGO semantic references"
  )
})

test_that("RRvGO has no silent reference or role", {
  expect_error(
    .echogo_resolve_semantic_reference(
      semantic_reference_role = "proxy", run_rrvgo = TRUE, caller = "test"
    ),
    "no semantic reference was declared"
  )
  expect_error(
    .echogo_resolve_semantic_reference(
      semantic_reference_orgdb = "org.Dm.eg.db",
      run_rrvgo = TRUE,
      caller = "test"
    ),
    "semantic_reference_role is missing"
  )
  expect_error(
    .echogo_resolve_semantic_reference(
      semantic_reference_orgdb = "org.Dm.eg.db",
      semantic_reference_role = "nearest_species",
      run_rrvgo = TRUE,
      caller = "test"
    ),
    "should be one of"
  )
})

test_that("high-level preferred and legacy semantic names pass the same reference", {
  captured <- list()
  testthat::local_mocked_bindings(
    echogo_preflight_species = function(species, ...) species,
    .echogo_preflight_semantic_dependencies = function(...) invisible(TRUE),
    run_echogo_pipeline = function(...) {
      captured[[length(captured) + 1L]] <<- list(...)
      list(files = list(), dirs = list())
    },
    .package = "EchoGO"
  )

  preferred <- run_full_echogo(
    species = "drerio",
    target_context = "drerio",
    semantic_reference_orgdb = "org.Dr.eg.db",
    semantic_reference_role = "target_reference",
    run_rrvgo = TRUE,
    make_report = FALSE,
    outdir = tempfile("preferred-config-")
  )
  expect_warning(
    legacy <- run_full_echogo(
      species = "drerio",
      target_context = "drerio",
      orgdb = "org.Dr.eg.db",
      semantic_reference_role = "target_reference",
      run_rrvgo = TRUE,
      make_report = FALSE,
      outdir = tempfile("legacy-config-")
    ),
    "compatibility alias"
  )

  expect_identical(captured[[1]]$semantic_reference_orgdb, "org.Dr.eg.db")
  expect_identical(captured[[2]]$semantic_reference_orgdb, "org.Dr.eg.db")
  expect_null(captured[[1]]$orgdb)
  expect_null(captured[[2]]$orgdb)
  expect_type(preferred, "list")
  expect_type(legacy, "list")
})

test_that("semantic API naming cannot change exact-term evidence or partitions", {
  evidence <- tibble::tibble(
    term_id = sprintf("GO:%07d", 1:4),
    ontology = "BP",
    evidence_profile = c(
      "TARGET_ONLY", "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT",
      "NO_PRIMARY_SUPPORT"
    ),
    primary_evidence = c(TRUE, TRUE, TRUE, FALSE),
    target_goseq_supported = c(TRUE, TRUE, FALSE, FALSE),
    alternative_context_support_n = c(0L, 1L, 1L, 0L),
    target_goseq_adjusted_p = c(0.01, 0.02, NA, NA)
  )
  original <- evidence

  preferred <- .echogo_resolve_semantic_reference(
    semantic_reference_orgdb = "org.Dr.eg.db",
    semantic_reference_role = "target_reference"
  )
  legacy <- suppressWarnings(.echogo_resolve_semantic_reference(
    orgdb = "org.Dr.eg.db",
    semantic_reference_role = "target_reference"
  ))
  preferred_parts <- .echogo_rrvgo_partition_inputs(evidence)
  legacy_parts <- .echogo_rrvgo_partition_inputs(evidence)

  expect_identical(evidence, original)
  expect_identical(preferred$orgdb, legacy$orgdb)
  expect_identical(preferred_parts, legacy_parts)
})

test_that("default-domain exploration remains opt-in", {
  expect_false(eval(formals(run_full_echogo)$run_exploratory_default_domain))
  expect_false(eval(formals(run_echogo_pipeline)$run_exploratory_default_domain))
  expect_false(eval(formals(run_gprofiler_cross_species)$run_exploratory_default_domain))
})

test_that("scaffold contains declarations but no arbitrary biological defaults", {
  root <- tempfile("echogo-scaffold-contract-")
  echogo_scaffold(root)
  config <- readLines(file.path(root, "input", "config.yml"), warn = FALSE)

  expect_true(any(grepl("^species: \\[\\]$", config)))
  expect_true(any(grepl("^target_context: null$", config)))
  expect_true(any(grepl("^run_exploratory_default_domain: false$", config)))
  expect_true(any(grepl("^run_rrvgo: false$", config)))
  expect_true(any(grepl("^semantic_reference_orgdb: null$", config)))
  expect_true(any(grepl("^semantic_reference_role: null$", config)))
  expect_false(any(grepl("hsapiens|mmusculus|drerio|org\\.[A-Za-z]+\\.eg\\.db", config)))
})

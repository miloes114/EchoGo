# Public configuration-contract helpers --------------------------------------

.echogo_resolve_enrichment_contexts <- function(species, argument_missing = FALSE,
                                                caller = "EchoGO") {
  if (isTRUE(argument_missing) || is.null(species) || !length(species)) {
    legacy <- getOption("EchoGO.default_species", NULL)
    if (!is.null(legacy) && length(legacy)) {
      warning(
        caller, " is using the user-configured EchoGO.default_species option as a ",
        "v0.1.3 compatibility fallback. Supply species = ... explicitly; EchoGO ",
        "does not choose annotation contexts for own-data analyses.",
        call. = FALSE
      )
      species <- legacy
    } else {
      stop(
        caller, " requires researcher-selected enrichment contexts. Supply ",
        "species = c(...) or species_expr = .... EchoGO does not choose a ",
        "biological context panel for own-data analyses.",
        call. = FALSE
      )
    }
  }

  species_names <- names(species)
  species <- as.character(species)
  if (!is.null(species_names)) names(species) <- species_names
  species <- species[!is.na(species) & nzchar(trimws(species))]
  if (!length(species)) {
    stop(caller, " requires at least one non-empty enrichment context.", call. = FALSE)
  }
  species[!duplicated(species)]
}

.echogo_resolve_target_declaration <- function(target_context,
                                               argument_missing = FALSE,
                                               caller = "EchoGO") {
  compatibility <- isTRUE(argument_missing) || is.null(target_context) || !length(target_context)
  if (compatibility) {
    warning(
      caller, " received no explicit target-context declaration. For v0.1.3 ",
      "compatibility this run is treated as having no target g:Profiler context. ",
      "Use target_context = NA_character_ to declare that state explicitly, or ",
      "supply a queried context code/label.",
      call. = FALSE
    )
    return(list(
      value = NULL,
      configured_value = NA_character_,
      state = "COMPATIBILITY_NO_TARGET"
    ))
  }

  if (length(target_context) != 1L) {
    stop(
      "target_context must be one queried context code/label or NA_character_ ",
      "for an explicit no-target declaration.",
      call. = FALSE
    )
  }

  if (is.na(target_context) || identical(toupper(trimws(as.character(target_context))), "NO_TARGET")) {
    return(list(
      value = NULL,
      configured_value = NA_character_,
      state = "NO_TARGET_DECLARED"
    ))
  }

  target_context <- trimws(as.character(target_context))
  if (!nzchar(target_context)) {
    stop(
      "target_context cannot be empty. Supply a queried context code/label or ",
      "NA_character_ for an explicit no-target declaration.",
      call. = FALSE
    )
  }

  list(
    value = target_context,
    configured_value = target_context,
    state = "TARGET_DECLARED"
  )
}

.echogo_normalize_semantic_orgdb <- function(x, argument) {
  if (is.null(x) || !length(x)) return(NULL)
  x <- as.character(x)
  x <- trimws(x[!is.na(x) & nzchar(trimws(x))])
  if (!length(x)) return(NULL)
  if (length(x) != 1L) {
    stop(argument, " must name exactly one semantic-reference OrgDb package.", call. = FALSE)
  }
  x
}

.echogo_resolve_semantic_reference <- function(
    semantic_reference_orgdb = NULL,
    orgdb = NULL,
    semantic_reference_role = NULL,
    run_rrvgo = TRUE,
    caller = "EchoGO",
    warn_legacy = TRUE
) {
  preferred <- .echogo_normalize_semantic_orgdb(
    semantic_reference_orgdb, "semantic_reference_orgdb"
  )
  legacy <- .echogo_normalize_semantic_orgdb(orgdb, "orgdb")

  if (!is.null(preferred) && !is.null(legacy) && !identical(preferred, legacy)) {
    stop(
      "Conflicting RRvGO semantic references: semantic_reference_orgdb = '",
      preferred, "' but legacy orgdb = '", legacy, "'. Supply one value or ",
      "make both identical.",
      call. = FALSE
    )
  }

  resolved <- preferred
  source <- if (!is.null(preferred)) "semantic_reference_orgdb" else NULL
  if (is.null(resolved) && !is.null(legacy)) {
    resolved <- legacy
    source <- "orgdb_compatibility_alias"
    if (isTRUE(warn_legacy)) {
      warning(
        caller, ": orgdb is a v0.1.4 compatibility alias for the RRvGO ",
        "semantic reference. Use semantic_reference_orgdb = '", resolved,
        "'. This database is not an enrichment context.",
        call. = FALSE
      )
    }
  }

  role <- NULL
  if (!is.null(semantic_reference_role) && length(semantic_reference_role)) {
    if (length(semantic_reference_role) != 1L || is.na(semantic_reference_role)) {
      stop(
        "semantic_reference_role must be exactly 'target_reference' or 'proxy'.",
        call. = FALSE
      )
    }
    role <- match.arg(
      as.character(semantic_reference_role),
      c("target_reference", "proxy")
    )
  }

  if (isTRUE(run_rrvgo) && is.null(resolved)) {
    stop(
      caller, " requested RRvGO but no semantic reference was declared. Supply ",
      "semantic_reference_orgdb with the target-organism OrgDb when suitable, ",
      "or a biologically defensible proxy. This OrgDb is used only for semantic ",
      "similarity and is not another enrichment context.",
      call. = FALSE
    )
  }
  if (isTRUE(run_rrvgo) && is.null(role)) {
    stop(
      caller, " requested RRvGO but semantic_reference_role is missing. Declare ",
      "semantic_reference_role = 'target_reference' or 'proxy'; EchoGO does not ",
      "infer this role from taxonomy or the enrichment context panel.",
      call. = FALSE
    )
  }

  list(orgdb = resolved, role = role, source = source)
}

.echogo_preflight_semantic_dependencies <- function(semantic_reference_orgdb) {
  required <- c("rrvgo", "GO.db", semantic_reference_orgdb)
  required <- unique(required[!is.na(required) & nzchar(required)])
  available <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
  if (all(available)) return(invisible(TRUE))

  missing <- required[!available]
  stop(
    "RRvGO semantic reduction requires installed package(s): ",
    paste(missing, collapse = ", "), ".\nInstall them explicitly, for example:\n",
    "  BiocManager::install(c(",
    paste(sprintf("'%s'", missing), collapse = ", "),
    "))\nEchoGO does not install optional semantic dependencies unattended.",
    call. = FALSE
  )
}

.onAttach <- function(libname, pkgname) {
  version <- tryCatch(
    as.character(utils::packageVersion(pkgname)),
    error = function(e) ""
  )
  version_label <- if (nzchar(version)) paste0(" v", version) else ""
  # Let cli detect console colour support (including NO_COLOR). In monochrome
  # sessions, and when this optional dependency is absent, keep the same layout.
  bold <- teal <- amber <- identity
  if (requireNamespace("cli", quietly = TRUE)) {
    bold <- cli::style_bold
    teal <- cli::make_ansi_style("#1294A5")
    amber <- cli::make_ansi_style("#C67500")
  }
  rule <- paste0("  ", strrep("-", 74L))

  startup <- c(
    "",
    paste0("  ", bold(paste0("Echo", teal("GO"))), version_label),
    "  Functional interpretation across annotation contexts",
    teal(rule),
    "",
    "  Explore how annotation context shapes functional interpretation.",
    "  Target GOseq stays primary; every GO term retains its source and context.",
    "",
    paste0("  ", bold(teal("EVIDENCE"))),
    "  Target only | Target + context | Context-derived hypotheses",
    "",
    paste0("  ", bold(amber("FULL OFFLINE DEMO"))),
    "    echogo_quickstart(run_demo = TRUE, full = TRUE)",
    "    Cached g:Profiler; builds semantic summaries, networks and HTML report.",
    "",
    paste0("  ", bold(teal("YOUR EXPERIMENT"))),
    sprintf("    %-30s # %s", 'echogo_scaffold("my_project")', "Create an input template"),
    sprintf("    %-30s # %s", "echogo_pick_species()", "Browse annotation contexts"),
    "    Inputs: DE results, tested-feature background, annotation and target GOseq.",
    "",
    paste0("  ", bold(teal("GUIDES"))),
    sprintf("    %-30s # %s", "echogo_help()", "Workflow and commands"),
    sprintf("    %-30s # %s", 'browseVignettes("EchoGO")', "Browse installed guides"),
    '    vignette("reference-based-inputs")',
    "",
    "  Reports: results/report/ inside your demo or analysis folder",
    '  Cite:    citation("EchoGO")',
    "  Quiet:   suppressPackageStartupMessages(library(EchoGO))",
    teal(rule)
  )

  packageStartupMessage(paste(startup, collapse = "\n"))

  # Offline indicator (optional but helpful)
  if (!getOption("EchoGO.taxonomy_online", TRUE)) {
    packageStartupMessage("Taxonomy enrichment offline: using cached or fallback ranks.\n")
  }
}


.onLoad <- function(...) {
  op <- options()

  # Suggest a per-user demo root and create it only when needed.
  demo_root_default <- tryCatch({
    if (requireNamespace("rappdirs", quietly = TRUE)) {
      rappdirs::user_data_dir("EchoGO")
    } else {
      file.path(tempdir(), "EchoGO_demo")
    }
  }, error = function(e) file.path(tempdir(), "EchoGO_demo"))

  op.echogo <- list(
    # ---- defaults ----
    EchoGO.legacy_aliases      = FALSE,
    # ---- online/offline knobs ----
    EchoGO.gprofiler_base      = "https://biit.cs.ut.ee/gprofiler",
    EchoGO.species_autoupdate  = TRUE,
    EchoGO.taxonomy_online     = TRUE,
    EchoGO.species_max_age_days= 30,
    EchoGO.net_timeout_sec     = 25,
    EchoGO.debug               = FALSE,

    # ---- Per-user demo target -----------------------------------------------
    EchoGO.demo_root           = demo_root_default
  )

  toset <- !(names(op.echogo) %in% names(op))
  if (any(toset)) options(op.echogo[toset])
  invisible()
}


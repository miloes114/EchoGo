.onAttach <- function(libname, pkgname) {

  banner <- "
EchoGO - matched organism-context enrichment
  "

  packageStartupMessage(
    paste0(
      banner,
      "\nWelcome to EchoGO, a toolkit for background-aware and exploratory functional enrichment.",
      "\nHow to cite: citation('EchoGO')",
      "\nQuick start: echogo_help(); echogo_quickstart(run_demo = TRUE)",
      "\nPick species interactively: echogo_pick_species()",
      "\nMore documentation: vignette('EchoGO_workflow'), vignette('EchoGO_interpretation')",
      "\nTo silence this message: suppressPackageStartupMessages(library(EchoGO))\n"
    )
  )

  # Show default OrgDb & species
  try({
    active_orgdb   <- getOption("EchoGO.default_orgdb", "org.Mm.eg.db")
    active_species <- getOption("EchoGO.default_species")

    packageStartupMessage(
      paste0(
        "\nDefault OrgDb: ", paste(active_orgdb, collapse = ", "),
        "\nDefault organism contexts: ", paste(active_species, collapse = ", "),
        "\nTip: change with options(EchoGO.default_orgdb=...), options(EchoGO.default_species=...)\n"
      )
    )
  }, silent = TRUE)

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
    EchoGO.default_species     = c("hsapiens","mmusculus","drerio","btaurus","rnorvegicus","dmelanogaster","ggallus"),
    EchoGO.default_orgdb       = "org.Mm.eg.db",

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


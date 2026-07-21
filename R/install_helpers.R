# --- R/install_helpers.R -----------------------------------------------------

#' Install GO.db / OrgDb packages into the active R library (no prompts)
#'
#' Installs Bioconductor annotation packages to the first active library, which
#' respects project-local libraries such as `renv`. If `pkgs = NULL`, it
#' installs any **missing** items from a standard OrgDb set plus GO.db. Returns
#' a logical vector (per package) of whether the package is now loadable.
#'
#' @param pkgs Character vector of package names (e.g., c("GO.db","org.Mm.eg.db")).
#'   If NULL, installs missing ones from the default set:
#'   c("GO.db","org.Hs.eg.db","org.Mm.eg.db","org.Dr.eg.db",
#'     "org.Dm.eg.db","org.Rn.eg.db","org.Ce.eg.db").
#' @param update Logical; update existing packages? Default FALSE for reproducibility.
#' @param lib Character scalar; library to install into. Defaults to
#'   `.libPaths()[1]`, including the active `renv` project library when present.
#' @return Named logical vector: TRUE if package is loadable after installation.
#' @export
echogo_install_orgdb <- function(
    pkgs   = NULL,
    update = FALSE,
    lib    = .libPaths()[1]
) {
  # Default set if not provided
  default_set <- c(
    "GO.db",
    "org.Hs.eg.db","org.Mm.eg.db","org.Dr.eg.db",
    "org.Dm.eg.db","org.Rn.eg.db","org.Ce.eg.db"
  )

  requested <- if (is.null(pkgs)) default_set else unique(as.character(pkgs))
  requested <- requested[!is.na(requested) & nzchar(requested)]

  installed <- vapply(requested, requireNamespace, logical(1), quietly = TRUE)
  to_install <- if (isTRUE(update)) requested else requested[!installed]

  # Nothing to do?
  if (!length(to_install)) {
    now <- vapply(requested, requireNamespace, logical(1), quietly = TRUE)
    names(now) <- requested
    return(now)
  }

  if (length(lib) != 1L || is.na(lib) || !nzchar(lib)) {
    stop("'lib' must be a single non-empty library path.", call. = FALSE)
  }

  # Ensure BiocManager + repositories, while preserving the active library.
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    utils::install.packages("BiocManager", lib = lib)
  }
  options(repos = BiocManager::repositories())

  # Install
  BiocManager::install(
    to_install,
    lib = lib,
    ask = FALSE,
    update = isTRUE(update)
  )

  # Report loadability after install
  out <- vapply(requested, requireNamespace, logical(1), quietly = TRUE)
  names(out) <- requested
  return(out)
}

#' Ensure OrgDb packages are available (install if missing)
#'
#' This checks that the requested OrgDb packages are available. If they are
#' missing and `auto_install = TRUE`, it calls [echogo_install_orgdb()] to
#' install them into the active library.
#'
#' @param pkgs Character vector of OrgDb names. If NULL, uses
#'   getOption(\"EchoGO.default_orgdb\", \"org.Mm.eg.db\").
#' @param auto_install Logical; if TRUE, will attempt installation
#'   for missing packages.
#' @return Character vector of packages that are loadable now (possibly a subset of `pkgs`).
#' @export
echogo_require_orgdb <- function(pkgs = NULL, auto_install = TRUE) {
  if (is.null(pkgs)) {
    pkgs <- getOption("EchoGO.default_orgdb", "org.Mm.eg.db")
  }

  pkgs <- as.character(pkgs)
  pkgs <- pkgs[nzchar(pkgs)]

  have <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (any(!have) && isTRUE(auto_install)) {
    message("Installing missing OrgDb: ", paste(pkgs[!have], collapse = ", "))
    echogo_install_orgdb(pkgs[!have], update = FALSE)
    have <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  }

  if (any(!have)) {
    warning(
      "Some OrgDb are still unavailable: ",
      paste(pkgs[!have], collapse = ", "),
      "\nTip: EchoGO::echogo_install_orgdb(c(",
      paste(sprintf("'%s'", pkgs[!have]), collapse = ", "),
      "))"
    )
  }

  pkgs[have]
}
#' Show BiocManager install commands for missing OrgDb packages
#'
#' Convenience helper that checks which OrgDb packages are installed
#' (via [echogo_list_orgdb()]) and prints ready-to-copy
#' BiocManager::install() commands for the missing ones.
#'
#' @param pkgs Character vector of OrgDb package names to check.
#'   Defaults to a common set: human, mouse, zebrafish, fly, rat, worm.
#'
#' @return Invisibly returns the data.frame from [echogo_list_orgdb()].
#' @export
echogo_install_orgdb_instructions <- function(
    pkgs = c("org.Hs.eg.db","org.Mm.eg.db","org.Dr.eg.db",
             "org.Dm.eg.db","org.Rn.eg.db","org.Ce.eg.db")) {

  status  <- echogo_list_orgdb(pkgs)
  missing <- status$package[!status$installed]

  if (!length(missing)) {
    cat("All requested OrgDb packages are installed.\n")
    return(invisible(status))
  }

  cat("The following OrgDb packages are missing:\n")
  cat("  ", paste(missing, collapse = ", "), "\n\n", sep = "")

  cat("Install them using BiocManager:\n\n")
  cat("  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n")
  for (pkg in missing) {
    cat("  BiocManager::install('", pkg, "')\n", sep = "")
  }

  invisible(status)
}


# --- R/install_helpers.R -----------------------------------------------------

#' Install GO.db / OrgDb packages into the user library (no prompts)
#'
#' Installs Bioconductor annotation packages to the **user** library to avoid
#' permission issues. If `pkgs = NULL`, it installs any **missing** items from a
#' standard OrgDb set plus GO.db. Returns a logical vector (per package) of
#' whether the package is now loadable.
#'
#' @param pkgs Character vector of package names (e.g., c("GO.db","org.Mm.eg.db")).
#'   If NULL, installs missing ones from the default set:
#'   c("GO.db","org.Hs.eg.db","org.Mm.eg.db","org.Dr.eg.db",
#'     "org.Dm.eg.db","org.Rn.eg.db","org.Ce.eg.db").
#' @param update Logical; update existing packages? Default FALSE for reproducibility.
#' @return Named logical vector: TRUE if package is loadable after installation.
#' @export
echogo_install_orgdb <- function(
    pkgs   = NULL,
    update = FALSE
) {
  # Default set if not provided
  default_set <- c(
    "GO.db",
    "org.Hs.eg.db","org.Mm.eg.db","org.Dr.eg.db",
    "org.Dm.eg.db","org.Rn.eg.db","org.Ce.eg.db"
  )

  if (is.null(pkgs)) {
    # Only install what is missing from the default set
    missing <- default_set[!vapply(default_set, requireNamespace, logical(1), quietly = TRUE)]
    pkgs <- if (length(missing)) missing else character(0)
  }

  # Nothing to do?
  if (!length(pkgs)) {
    now <- vapply(default_set, requireNamespace, logical(1), quietly = TRUE)
    names(now) <- default_set
    return(now)
  }

  # Ensure BiocManager + repos + user lib path
  if (!requireNamespace("BiocManager", quietly = TRUE)) utils::install.packages("BiocManager")
  options(repos = BiocManager::repositories())
  dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

  # Install
  BiocManager::install(pkgs, ask = FALSE, update = isTRUE(update))

  # Report loadability after install
  out <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  names(out) <- pkgs
  return(out)
}

#' Ensure OrgDb packages are available (install if missing)
#'
#' This checks that the requested OrgDb packages are available. If they are
#' missing and `auto_install = TRUE`, it calls [echogo_install_orgdb()] to
#' install them into the user library.
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


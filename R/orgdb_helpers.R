#' List available OrgDb annotation packages
#'
#' Shows which organism annotation databases (OrgDb) are installed and which are missing.
#'
#' @param pkgs Character vector of OrgDb package names to check.
#'   Defaults to a common set: human, mouse, zebrafish, fly, rat, worm.
#' @return A data.frame with columns: package, installed (TRUE/FALSE), and loadable (TRUE/FALSE).
#' @export
#' @examples
#' echogo_list_orgdb()
#' echogo_list_orgdb(c("org.Hs.eg.db", "org.Mm.eg.db"))
echogo_list_orgdb <- function(
    pkgs = c("org.Hs.eg.db","org.Mm.eg.db","org.Dr.eg.db",
             "org.Dm.eg.db","org.Rn.eg.db","org.Ce.eg.db")) {

  res <- data.frame(
    package   = pkgs,
    installed = vapply(pkgs, function(p) requireNamespace(p, quietly = TRUE), logical(1)),
    stringsAsFactors = FALSE
  )
  res$loadable <- vapply(pkgs, function(p) {
    if (!res$installed[res$package == p]) return(FALSE)
    ok <- tryCatch({
      suppressPackageStartupMessages(library(p, character.only = TRUE))
      TRUE
    }, error = function(e) FALSE)
    ok
  }, logical(1))

  rownames(res) <- NULL
  return(res)
}

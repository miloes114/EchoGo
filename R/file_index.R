#' Write a file index for an EchoGO output directory
#'
#' Recursively scans \code{dir} and writes \code{__file_index.csv} there.
#' Returns the tibble (invisibly).
#'
#' @param dir Path to the EchoGO output directory (e.g. "echogo_out_test").
#' @param index_name File name for the index (default "__file_index.csv").
#' @return Invisibly, a tibble with columns rel_path and full_path.
#' @export
echogo_write_index <- function(dir, index_name = "__file_index.csv") {
  stopifnot(is.character(dir), length(dir) == 1)
  dir_abs <- normalizePath(dir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(dir_abs)) {
    dir.create(dir_abs, recursive = TRUE, showWarnings = FALSE)
  }
  rels <- list.files(dir_abs, recursive = TRUE, all.files = TRUE,
                     include.dirs = FALSE, no.. = TRUE)
  idx <- tibble::tibble(
    rel_path  = gsub("\\\\", "/", rels),
    full_path = gsub("\\\\", "/", normalizePath(file.path(dir_abs, rels),
                                                winslash = "/", mustWork = FALSE))
  )
  out <- file.path(dir_abs, index_name)
  readr::write_csv(idx, out)
  message("echogo_write_index(): wrote ", out)
  invisible(idx)
}

#' Read the EchoGO file index (optionally require it to exist)
#'
#' @param dir Path to the EchoGO output directory.
#' @param index_name Index file name (default "__file_index.csv").
#' @param create_if_missing Logical; if TRUE and index is missing, create it.
#'                          Default FALSE to avoid writing from Rmd when knitted.
#' @return A tibble with rel_path and full_path (possibly empty).
#' @export
echogo_read_index <- function(dir, index_name = "__file_index.csv",
                              create_if_missing = FALSE) {
  dir_abs <- normalizePath(dir, winslash = "/", mustWork = FALSE)
  idx_path <- file.path(dir_abs, index_name)
  if (!file.exists(idx_path) && isTRUE(create_if_missing)) {
    echogo_write_index(dir_abs, index_name = index_name)
  }
  if (file.exists(idx_path)) {
    return(suppressWarnings(readr::read_csv(idx_path, show_col_types = FALSE)))
  }
  tibble::tibble(rel_path = character(), full_path = character())
}

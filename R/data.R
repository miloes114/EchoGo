#' Curated taxonomy fallback for EchoGO species metadata
#'
#' A data frame keyed by NCBI taxonomy ID with available taxonomic ranks from
#' superkingdom through genus. EchoGO uses it when online taxonomy lookup is
#' unavailable.
#'
#' @format A data frame with one row per NCBI taxonomy ID and columns including
#'   `ncbi`, `superkingdom`, `kingdom`, `phylum`, `class`, `order`, `family`,
#'   and `genus`.
#' @source Curated from NCBI Taxonomy for the packaged g:Profiler fallback set.
"echogo_taxonomy_fallback"

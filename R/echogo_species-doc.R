#' g:Profiler species fallback (packaged)
#'
#' A cached list of g:Profiler organisms used when no internet connection
#' is available. Generated via \code{echogo_update_species_cache()}.
#'
#' @format A tibble/data.frame with columns:
#' \describe{
#'   \item{organism}{Character. g:Profiler organism ID (e.g., "hsapiens").}
#'   \item{name}{Character. Display/scientific name.}
#'   \item{ncbi}{Integer. NCBI taxonomy ID (if available).}
#'   \item{alias}{Character. Optional aliases (semicolon-separated).}
#' }
#' @source Retrieved from g:Profiler API
#'   (\code{/api/util/organisms_list}).
"echogo_species"

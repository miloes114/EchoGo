#' Freeze demo results for offline vignettes
#'
#' Copies a completed local EchoGO demo run into
#' `inst/extdata/echogo_demo_results/` so that vignettes can be rendered
#' offline and independent of user-specific paths or environments.
#'
#' This helper expects that the user has already executed:
#' `echogo_quickstart(run_demo = TRUE)`, which generates a full demo
#' directory in the location returned by:
#' `getOption("EchoGO.demo_root", rappdirs::user_data_dir("EchoGO"))`.
#'
#' The function performs the following actions:
#' \enumerate{
#'   \item Locates the local demo results directory.
#'   \item Creates the destination folder under `inst/extdata/`.
#'   \item Removes any existing contents (but keeps the folder).
#'   \item Copies all files and subdirectories recursively.
#'   \item Warns if any files could not be copied.
#' }
#'
#' No values are returned; the function is used for package maintenance
#' and is not exported.
#'
#' @section Used in:
#' Package build scripts (`data-raw/`) to ensure reproducible
#' demo datasets for vignettes and examples.
#'
#' @keywords internal
#' @noRd

# data-raw/freeze_demo_results.R
# Freeze a local demo run into inst/extdata/echogo_demo_results for offline vignettes.

if (!requireNamespace("rappdirs", quietly = TRUE)) {
  stop("Please install 'rappdirs' to locate the demo directory.")
}

src <- file.path(
  getOption("EchoGO.demo_root", rappdirs::user_data_dir("EchoGO")),
  "echogo_demo", "results"
)
dst <- file.path("inst", "extdata", "echogo_demo_results")

if (!dir.exists(src)) {
  stop("Run echogo_quickstart(run_demo = TRUE) once to produce demo results at:\n", src)
}

dir.create(dst, recursive = TRUE, showWarnings = FALSE)

# Clean destination (keep folder but remove contents)
unlink(list.files(dst, all.files = TRUE, full.names = TRUE, no.. = TRUE),
       recursive = TRUE, force = TRUE)

ok <- file.copy(list.files(src, full.names = TRUE, no.. = TRUE),
                dst, recursive = TRUE, overwrite = TRUE)

if (!all(ok)) warning("Some demo results could not be copied.")
message("[EchoGO] Frozen demo results at: ",
        normalizePath(dst, winslash = "/"))

#' Build complete offline taxonomy fallback (superkingdom → genus)
#'
#' Generates a full NCBI Taxonomy fallback table for **all species with a
#' valid NCBI Taxonomy ID** in the EchoGO species metadata. The resulting
#' object, `echogo_taxonomy_fallback`, contains taxonomic ranks from
#' *superkingdom* to *genus* and is saved into
#' `data/echogo_taxonomy_fallback.rda`.
#'
#' This script is designed for maintainers and supports **resumable
#' execution** via a progress cache stored in `data-raw/_taxonomy_progress.rds`.
#' This makes it suitable for large multi-day fetches across many species.
#'
#' @section Workflow:
#' \enumerate{
#'   \item Loads the full species table using
#'   `echogo_gprofiler_species()`, preferring fresh retrieval.
#'
#'   \item Extracts all unique NCBI Taxonomy IDs.
#'
#'   \item Applies a polite NCBI query rate based on presence of an
#'   `ENTREZ_KEY` (≈9 req/s with key, ≈3 req/s without).
#'
#'   \item Fetches lineage information for each species using
#'   \pkg{taxize} (`classification()`), storing intermediate results
#'   in a persistent cache that can survive interruptions.
#'
#'   \item Periodically checkpoints the cache (every 100 species).
#'
#'   \item After completion, compiles all cached lineages into a
#'   wide taxonomy table with ranks:
#'   \code{superkingdom, kingdom, phylum, class, order, family, genus}.
#'
#'   \item Saves the final object as
#'   `data/echogo_taxonomy_fallback.rda`, using \pkg{usethis} when
#'   available.
#' }
#'
#' @details
#' This dataset provides complete offline taxonomic resolution for all
#' supported species, enabling EchoGO to operate without internet access,
#' during server outages, or when API limits are reached. Compared to the
#' minimal fallback (which covers only core model species), this script
#' produces a full fallback covering the entire species table.
#'
#' The progress cache is intentionally not removed by default so that the
#' script may be rerun incrementally. Maintainers may uncomment the cleanup
#' line at the end to force rebuilding from scratch.
#'
#' @note
#' This script is not part of the user-facing API. It exists exclusively for
#' package preparation under `inst/scripts/` and should be run manually when
#' species tables or taxonomic data are updated.
#'
#' @seealso
#' \code{\link[taxize]{classification}}
#' \code{\link{echogo_gprofiler_species}}
#' `data-raw/build_taxonomy_fallback_minimal.R` for the lightweight fallback.
#'
#' @keywords internal data-raw taxonomy fallback NCBI resumable
#' @noRd


# inst/scripts/build_taxonomy_fallback.R
# Build a full offline taxonomy fallback (superkingdom..genus) for ALL species
# with an NCBI taxonomy id in the EchoGO species table.
# Output: data/echogo_taxonomy_fallback.rda
# Safe to resume if interrupted (uses a progress cache in data-raw/).

suppressPackageStartupMessages({
  has_taxize <- requireNamespace("taxize", quietly = TRUE)
  has_usethis <- requireNamespace("usethis", quietly = TRUE)
})
if (!has_taxize) stop("Please install 'taxize' first: install.packages('taxize')")

# ---- Load species table (prefer fresh, else fallback) -----------------------
sp <- tryCatch(EchoGO::echogo_gprofiler_species(refresh = TRUE),
               error = function(e) EchoGO::echogo_gprofiler_species(refresh = FALSE))

ncbi_ids <- unique(sp$ncbi[!is.na(sp$ncbi)])
if (!length(ncbi_ids)) stop("No NCBI taxonomy IDs found in species table.")

# ---- Rate limit based on ENTrez key presence --------------------------------
use_key   <- nzchar(Sys.getenv("ENTREZ_KEY", ""))
sleep_sec <- if (use_key) 0.11 else 0.35  # ~9/s with key; ~3/s without
message("[EchoGO] ENTREZ_KEY present: ", use_key, " | rate ~", round(1/sleep_sec, 1), " req/s")
message("[EchoGO] Total NCBI IDs to fetch: ", length(ncbi_ids))

# ---- Resumable progress cache -----------------------------------------------
dir.create("data-raw", showWarnings = FALSE)
progress_file <- file.path("data-raw", "_taxonomy_progress.rds")

lin_cache <- if (file.exists(progress_file)) {
  message("[EchoGO] Resuming from: ", progress_file)
  readRDS(progress_file)
} else {
  new.env(parent = emptyenv())
}

# Helper: store lineage table (name,rank) for a taxid as data.frame
.cache_put <- function(tx, cl_df) assign(as.character(tx), cl_df, envir = lin_cache)
.cache_get <- function(tx)       get(as.character(tx), envir = lin_cache, inherits = FALSE)
.cache_has <- function(tx)       exists(as.character(tx), envir = lin_cache, inherits = FALSE)
.cache_keys <- function()        ls(envir = lin_cache)

# ---- Lineage fetcher (unwrap taxize list result) ----------------------------
ranks <- c("superkingdom","kingdom","phylum","class","order","family","genus")

fetch_lineage_df <- function(tx) {
  cl <- try(taxize::classification(tx, db = "ncbi"), silent = TRUE)
  if (inherits(cl, "try-error")) return(NULL)
  if (is.list(cl) && length(cl) > 0) cl <- cl[[1]]  # unwrap single-result list
  if (!is.data.frame(cl) || !all(c("name","rank") %in% names(cl))) return(NULL)
  cl[, c("name","rank")]
}

# ---- Iterate with periodic checkpointing ------------------------------------
todo <- setdiff(ncbi_ids, as.integer(.cache_keys()))
n_total <- length(ncbi_ids)
n_todo  <- length(todo)
message("[EchoGO] Already cached: ", n_total - n_todo, " | To fetch: ", n_todo)

if (n_todo > 0) {
  for (i in seq_along(todo)) {
    tx <- todo[i]
    if (i %% 50 == 1) {
      message(sprintf("[EchoGO] %d / %d (%.1f%%)", i, n_todo, 100 * i / max(1, n_todo)))
    }
    cl_df <- fetch_lineage_df(tx)
    if (!is.null(cl_df)) .cache_put(tx, cl_df)
    if (sleep_sec > 0) Sys.sleep(sleep_sec)
    # checkpoint every 100
    if (i %% 100 == 0) saveRDS(lin_cache, progress_file)
  }
  # final checkpoint
  saveRDS(lin_cache, progress_file)
}

# ---- Compile cache into the final wide table --------------------------------
keys <- as.integer(.cache_keys())
if (!length(keys)) stop("No taxonomy rows fetched; check network/ENTREZ_KEY and try again.")

as_row <- function(tx) {
  cl_df <- .cache_get(tx)
  out <- rep(NA_character_, length(ranks)); names(out) <- ranks
  if (is.data.frame(cl_df) && nrow(cl_df)) {
    m <- match(ranks, cl_df$rank); hit <- which(!is.na(m))
    if (length(hit)) out[hit] <- cl_df$name[m[hit]]
  }
  data.frame(ncbi = as.integer(tx), t(out), check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)
}

rows <- lapply(keys, as_row)
fallback <- do.call(rbind, rows)

# Deduplicate just in case (keep first)
fallback <- fallback[!duplicated(fallback$ncbi), , drop = FALSE]

# ---- Save packaged dataset ---------------------------------------------------
dir.create("data", showWarnings = FALSE)
echogo_taxonomy_fallback <- fallback
if (has_usethis) {
  usethis::use_data(echogo_taxonomy_fallback, overwrite = TRUE, compress = "xz")
} else {
  save(echogo_taxonomy_fallback, file = "data/echogo_taxonomy_fallback.rda", compress = "xz")
}
message("[EchoGO] Saved data/echogo_taxonomy_fallback.rda with ", nrow(echogo_taxonomy_fallback), " rows.")

# ---- Clean up progress (optional; keep if you want resume for future updates)
# unlink(progress_file)

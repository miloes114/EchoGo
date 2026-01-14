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

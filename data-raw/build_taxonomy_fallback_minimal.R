#' Build minimal offline taxonomy fallback (superkingdom → genus)
#'
#' Constructs a compact NCBI Taxonomy fallback table for a curated set of
#' reference organisms used throughout EchoGO. The resulting object,
#' `echogo_taxonomy_fallback`, contains taxonomy ranks from
#' *superkingdom* to *genus* and is saved to `data/echogo_taxonomy_fallback.rda`
#' for use when online taxonomy queries are unavailable.
#'
#' The script performs the following steps:
#' \enumerate{
#'   \item Determines the set of core organism codes using
#'   `inst/extdata/species_tags.yaml` when available, falling back to a
#'   predefined list of key model species.
#'
#'   \item Retrieves the corresponding NCBI Taxonomy IDs via
#'   `echogo_gprofiler_species()`.
#'
#'   \item Queries NCBI (via \pkg{taxize}) for each species’
#'   lineage, extracting only the ranks:
#'   `superkingdom`, `kingdom`, `phylum`, `class`, `order`,
#'   `family`, and `genus`.
#'
#'   \item Applies a polite query rate depending on the presence of
#'   an `ENTREZ_KEY`.
#'
#'   \item Deduplicates results and writes the final object to
#'   `data/echogo_taxonomy_fallback.rda` using \pkg{usethis} when available.
#' }
#'
#' This fallback is used internally whenever online taxonomy queries fail
#' (e.g., no internet, API rate limits, missing NCBI credentials), ensuring
#' that EchoGO retains lightweight, deterministic taxonomy resolution for
#' common reference species.
#'
#' @details
#' Only minimal taxonomy ranks are stored to keep package size small while
#' providing sufficient resolution for taxonomy-aware features (e.g., species
#' panel summaries, rarefaction reporting, diagnostic checks).
#'
#' @note
#' This script is intended for package maintainers and is not part of the
#' public API. It is placed under `data-raw/` and executed manually when
#' taxonomy fallback data need to be refreshed.
#'
#' @seealso
#' \code{\link[taxize]{classification}} for taxonomy retrieval,
#' \code{\link{echogo_gprofiler_species}} for species metadata.
#'
#' @keywords internal data-raw taxonomy fallback
#' @noRd


# data-raw/build_taxonomy_fallback_minimal.R
# Build a compact offline taxonomy fallback (superkingdom..genus) for a curated set.
# Writes: data/echogo_taxonomy_fallback.rda

suppressPackageStartupMessages({
  has_taxize <- requireNamespace("taxize", quietly = TRUE)
  has_yaml   <- requireNamespace("yaml",   quietly = TRUE)
})
if (!has_taxize) stop("Please install 'taxize' first: install.packages('taxize')")

# --- Choose the core organism IDs -------------------------------------------
get_core_ids <- function() {
  # Prefer tags from YAML if present
  yaml_path <- system.file("extdata", "species_tags.yaml", package = "EchoGO")
  if (has_yaml && file.exists(yaml_path)) {
    tags <- tryCatch(yaml::read_yaml(yaml_path), error = function(e) NULL)
    if (is.list(tags) && length(tags)) {
      panels <- intersect(names(tags), c("AnimalModels","PlantModels","FungalModels","CoreModels"))
      ids <- unique(unlist(unname(tags[panels])))
      ids <- ids[!is.na(ids) & nzchar(ids)]
      return(sort(unique(ids)))
    }
  }
  # Fallback hard-coded core set (edit if needed)
  c(
    "hsapiens","mmusculus","rnorvegicus","drerio","dmelanogaster","celegans",
    "ggallus","xtropicalis",
    "athaliana","osativa","zmays","slycopersicum",  # tomato id often 'slycopersicum'
    "scerevisiae","spombe"
  )
}

core_ids <- get_core_ids()

# --- Load species table & map to NCBI IDs -----------------------------------
sp <- tryCatch(EchoGO::echogo_gprofiler_species(refresh = TRUE),
               error = function(e) EchoGO::echogo_gprofiler_species(refresh = FALSE))

seed <- sp[sp$organism %in% core_ids & !is.na(sp$ncbi), c("organism","ncbi")]
missing_ids <- setdiff(core_ids, seed$organism)
if (length(missing_ids)) {
  message("[EchoGO] Note: these core IDs are not present (or lack NCBI) in current species table: ",
          paste(missing_ids, collapse = ", "))
}
if (!nrow(seed)) stop("No core species with valid NCBI IDs; cannot build fallback.")

# --- Configure polite rate (key optional) -----------------------------------
use_key   <- nzchar(Sys.getenv("ENTREZ_KEY", ""))
sleep_sec <- if (use_key) 0.11 else 0.35  # ~9/s with key; ~3/s without
message("[EchoGO] ENTREZ_KEY present: ", use_key, "  |  rate ~", round(1/sleep_sec, 1), " req/s")

# --- Fetch lineages ----------------------------------------------------------
ranks <- c("superkingdom","kingdom","phylum","class","order","family","genus")

get_lineage <- function(tx) {
  cl <- try(taxize::classification(tx, db = "ncbi"), silent = TRUE)
  if (inherits(cl, "try-error")) return(NULL)
  if (is.list(cl) && length(cl) > 0) cl <- cl[[1]]  # unwrap list
  if (!is.data.frame(cl) || !all(c("name","rank") %in% names(cl))) return(NULL)

  out <- setNames(as.list(rep(NA_character_, length(ranks))), ranks)
  m <- match(ranks, cl$rank); hit <- which(!is.na(m))
  if (length(hit)) out[hit] <- as.list(cl$name[m[hit]])
  data.frame(ncbi = as.integer(tx), out, check.names = FALSE, stringsAsFactors = FALSE)
}

message("[EchoGO] Fetching taxonomy for ", nrow(seed), " core species...")
res <- vector("list", nrow(seed))
for (i in seq_len(nrow(seed))) {
  if (i %% 20 == 1) message(sprintf("  - %d / %d", i, nrow(seed)))
  res[[i]] <- get_lineage(seed$ncbi[i])
  if (sleep_sec > 0) Sys.sleep(sleep_sec)
}

fallback <- do.call(rbind, Filter(Negate(is.null), res))
if (is.null(fallback) || !nrow(fallback)) {
  stop("No taxonomy rows fetched; check network/ENTREZ_KEY and try again.")
}
fallback <- fallback[!duplicated(fallback$ncbi), , drop = FALSE]

# --- Save to data/ as echogo_taxonomy_fallback.rda ---------------------------
dir.create("data", showWarnings = FALSE)
echogo_taxonomy_fallback <- fallback

if (requireNamespace("usethis", quietly = TRUE)) {
  usethis::use_data(echogo_taxonomy_fallback, overwrite = TRUE, compress = "xz")
} else {
  save(echogo_taxonomy_fallback, file = "data/echogo_taxonomy_fallback.rda", compress = "xz")
}
message("[EchoGO] Saved data/echogo_taxonomy_fallback.rda with ",
        nrow(echogo_taxonomy_fallback), " rows.")

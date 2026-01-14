#' @title Resolve reference-based RNA-seq inputs for EchoGO
#' @description
#' Detects the standard reference-based input layout:
#'   - allcounts_table.txt            (background counts)
#'   - dge_*.csv                      (DESeq2 results for one contrast)
#'   - *.GOseq.enriched.tsv           (GOseq over-represented categories)
#'   - Trout_eggNOG_for_EchoGO.tsv    (Trinotate-like eggNOG table)
#'   - config.yml                     (species, orgdb, report title)
#'
#' and returns a list of normalized paths for downstream functions.
#'
#' @param input_dir Input folder (typically ".../dge_CONTRAST/input").
#' @return A named list with elements:
#'   \code{goseq_file}, \code{trinotate_file}, \code{de_file},
#'   \code{count_matrix_file}, \code{config_file}, \code{contrast}.
#' @export
echogo_resolve_reference_inputs <- function(input_dir) {
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  files <- list.files(input_dir, full.names = TRUE)

  # GOseq enriched (required)
  goseq_file <- grep("GOseq\\.enriched\\.", files, value = TRUE)
  if (!length(goseq_file)) {
    stop("Could not find a GOseq enriched file (pattern '*GOseq.enriched.*') in: ", input_dir)
  }
  if (length(goseq_file) > 1L) {
    message("Multiple GOseq enriched files found; using first: ", basename(goseq_file[1]))
  }
  goseq_file <- goseq_file[1]

  # Trinotate-like annotation (required)
  trinotate_file <- grep("Trout_eggNOG_for_EchoGO|Trinotate\\.", files, value = TRUE)
  if (!length(trinotate_file)) {
    stop("Could not find Trinotate-like annotation (Trout_eggNOG_for_EchoGO*.tsv or Trinotate.*) in: ", input_dir)
  }
  if (length(trinotate_file) > 1L) {
    message("Multiple Trinotate-like files found; using first: ", basename(trinotate_file[1]))
  }
  trinotate_file <- trinotate_file[1]

  # DESeq2 table (optional but useful for logging / fold estimation)
  de_file <- grep("^.+dge_.*\\.csv$", files, value = TRUE)
  if (!length(de_file)) {
    de_file <- grep("dge_.*\\.csv$", files, value = TRUE)
  }
  if (!length(de_file)) {
    message("No DESeq2 CSV detected (dge_*.csv); 'de_file' will be NA.")
    de_file <- NA_character_
  } else {
    if (length(de_file) > 1L) {
      message("Multiple dge_*.csv files found; using first: ", basename(de_file[1]))
    }
    de_file <- de_file[1]
  }

  # Count matrix (optional; used only if you want EchoGO to recompute foldEnrichment)
  count_matrix_file <- grep("allcounts_table|gene.counts.matrix", files, value = TRUE)
  if (!length(count_matrix_file)) {
    message("No count matrix found (allcounts_table.txt / gene.counts.matrix.*).",
            " Fold enrichment in GOseq will be used as-is or rederived from numDEInCat/numInCat.")
    count_matrix_file <- NA_character_
  } else {
    if (length(count_matrix_file) > 1L) {
      message("Multiple count matrices found; using first: ", basename(count_matrix_file[1]))
    }
    count_matrix_file <- count_matrix_file[1]
  }

  # Config file (optional but recommended)
  config_file <- grep("config\\.yml$", files, value = TRUE)
  if (!length(config_file)) {
    config_file <- NA_character_
  } else {
    config_file <- config_file[1]
  }

  # Infer contrast from GOseq enriched name (e.g. "dge_RRf_BBf")
  contrast <- sub("\\.GOseq\\.enriched.*$", "", basename(goseq_file))

  list(
    input_dir         = input_dir,
    goseq_file        = goseq_file,
    trinotate_file    = trinotate_file,
    de_file           = de_file,
    count_matrix_file = count_matrix_file,
    config_file       = config_file,
    contrast          = contrast
  )
}

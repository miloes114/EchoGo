#' Build a compact EchoGO demo dataset from a completed project
#'
#' This helper constructs the demo input files shipped with EchoGO
#' (counts, DE results, GOseq enrichment and a synthetic Trinotate table)
#' starting from a real project that has already been processed by
#' \code{run_full_echogo()}.
#'
#' The function:
#' \itemize{
#'   \item Reads the consensus workbook and GOseq enrichment table.
#'   \item Extracts real orthology-based symbols from g:Profiler results.
#'   \item Synthesizes a minimal \code{Trinotate_demo.tsv} mapping
#'         TRINITY transcript IDs to symbols (Metazoa-tagged).
#'   \item Rewrites GOseq \code{gene_ids} so they only contain mapped
#'         transcripts (critical for a self-consistent demo).
#'   \item Filters DE results and counts to the same mapped transcript set.
#' }
#'
#' @param src_results Path to a completed EchoGO results directory. This is
#'   expected to contain at least the subdirectories \file{goseq/} and
#'   \file{consensus/} with the standard EchoGO file names.
#' @param src_input Path to the original input directory used for the full
#'   run (typically the \code{input/} created by \code{echogo_scaffold()}),
#'   containing the DE results and count matrix.
#' @param demo_dir Output directory where the distilled demo files will be
#'   written. The directory is created if it does not exist.
#' @param N_GOSEQ Maximum number of GOseq terms to include in the demo
#'   enrichment table.
#' @param N_DE Maximum number of differential expression rows to keep in
#'   \file{DE_results_demo.tsv}.
#' @param N_COUNTS Maximum number of rows to keep in the demo count matrix.
#' @param N_SAMPLES Maximum number of sample columns to keep in the demo
#'   count matrix.
#' @param N_MAP Maximum number of TRINITY transcript IDs to map to gene
#'   symbols in the synthetic Trinotate table.
#'
#' @return Invisibly returns \code{normalizePath(demo_dir)}. The function is
#'   called for its side-effects of writing the following files:
#'   \itemize{
#'     \item \file{GOseq_enrichment_demo.csv}
#'     \item \file{Trinotate_demo.tsv}
#'     \item \file{DE_results_demo.tsv} (if a DE source is found)
#'     \item \file{counts_demo.tsv} (if a counts source is found)
#'     \item \file{README_demo.txt}
#'   }
#'
#' @details
#' This function is mainly intended for maintainers to regenerate the
#' frozen demo data under \file{inst/extdata/echogo_demo} from a larger
#' “real” project. It assumes the standard EchoGO naming conventions:
#' \itemize{
#'   \item \file{consensus/consensus_enrichment_results_with_and_without_bg.xlsx}
#'   \item \file{goseq/GOseq_enrichment_full_annotated.csv}
#'   \item An input count matrix (\code{*count*matrix*} or \code{*counts*.tsv})
#'   \item One or more DE result tables matching \code{DE_*.tsv} or
#'         \code{*DE_results*subset*.tsv}
#' }
#'
#' If some of these inputs are not found, the function throws an error for
#' the core consensus/GOseq files and silently skips optional outputs
#' (DE/counts) with a warning.
#'
#' @examples
#' \dontrun{
#' echogo_build_demo(
#'   src_results = "echogo_project/results",
#'   src_input   = "echogo_project/input",
#'   demo_dir    = "inst/extdata/echogo_demo"
#' )
#' }
#'
#' @keywords internal
echogo_build_demo <- function(
    src_results = "echogo_project/results",
    src_input   = "echogo_project/input",
    demo_dir    = file.path("inst", "extdata", "echogo_demo"),
    N_GOSEQ     = 80,
    N_DE        = 400,
    N_COUNTS    = 250,
    N_SAMPLES   = 8,
    N_MAP       = 400
) {

  # -------------------- HELPERS ---------------------------------------------

  #' Internal infix: fallback for NULL / empty values
  #'
  #' Returns \code{b} when \code{a} is \code{NULL}, length zero, or entirely
  #' \code{NA}; otherwise returns \code{a}.
  #'
  #' @param a Primary value.
  #' @param b Fallback value.
  #'
  #' @noRd
  `%||%` <- function(a, b) {
    if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
  }

  #' Pick the first matching file from a set of glob patterns
  #'
  #' @param ... One or more file paths or glob patterns (passed to
  #'   \code{Sys.glob()}).
  #' @return A single character string with the first existing match, or
  #'   \code{NA_character_} if no match is found.
  #'
  #' @noRd
  pick_one <- function(...) {
    pats <- unlist(list(...))
    hit  <- unlist(lapply(pats, Sys.glob))
    if (length(hit)) hit[[1]] else NA_character_
  }

  #' Split a comma-separated vector into a clean character vector
  #'
  #' @param x Character vector with comma-separated entries.
  #' @return A character vector with trimmed, non-empty elements.
  #'
  #' @noRd
  split_list <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(character(0))
    out <- unlist(strsplit(paste(x, collapse = ","), ",\\s*"))
    out <- trimws(out)
    out[nzchar(out)]
  }

  #' Repair malformed count columns from DE/Count matrices
  #'
  #' Some external pipelines produce count matrix columns that contain
  #' tab-separated values inside a single field. This helper splits such
  #' columns into two numeric columns with suffixes \code{"_1"} and
  #' \code{"_2"} and parses simple numeric columns as doubles.
  #'
  #' @param df A data frame or tibble with count data.
  #' @return A tibble with only numeric/sample columns and repaired names.
  #'
  #' @noRd
  repair_count_columns <- function(df) {
    fix_one <- function(v) {
      if (!is.character(v)) return(v)
      if (!any(grepl("\\t", v, fixed = TRUE))) {
        return(suppressWarnings(readr::parse_double(v)))
      }
      m <- stringr::str_split_fixed(v, "\\t", 2)
      tibble::tibble(
        V1 = suppressWarnings(readr::parse_double(m[, 1])),
        V2 = suppressWarnings(readr::parse_double(m[, 2]))
      )
    }

    out <- list()
    for (nm in names(df)) {
      col   <- df[[nm]]
      fixed <- fix_one(col)
      if (is.data.frame(fixed)) {
        out[[paste0(nm, "_1")]] <- fixed$V1
        out[[paste0(nm, "_2")]] <- fixed$V2
      } else {
        out[[nm]] <- fixed
      }
    }
    tibble::as_tibble(out, .name_repair = "unique")
  }

  # -------------------- ENSURE OUTPUT DIR -----------------------------------
  dir.create(demo_dir, recursive = TRUE, showWarnings = FALSE)

  # -------------------- LOCATE INPUTS FROM THE FULL RUN ---------------------
  cons_xlsx <- pick_one(file.path(src_results, "consensus",
                                  "consensus_enrichment_results_with_and_without_bg.xlsx"))
  goseq_csv <- pick_one(file.path(src_results, "goseq",
                                  "GOseq_enrichment_full_annotated.csv"))

  if (is.na(cons_xlsx) || is.na(goseq_csv)) {
    stop(
      "Could not find consensus and/or GOseq outputs under ",
      normalizePath(src_results, winslash = "/")
    )
  }

  cts_in <- pick_one(
    file.path(src_input, "*count*matrix*"),
    file.path(src_input, "*counts*.tsv"),
    file.path(src_input, "*counts*.csv")
  )
  de_in <- pick_one(
    file.path(src_input, "DE_*.tsv"),
    file.path(src_input, "*DE_results*subset*.tsv")
  )
  goseq_in <- pick_one(
    file.path(src_input, "*DE_results*.DE.subset.GOseq.enriched"),
    file.path(src_input, "*GOseq*enrich*")
  )

  # -------------------- READ FULL RESULTS -----------------------------------
  cons <- readxl::read_xlsx(cons_xlsx)

  # Real symbols from g:Profiler columns in consensus
  gp_syms <- unique(c(
    split_list(cons$genes_gprofiler_bg   %||% character(0)),
    split_list(cons$genes_gprofiler_nobg %||% character(0))
  ))
  gp_syms <- gp_syms[nzchar(gp_syms)]

  # GOseq annotated (prefer; may or may not include gene_ids)
  goseq_full <- readr::read_csv(goseq_csv, show_col_types = FALSE)

  goseq_pick <- goseq_full %>%
    dplyr::mutate(.rank = rank(over_represented_FDR,
                               ties.method = "first", na.last = "keep")) %>%
    dplyr::arrange(.rank) %>%
    dplyr::slice_head(n = N_GOSEQ) %>%
    dplyr::select(
      category,
      term,
      ontology,
      numDEInCat,
      numInCat,
      over_represented_FDR,
      gene_ids = dplyr::any_of("gene_ids")
    )

  # If annotated file lacks gene_ids, fall back to raw enriched input
  if (!"gene_ids" %in% names(goseq_pick) || all(is.na(goseq_pick$gene_ids))) {
    if (!is.na(goseq_in)) {
      goseq_raw <- suppressMessages(readr::read_delim(
        goseq_in,
        delim = if (grepl("\\.csv$", goseq_in, TRUE)) "," else "\t",
        show_col_types = FALSE
      ))
      key_cols <- c(
        "category", "term", "ontology",
        "numDEInCat", "numInCat", "over_represented_FDR", "gene_ids"
      )
      common <- intersect(key_cols, names(goseq_raw))
      goseq_pick <- goseq_raw %>%
        dplyr::arrange(over_represented_FDR) %>%
        dplyr::slice_head(n = N_GOSEQ) %>%
        dplyr::select(dplyr::all_of(common))
    } else {
      warning("No gene_ids available for GOseq demo; mapping may be limited.")
    }
  }

  # -------------------- SYNTHETIC TRINOTATE & MAPPING -----------------------
  demo_transcripts <- split_list(goseq_pick$gene_ids)
  demo_transcripts <- unique(demo_transcripts)
  if (!length(demo_transcripts)) {
    stop("No TRINITY transcripts in GOseq gene_ids.")
  }

  max_map <- min(length(demo_transcripts),
                 max(length(gp_syms), 1L), N_MAP)
  if (max_map < 5) {
    warning("Very few symbols found in consensus; g:Profiler may be sparse in the demo.")
  }

  map_tbl <- tibble::tibble(
    transcript_id = demo_transcripts[seq_len(max_map)],
    gene_symbol   = (gp_syms %||% "TP53")[
      (seq_len(max_map) - 1) %% length(gp_syms %||% "TP53") + 1L
    ]
  )
  mapped_ids <- map_tbl$transcript_id

  # Trinotate demo (meets loader’s heuristic)
  tri_demo <- map_tbl %>%
    dplyr::transmute(
      gene_id              = transcript_id,
      transcript_id        = transcript_id,
      sprot_Top_BLASTX_hit = paste0(gene_symbol, "^Metazoa"),
      EggNM.Preferred_name = gene_symbol,
      EggNM.max_annot_lvl  = 33208L
    )
  readr::write_tsv(tri_demo, file.path(demo_dir, "Trinotate_demo.tsv"))

  # -------------------- REWRITE GOseq gene_ids TO MAPPED SET -----------------
  if ("gene_ids" %in% names(goseq_pick)) {
    goseq_pick <- goseq_pick %>%
      dplyr::mutate(
        gene_ids = vapply(gene_ids, function(g) {
          ids <- split_list(g)
          ids <- ids[ids %in% mapped_ids]
          paste(ids, collapse = ", ")
        }, character(1))
      ) %>%
      dplyr::filter(!is.na(gene_ids) & nzchar(gene_ids))
  }
  goseq_out <- goseq_pick %>%
    dplyr::select(
      category,
      term,
      ontology,
      numDEInCat,
      numInCat,
      over_represented_FDR,
      gene_ids = dplyr::any_of("gene_ids")
    )
  readr::write_csv(goseq_out, file.path(demo_dir, "GOseq_enrichment_demo.csv"))

  # -------------------- DE (filtered to mapped transcripts) ------------------
  if (!is.na(de_in)) {
    de_full <- suppressMessages(readr::read_delim(
      de_in,
      delim = if (grepl("\\.csv$", de_in, TRUE)) "," else "\t",
      show_col_types = FALSE
    ))
    if (!"gene_id" %in% names(de_full)) {
      cand <- intersect(c("id", "gene", "GeneID", "X1"), names(de_full))
      if (length(cand)) {
        names(de_full)[match(cand[1], names(de_full))] <- "gene_id"
      }
    }
    de_demo <- de_full %>%
      dplyr::filter(gene_id %in% mapped_ids) %>%
      dplyr::slice_head(n = N_DE)
    readr::write_tsv(de_demo, file.path(demo_dir, "DE_results_demo.tsv"))
  }

  # -------------------- COUNTS (filtered to mapped transcripts) --------------
  if (!is.na(cts_in)) {
    cts_full <- suppressMessages(readr::read_delim(
      cts_in,
      delim = if (grepl("\\.csv$", cts_in, TRUE)) "," else "\t",
      show_col_types = FALSE
    ))
    if (!"gene_id" %in% names(cts_full)) {
      names(cts_full)[1] <- "gene_id"
    }

    cts_demo <- cts_full %>%
      dplyr::filter(gene_id %in% mapped_ids) %>%
      dplyr::slice_head(n = N_COUNTS)
    cts_demo <- repair_count_columns(cts_demo)
    sample_cols <- setdiff(names(cts_demo), "gene_id")
    keep_samps  <- head(sample_cols, min(N_SAMPLES, length(sample_cols)))
    cts_demo    <- cts_demo[, c("gene_id", keep_samps), drop = FALSE]

    readr::write_tsv(cts_demo, file.path(demo_dir, "counts_demo.tsv"))
  }

  # -------------------- README ----------------------------------------------
  readme <- c(
    "EchoGO demo dataset (distilled from a completed full run)",
    "---------------------------------------------------------",
    "Files:",
    " - GOseq_enrichment_demo.csv : GOseq subset with gene_ids RESTRICTED to mapped transcripts",
    " - Trinotate_demo.tsv        : TRINITY -> SYMBOL mapping (synthetic, Metazoa-tagged)",
    " - DE_results_demo.tsv       : DE subset for mapped transcripts (if DE source found)",
    " - counts_demo.tsv           : Counts subset for mapped transcripts (if counts source found)",
    "",
    "Run the demo:",
    " demo_dir <- system.file('extdata','echogo_demo', package='EchoGO')",
    " echogo_quickstart(run_demo = TRUE)"
  )
  writeLines(readme, file.path(demo_dir, "README_demo.txt"))

  demo_dir_norm <- normalizePath(demo_dir, winslash = "/")
  message("[EchoGO] Demo distilled to: ", demo_dir_norm)
  invisible(demo_dir_norm)
}

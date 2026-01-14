# Build a small, self-consistent demo from a completed full run (HOL03_test).
# Writes/overwrites: E:/.../EchoGO/inst/extdata/echogo_demo/
#
# Key points:
# - Pull real symbols from your consensus (g:Profiler genes)
# - Synthesize Trinotate_demo.tsv mapping TRINITY -> SYMBOL (Metazoa-tagged)
# - REWRITE GOseq gene_ids to only mapped transcripts (critical!)
# - Filter DE and counts to same mapped set

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

`%||%` <- function(a, b) { if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a }

# -------------------- CONFIG -------------------------------------------------
# Your completed run folder (assumes outputs are in HOL03_test/results/)
src_project <- "E:/Gigascience submission/EchoGo_Gigascience_submission/EchoGO/HOL03_test"
src_results <- file.path(src_project, "results")  # has goseq/, gprofiler/, consensus/
src_input   <- file.path(src_project, "input")    # original inputs

# Replace the existing shipped demo here:
demo_dir <- "E:/Gigascience submission/EchoGo_Gigascience_submission/EchoGO/inst/extdata/echogo_demo"

# overwrite the existing demo files:
WIPE_DEMO_DIR <- TRUE

N_GOSEQ   <- 80     # GO terms to keep
N_DE      <- 400    # DE rows
N_COUNTS  <- 250    # count rows
N_SAMPLES <- 8      # sample columns
N_MAP     <- 400    # max transcript->symbol mappings to synthesize

# -------------------- HELPERS ------------------------------------------------
detect_delim <- function(path, n = 5) {
  x <- readLines(path, n = n, warn = FALSE)
  x <- x[nzchar(x)]
  if (!length(x)) return("\t")
  seps <- c("\t", ",", ";", "|")
  counts <- vapply(seps, function(s) sum(stringr::str_count(x, stringr::fixed(s))), numeric(1))
  seps[[which.max(counts)]] %||% "\t"
}

read_table_auto <- function(path) {
  delim <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else detect_delim(path)
  suppressMessages(readr::read_delim(
    file = path,
    delim = delim,
    show_col_types = FALSE,
    progress = FALSE,
    trim_ws = TRUE,
    comment = "#"
  ))
}

ensure_gene_id <- function(df) {
  # make names safe + unique
  names(df) <- make.names(names(df), unique = TRUE)

  # already good
  if ("gene_id" %in% names(df)) return(df)

  # common alternatives (case-insensitive-ish because of make.names)
  candidates <- c(
    "gene", "Gene", "GeneID", "geneid", "id", "ID",
    "Row.names", "row.names", "rownames", "transcript_id", "Transcript_ID",
    "X1", "X", "...1", "V1"
  )

  hit <- intersect(candidates, names(df))
  if (length(hit) >= 1) {
    names(df)[match(hit[1], names(df))] <- "gene_id"
    return(df)
  }

  # last resort: assume first column is the feature ID
  names(df)[1] <- "gene_id"
  df
}

pick_one <- function(...) {
  pats <- unlist(list(...))
  hit  <- unlist(lapply(pats, Sys.glob))
  if (length(hit)) hit[[1]] else NA_character_
}

split_list <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(character(0))
  out <- unlist(strsplit(paste(x, collapse=","), ",\\s*"))
  out <- trimws(out)
  out[nzchar(out)]
}

repair_count_columns <- function(df) {
  fix_one <- function(v) {
    if (!is.character(v)) return(v)
    if (!any(grepl("\t", v, fixed = TRUE))) return(suppressWarnings(readr::parse_double(v)))
    m <- str_split_fixed(v, "\t", 2)
    tibble(
      V1 = suppressWarnings(readr::parse_double(m[,1])),
      V2 = suppressWarnings(readr::parse_double(m[,2]))
    )
  }
  out <- list()
  for (nm in names(df)) {
    col <- df[[nm]]
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

# -------------------- PREP DEMO DIR ------------------------------------------
dir.create(demo_dir, recursive = TRUE, showWarnings = FALSE)

if (isTRUE(WIPE_DEMO_DIR)) {
  old <- list.files(demo_dir, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  if (length(old)) unlink(old, recursive = TRUE, force = TRUE)
}

# -------------------- LOCATE OUTPUTS FROM THE FULL RUN -----------------------
cons_xlsx <- pick_one(file.path(src_results, "consensus", "consensus_enrichment_results_with_and_without_bg.xlsx"))
goseq_csv <- pick_one(file.path(src_results, "goseq", "GOseq_enrichment_full_annotated.csv"))

if (is.na(cons_xlsx) || is.na(goseq_csv)) {
  stop(
    "Could not find consensus and/or goseq outputs under:\n  ",
    normalizePath(src_results, winslash = "/"),
    "\nExpected:\n  consensus/consensus_enrichment_results_with_and_without_bg.xlsx\n  goseq/GOseq_enrichment_full_annotated.csv"
  )
}

# -------------------- LOCATE RAW INPUTS (HOL03_test/input) -------------------
# Your files are extensionless; we treat them as TSV.
cts_in <- pick_one(
  file.path(src_input, "*count_matrix*"),
  file.path(src_input, "*count*matrix*"),
  file.path(src_input, "*counts*.tsv"),
  file.path(src_input, "*counts*.csv")
)

de_in <- pick_one(
  file.path(src_input, "*DE_results*DE.subset*"),
  file.path(src_input, "*DE_results*subset*"),
  file.path(src_input, "DE_*.tsv")
)

goseq_in <- pick_one(
  file.path(src_input, "*DE.subset.GOseq.enriched*"),
  file.path(src_input, "*GOseq*enrich*")
)

# -------------------- READ FULL RESULTS -------------------------------------
cons <- readxl::read_xlsx(cons_xlsx)

# Real symbols from g:Profiler columns in consensus
gp_syms <- unique(c(
  split_list(cons$genes_gprofiler_bg %||% character(0)),
  split_list(cons$genes_gprofiler_nobg %||% character(0))
))
gp_syms <- gp_syms[nzchar(gp_syms)]

goseq_full <- readr::read_csv(goseq_csv, show_col_types = FALSE)

goseq_pick <- goseq_full %>%
  mutate(.rank = rank(over_represented_FDR, ties.method = "first", na.last = "keep")) %>%
  arrange(.rank) %>%
  slice_head(n = N_GOSEQ) %>%
  select(
    category, term, ontology, numDEInCat, numInCat, over_represented_FDR,
    gene_ids = any_of("gene_ids")
  )

# If annotated file lacks gene_ids, fall back to raw GOseq enriched input
if (!"gene_ids" %in% names(goseq_pick) || all(is.na(goseq_pick$gene_ids))) {
  if (!is.na(goseq_in)) {
    goseq_raw <- suppressMessages(readr::read_delim(
      goseq_in,
      delim = if (grepl("\\.csv$", goseq_in, ignore.case = TRUE)) "," else "\t",
      show_col_types = FALSE
    ))
    key_cols  <- c("category","term","ontology","numDEInCat","numInCat","over_represented_FDR","gene_ids")
    common    <- intersect(key_cols, names(goseq_raw))
    goseq_pick <- goseq_raw %>%
      arrange(over_represented_FDR) %>%
      slice_head(n = N_GOSEQ) %>%
      select(all_of(common))
  } else {
    warning("No gene_ids available for GOseq demo; mapping may be limited.")
  }
}

# -------------------- SYNTHETIC TRINOTATE & MAPPING -------------------------
demo_transcripts <- split_list(goseq_pick$gene_ids)
demo_transcripts <- unique(demo_transcripts)
if (!length(demo_transcripts)) stop("No TRINITY transcripts found in GOseq gene_ids.")

max_map <- min(length(demo_transcripts), max(length(gp_syms), 1L), N_MAP)
if (max_map < 5) warning("Very few symbols found in consensus; g:Profiler may be sparse in the demo.")

fallback_syms <- gp_syms %||% "TP53"

map_tbl <- tibble(
  transcript_id = demo_transcripts[seq_len(max_map)],
  gene_symbol   = fallback_syms[(seq_len(max_map) - 1) %% length(fallback_syms) + 1L]
)
mapped_ids <- map_tbl$transcript_id

# Minimal Trinotate-like demo mapping (meets EchoGO loader heuristics)
tri_demo <- map_tbl %>%
  transmute(
    gene_id               = transcript_id,
    transcript_id         = transcript_id,
    sprot_Top_BLASTX_hit  = paste0(gene_symbol, "^Metazoa"),
    EggNM.Preferred_name  = gene_symbol,
    EggNM.max_annot_lvl   = 33208L
  )
readr::write_tsv(tri_demo, file.path(demo_dir, "Trinotate_demo.tsv"))

# -------------------- REWRITE GOseq gene_ids TO MAPPED SET -------------------
if ("gene_ids" %in% names(goseq_pick)) {
  goseq_pick <- goseq_pick %>%
    mutate(
      gene_ids = vapply(gene_ids, function(g) {
        ids <- split_list(g)
        ids <- ids[ids %in% mapped_ids]
        paste(ids, collapse = ", ")
      }, character(1))
    ) %>%
    filter(!is.na(gene_ids) & nzchar(gene_ids))
}

goseq_out <- goseq_pick %>%
  select(category, term, ontology, numDEInCat, numInCat, over_represented_FDR, gene_ids = any_of("gene_ids"))
readr::write_tsv(goseq_out, file.path(demo_dir, "GOseq_enrichment_demo.tsv"))

# (optional but recommended) ensure no stale CSV remains
csv_legacy <- file.path(demo_dir, "GOseq_enrichment_demo.tsv")
if (file.exists(csv_legacy)) file.remove(csv_legacy)


# -------------------- DE (filtered to mapped transcripts) -------------------
if (!is.na(de_in)) {

  de_full <- read_table_auto(de_in) %>% ensure_gene_id()

  # optional sanity print (remove later)
  message("[EchoGO] DE columns: ", paste(head(names(de_full), 12), collapse = ", "))

  de_demo <- de_full %>%
    dplyr::filter(.data$gene_id %in% mapped_ids) %>%
    dplyr::slice_head(n = N_DE)

  readr::write_tsv(de_demo, file.path(demo_dir, "DE_results_demo.tsv"))

} else {
  warning("DE input not found under src_input; skipping DE_results_demo.tsv")
}


# -------------------- COUNTS (filtered to mapped transcripts) ----------------
if (!is.na(cts_in)) {

  cts_full <- read_table_auto(cts_in) %>% ensure_gene_id()

  message("[EchoGO] Counts columns: ", paste(head(names(cts_full), 12), collapse = ", "))

  cts_demo <- cts_full %>%
    dplyr::filter(.data$gene_id %in% mapped_ids) %>%
    dplyr::slice_head(n = N_COUNTS)

  cts_demo <- repair_count_columns(cts_demo)

  sample_cols <- setdiff(names(cts_demo), "gene_id")
  keep_samps  <- head(sample_cols, min(N_SAMPLES, length(sample_cols)))
  cts_demo    <- cts_demo[, c("gene_id", keep_samps), drop = FALSE]

  readr::write_tsv(cts_demo, file.path(demo_dir, "counts_demo.tsv"))

} else {
  warning("Counts input not found under src_input; skipping counts_demo.tsv")
}

# -------------------- README -------------------------------------------------
readme <- c(
  "EchoGO demo dataset (distilled from HOL03_test full run)",
  "--------------------------------------------------------",
  paste0("Source results: ", normalizePath(src_results, winslash = "/")),
  paste0("Source input:   ", normalizePath(src_input, winslash = "/")),
  "",
  "Files:",
  " - GOseq_enrichment_demo.tsv : GOseq subset with gene_ids RESTRICTED to mapped transcripts",
  " - Trinotate_demo.tsv        : TRINITY -> SYMBOL mapping (synthetic, Metazoa-tagged)",
  " - DE_results_demo.tsv       : DE subset for mapped transcripts (if DE source found)",
  " - counts_demo.tsv           : Counts subset for mapped transcripts (if counts source found)",
  "",
  "Run the demo:",
  " demo_dir <- system.file('extdata','echogo_demo', package='EchoGO')",
  " echogo_quickstart(run_demo = TRUE)"
)
writeLines(readme, file.path(demo_dir, "README_demo.txt"))

message("[EchoGO] Demo distilled to: ", normalizePath(demo_dir, winslash = "/"))
message("[EchoGO] Inputs detected:")
message("  counts: ", cts_in)
message("  DE:     ", de_in)
message("  GOseq:  ", goseq_in)
message("[EchoGO] Outputs used:")
message("  consensus: ", cons_xlsx)
message("  goseq:     ", goseq_csv)

chk <- readr::read_tsv(file.path(demo_dir, "GOseq_enrichment_demo.tsv"), show_col_types = FALSE)
stopifnot(!any(grepl("^V\\d+$", names(chk))))
stopifnot(all(c("category","term","ontology","numDEInCat","numInCat","over_represented_FDR","gene_ids") %in% names(chk)))

#' @name load_and_annotate_goseq
#' @title Load and annotate GOseq enrichment results
#' @description
#' Loads GOseq results, maps gene names using Trinotate annotations (Metazoa-biased heuristic),
#' computes fold enrichment and GO term depth, and exports cleaned results.
#'
#' @param goseq_file Path to the GOseq enrichment file (TSV/CSV with columns including
#'   \code{category}, \code{term}, \code{ontology}, \code{numDEInCat}, \code{numInCat},
#'   \code{over_represented_FDR}, and \code{gene_ids}).
#' @param trinotate_file Path to the Trinotate report (TSV/XLS; minimally needs
#'   \code{transcript_id}, and ideally \code{sprot_Top_BLASTX_hit}, \code{EggNM.Preferred_name},
#'   \code{EggNM.max_annot_lvl}).
#' @param de_file Path to the DE results file (filtered DE list used by GOseq).
#' @param count_matrix_file Path to the full count matrix (rows = background size).
#' @param output_dir Output folder to save results (default: \code{"goseq"}).
#'   If \code{options(EchoGO.legacy_aliases)=TRUE}, a legacy mirror is also written to
#'   \code{"orthology_based_enrichment_support"}.
#' @return A data.frame of enriched GO terms with annotations and depth.
#' @export
#' @importFrom dplyr %>% select filter mutate arrange rename left_join distinct transmute
#' @importFrom stringr str_detect
#' @importFrom readr write_csv
#' @importFrom openxlsx write.xlsx
#' @importFrom utils read.delim
load_and_annotate_goseq <- function(
    goseq_file,
    trinotate_file,
    de_file,
    count_matrix_file,
    output_dir = "goseq"
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  .mk <- function(...) { p <- file.path(...); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
  .mirror_legacy <- function(src_root, legacy_root) {
    if (!nzchar(src_root) || !nzchar(legacy_root) || !dir.exists(src_root)) return(invisible(NULL))
    dir.create(legacy_root, recursive = TRUE, showWarnings = FALSE)
    src_files <- list.files(src_root, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
    for (f in src_files) {
      if (dir.exists(f)) next
      rel <- sub(paste0("^", gsub("\\\\","\\\\\\\\", normalizePath(src_root, winslash="/", mustWork=FALSE))), "", normalizePath(f, winslash="/", mustWork=FALSE))
      rel <- sub("^[/\\\\]", "", rel)
      dest <- file.path(legacy_root, rel)
      dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
      file.copy(f, dest, overwrite = TRUE)
    }
    cat("This legacy folder mirrors: ", basename(src_root), "/\n", file = file.path(legacy_root, "__moved_to.txt"))
    invisible(NULL)
  }

  # Respect global toggle for legacy mirrors
  legacy_on <- isTRUE(getOption("EchoGO.legacy_aliases", FALSE))

  # Anchor relative output_dir to active results dir if provided
  if (!grepl("^([A-Za-z]:)?[\\/]", output_dir)) {
    hinted <- getOption("EchoGO.active_results_dir", NULL)
    if (!is.null(hinted) && nzchar(hinted)) {
      if (identical(output_dir, "orthology_based_enrichment_support")) {
        output_dir <- file.path(hinted, "goseq")
      } else {
        output_dir <- file.path(hinted, output_dir)
      }
    }
  }

  # ---- I/O setup (canonical) ----
  output_dir <- normalizePath(.mk(output_dir), winslash = "/", mustWork = FALSE)

  # ---- Load GOseq enrichment ----
  df <- utils::read.delim(
    goseq_file,
    sep = if (grepl("\\.csv$", goseq_file, ignore.case = TRUE)) "," else "\t",
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!"category" %in% names(df)) stop("GOseq file missing required column: 'category'")
  if (!"term" %in% names(df))     stop("GOseq file missing required column: 'term'")
  if (!"ontology" %in% names(df)) stop("GOseq file missing required column: 'ontology'")
  if (!"numDEInCat" %in% names(df) || !"numInCat" %in% names(df))
    stop("GOseq file must include 'numDEInCat' and 'numInCat'")
  if (!"over_represented_FDR" %in% names(df))
    stop("GOseq file must include 'over_represented_FDR'")
  if (!"gene_ids" %in% names(df))
    warning("GOseq file has no 'gene_ids' column; gene name mapping will be empty.")

  df$clean_go_term <- trimws(df$category)

  # ---- Load Trinotate and derive primary names (safe handling of fields) ----
  tri <- utils::read.delim(
    trinotate_file,
    sep = if (grepl("\\.csv$", trinotate_file, ignore.case = TRUE)) "," else "\t",
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!nrow(tri)) {
    tri <- tri[0, , drop = FALSE]
  }
  tri[tri == "."] <- NA

  # blast taxonomy (if present)
  if ("sprot_Top_BLASTX_hit" %in% names(tri)) {
    tri$blast_taxonomy <- vapply(
      strsplit(as.character(tri$sprot_Top_BLASTX_hit), "\\^", fixed = FALSE),
      function(x) if (length(x)) tail(x, 1) else NA_character_, character(1)
    )
  } else {
    tri$blast_taxonomy <- rep(NA_character_, nrow(tri))
  }

  # Metazoan heuristic (EggNOG rank OR BLASTX taxonomy mentions Metazoa)
  animal_taxid_cutoff <- 33208
  egg_col  <- "EggNM.max_annot_lvl"
  pref_col <- "EggNM.Preferred_name"
  is_animal_eggnog <- egg_col %in% names(tri) &
    !is.na(suppressWarnings(as.numeric(tri[[egg_col]]))) &
    suppressWarnings(as.numeric(tri[[egg_col]])) >= animal_taxid_cutoff
  is_animal_blastx <- stringr::str_detect(tri$blast_taxonomy %||% "", "Metazoa")
  tri_animal <- tri[is_animal_eggnog | is_animal_blastx, , drop = FALSE]

  # Primary display name preference (length-safe assigns)
  if (!"transcript_id" %in% names(tri_animal)) {
    tri_animal$transcript_id <- rep(NA_character_, nrow(tri_animal))
  }
  tri_animal$primary_name <- rep(NA_character_, nrow(tri_animal))
  if ("sprot_Top_BLASTX_hit" %in% names(tri_animal)) {
    base <- sub("\\^.*", "", tri_animal$sprot_Top_BLASTX_hit)
    tri_animal$primary_name <- sub("_.*", "", base)
  }
  if (pref_col %in% names(tri_animal)) {
    use_pref <- is.na(tri_animal$primary_name) & !is.na(tri_animal[[pref_col]])
    tri_animal$primary_name[use_pref] <- tri_animal[[pref_col]][use_pref]
  }
  tri_animal$primary_name[is.na(tri_animal$primary_name)] <- tri_animal$transcript_id[is.na(tri_animal$primary_name)]

  transcript_map <- dplyr::select(tri_animal, transcript_id, primary_name) %>% dplyr::distinct()

  # ---- Map gene IDs to names (robust against missing/NA 'gene_ids') ----
  if ("gene_ids" %in% names(df)) {
    df$gene_names <- vapply(df$gene_ids, function(glist) {
      if (is.na(glist) || !nzchar(glist)) return("")
      ids <- unique(strsplit(as.character(glist), ",\\s*")[[1]])
      if (!length(ids)) return("")
      mapped <- transcript_map[transcript_map$transcript_id %in% ids, , drop = FALSE]
      paste(unique(stats::na.omit(mapped$primary_name)), collapse = ", ")
    }, character(1))
  } else {
    df$gene_names <- ""
  }

  # ---- Load DE + background and compute fold enrichment ----
  de_results   <- utils::read.delim(de_file, stringsAsFactors = FALSE)
  count_matrix <- utils::read.delim(count_matrix_file, stringsAsFactors = FALSE, check.names = FALSE)
  totalDE <- nrow(de_results)
  totalBG <- nrow(count_matrix)
  if (totalDE <= 0 || totalBG <= 0) warning("DE or background size is zero; foldEnrichment will be NA/Inf.")
  df$foldEnrichment <- with(df, (numDEInCat / max(totalDE, 1)) / (numInCat / max(totalBG, 1)))

  # ---- Remove empty-name rows (unchanged behavior) ----
  df <- df %>% dplyr::filter(!is.na(.data$gene_names), trimws(.data$gene_names) != "")

  # ---- Add GO term depth (prefer GO.db, fallback to OBO) ----
  normalize_goid <- function(x) {
    x <- as.character(x)
    m <- regmatches(x, regexpr("GO:\\d{7}", x))
    ifelse(nzchar(m), m, NA_character_)
  }

  goids <- normalize_goid(df$clean_go_term)
  depth_vec <- rep(NA_integer_, length(goids))
  used <- "none"

  if (requireNamespace("GO.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
    # FIX: avoid passing NA keys to mget()
    valid <- !is.na(goids)
    keys  <- unique(goids[valid])
    safe_mget <- function(keys, map) {
      if (!length(keys)) return(setNames(vector("list", 0), character(0)))
      AnnotationDbi::mget(keys, map, ifnotfound = NA)
    }
    anc_bp_lu <- safe_mget(keys, GO.db::GOBPANCESTOR)
    anc_mf_lu <- safe_mget(keys, GO.db::GOMFANCESTOR)
    anc_cc_lu <- safe_mget(keys, GO.db::GOCCANCESTOR)

    depth_vec <- vapply(seq_along(goids), function(i) {
      gi <- goids[i]
      if (is.na(gi)) return(NA_integer_)
      ai_raw <- c(
        unlist(anc_bp_lu[[gi]], use.names = FALSE),
        unlist(anc_mf_lu[[gi]], use.names = FALSE),
        unlist(anc_cc_lu[[gi]], use.names = FALSE)
      )
      ai <- ai_raw[!is.na(ai_raw) & ai_raw != "all"]
      if (!length(ai)) 0L else length(unique(ai))
    }, integer(1))

    used <- "GO.db"
  } else {
    # Fallback: ontologyIndex with bundled/cached OBO
    go_obo <- system.file("extdata","go-basic.obo", package = "EchoGO")
    if (go_obo == "") {
      go_obo <- file.path(tempdir(), "go-basic.obo")
      if (!file.exists(go_obo)) {
        message("Caching GO OBO locally for depth: ", go_obo)
        utils::download.file("http://purl.obolibrary.org/obo/go.obo",
                             destfile = go_obo, mode = "wb", quiet = TRUE)
      }
    }
    if (requireNamespace("ontologyIndex", quietly = TRUE) && file.exists(go_obo)) {
      go_ont <- ontologyIndex::get_ontology(go_obo, extract_tags = "minimal")
      depth_vec <- vapply(goids, function(term_id) {
        if (!is.na(term_id) && term_id %in% go_ont$id) {
          length(ontologyIndex::get_ancestors(go_ont, term_id))
        } else NA_integer_
      }, integer(1))
      used <- "ontologyIndex"
    }
  }

  df$depth <- depth_vec
  message(switch(used,
                 "GO.db"         = "Using GO.db for GO depth.",
                 "ontologyIndex" = "Using bundled/cached go-basic.obo for GO depth (ontologyIndex).",
                 "none"          = "No GO depth backend available; depth set to NA."
  ))

  # ---- Ensure 'all_genes' for downstream modules ----
  df$all_genes <- df$gene_names

  # ---- Export (canonical) ----
  readr::write_csv(df, file.path(output_dir, "GOseq_enrichment_full_annotated.csv"))
  openxlsx::write.xlsx(df, file.path(output_dir, "GOseq_enrichment_full_annotated.xlsx"), overwrite = TRUE)

  supp_table <- df %>%
    dplyr::select(category, term, ontology, numDEInCat, numInCat,
                  foldEnrichment, over_represented_FDR, depth, gene_names)
  readr::write_csv(supp_table, file.path(output_dir, "GO_enrichment_supplementary_clean.csv"))
  openxlsx::write.xlsx(supp_table, file.path(output_dir, "GO_enrichment_supplementary_clean.xlsx"), overwrite = TRUE)

  # ---- Optional legacy mirror (only if enabled) ----
  if (legacy_on) {
    legacy_dir <- file.path(dirname(output_dir), "orthology_based_enrichment_support")
    .mirror_legacy(src_root = output_dir, legacy_root = legacy_dir)
  }

  return(df)
}

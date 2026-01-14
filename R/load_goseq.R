#' @name load_and_annotate_goseq
#' @title Load and annotate GOseq enrichment results
#' @description
#' Loads GOseq results, maps gene names using Trinotate annotations (Metazoa-biased heuristic),
#' computes fold enrichment and GO term depth, and exports cleaned results.
#' #' In reference-based RNA-seq workflows, GOseq outputs may already contain gene identifiers
#' that are gene symbols/names (not transcript IDs). If the identifiers in \code{gene_ids}
#' do not match the annotation \code{transcript_id} column, EchoGO will treat \code{gene_ids}
#' as already-usable names and use them directly.
#'
#' @param goseq_file Path to the GOseq enrichment file (TSV/CSV with columns including
#'   \code{category}, \code{term}, \code{ontology}, \code{numDEInCat}, \code{numInCat},
#'   \code{over_represented_FDR}, and \code{gene_ids}).
#' @param trinotate_file Path to the Trinotate report (TSV/XLS; minimally needs
#'   \code{transcript_id}, and ideally \code{sprot_Top_BLASTX_hit}, \code{EggNM.Preferred_name},
#'   \code{EggNM.max_annot_lvl}).
#' @param de_file Path to the DE results file (filtered DE list used by GOseq).
#'   In reference-based mode this is optional (used mainly for totals / logging);
#'   if provided, it should correspond to the same contrast as \code{goseq_file}.
#' @param count_matrix_file Path to the full count matrix (rows = background size).
#'   Optional; used to estimate background size for fold enrichment. If missing, EchoGO
#'   may approximate fold enrichment from GOseq columns or leave it as-is depending on inputs.
#' @param output_dir Output folder to save results (default: \code{"goseq"}).
#'   If \code{options(EchoGO.legacy_aliases)=TRUE}, a legacy mirror is also written to
#'   \code{"orthology_based_enrichment_support"}.
#' @return A data.frame of enriched GO terms with annotations and depth.
#' #' The returned data.frame also carries an attribute \code{echogo_bg_universe} when it can be
#' derived from Trinotate/EggNOG primary names (used downstream as a background universe).
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

  # ---- Robust readers -------------------------------------------------

  .detect_delim <- function(file, candidates = c("\t",";",","), n = 50L, expected = NULL) {

    # 1) Header-based scoring (most reliable for GOseq)
    hdr <- tryCatch(readLines(file, n = 1L, warn = FALSE), error = function(e) "")
    hdr <- enc2utf8(hdr)

    score_header <- function(sep) {
      parts <- strsplit(hdr, sep, fixed = TRUE)[[1]]
      parts <- gsub("^\ufeff", "", parts)          # BOM
      parts <- gsub("\u00A0", " ", parts, fixed=TRUE)
      parts <- trimws(parts)
      parts_lc <- tolower(parts)

      # if separator doesn't split header, this sep is wrong
      if (length(parts_lc) < 2) return(-Inf)

      if (!is.null(expected)) {
        sum(parts_lc %in% tolower(expected))
      } else {
        0
      }
    }

    if (!is.null(expected)) {
      hs <- vapply(candidates, score_header, numeric(1))
      best_h <- max(hs, na.rm = TRUE)
      if (is.finite(best_h) && best_h > 0) {
        # pick the candidate with best header match (tie -> first)
        return(candidates[which.max(hs)])
      }
    }

    # 2) Fallback: count.fields heuristic (for generic files)
    best <- candidates[1]
    best_score <- -Inf

    for (sep in candidates) {
      cf <- tryCatch(
        utils::count.fields(file, sep = sep, quote = "", comment.char = "", skipNul = TRUE),
        error = function(e) integer(0)
      )
      if (!length(cf)) next
      cf <- cf[seq_len(min(length(cf), n))]
      med <- stats::median(cf)
      if (is.na(med) || med < 2) next

      # penalize crazy-wide parses (e.g., splitting gene lists on comma)
      penalty <- if (med > 50) 50 else 0
      score <- (med - stats::sd(cf)) - penalty

      if (is.finite(score) && score > best_score) {
        best_score <- score
        best <- sep
      }
    }

    best
  }


  .read_delim_robust <- function(file, candidates = c("\t",";",","), dec=".", expected = NULL, ...) {
    sep <- .detect_delim(file, candidates = candidates, expected = expected)

    # Special case: csv2-style numbers (sep=';' and decimals use ',')
    if (identical(sep, ";")) {
      first <- readLines(file, n = 2L, warn = FALSE)
      # if we see patterns like "0,123" and NOT "0.123", assume dec=","
      if (any(grepl("\\d,\\d", first)) && !any(grepl("\\d\\.\\d", first))) dec <- ","
    }

    # First try: base read.delim with permissive settings
    out <- tryCatch(
      utils::read.delim(
        file,
        sep = sep,
        dec = dec,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        quote = "",
        fill = TRUE,
        comment.char = "",
        ...
      ),
      error = function(e) NULL
    )
    if (!is.null(out)) return(out)

    # Second try: data.table::fread (more forgiving) if available
    if (requireNamespace("data.table", quietly = TRUE)) {
      out2 <- tryCatch(
        data.table::fread(
          file,
          sep = sep,
          dec = dec,
          data.table = FALSE,
          fill = TRUE,
          quote = "",
          showProgress = FALSE
        ),
        error = function(e) NULL
      )
      if (!is.null(out2)) return(out2)
    }

    stop("Failed to read file robustly: ", file)
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

  # ---- Load GOseq enrichment (ROBUST) ----
  df <- .read_delim_robust(
    goseq_file,
    candidates = c("\t",";",","),
    expected = c(
      "category","term","ontology","numdeincat","numincat",
      "over_represented_fdr","gene_ids",
      "over_represented_pvalue","under_represented_pvalue",
      "under_represented_fdr","go_term"
    )
  )
  message("GOseq columns: ", paste(names(df), collapse = " | "))
  message("Classes: ", paste(sapply(df[, c("numDEInCat","numInCat","over_represented_FDR")], class), collapse = ", "))

  # ---- Sanitize GOseq column names (BOM/whitespace/duplicates) ----
  names(df) <- gsub("^\ufeff", "", names(df))  # remove UTF-8 BOM if present
  names(df) <- trimws(names(df))              # remove leading/trailing spaces
  names(df) <- make.unique(names(df))         # avoid duplicate names (e.g. 'term' + ' term')

  # ---- Normalize GOseq column names (case-insensitive + synonyms) ----
  names_lc <- tolower(names(df))

  .pick_col <- function(cands) {
    w <- match(tolower(cands), names_lc)
    w <- w[!is.na(w)][1]
    if (is.na(w)) return(NULL)
    names(df)[w]
  }

  # If 'category' is missing, try re-reading GOseq with each candidate sep
  cat_col <- .pick_col(c("category","clean_go_term","go_term","go","goid","term_id"))
  if (is.null(cat_col)) {
    for (sep_try in c("\t",";",",")) {
      tmp <- tryCatch(
        utils::read.delim(
          goseq_file,
          sep = sep_try,
          stringsAsFactors = FALSE,
          check.names = FALSE,
          quote = "",
          fill = TRUE,
          comment.char = ""
        ),
        error = function(e) NULL
      )
      if (is.null(tmp) || !nrow(tmp)) next
      tmp_names_lc <- tolower(names(tmp))
      if (any(tmp_names_lc %in% c("category","clean_go_term","go_term","go","goid","term_id"))) {
        df <- tmp
        names_lc <- tolower(names(df))
        break
      }
    }
  }

  # Re-pick after possible re-read
  cat_col <- .pick_col(c("category","clean_go_term","go_term","go","goid","term_id"))
  term_col <- .pick_col(c("term","name","description","term_name"))
  ont_col  <- .pick_col(c("ontology","ont"))
  nde_col  <- .pick_col(c("numdeincat","num_de_in_cat","numdeincategory","num_sig_in_cat","n_de"))
  nin_col  <- .pick_col(c("numincat","num_in_cat","numincategory","n_in_cat"))
  fdr_col  <- .pick_col(c("over_represented_fdr","over_represented_fdr_adj","fdr","padj","qvalue"))
  gid_col  <- .pick_col(c("gene_ids","gene_id","genes","geneids","genes_in_cat"))

  # Create canonical columns expected downstream
  if (!is.null(cat_col) && cat_col != "category") df$category <- df[[cat_col]]
  if (is.null(term_col)) df$term <- as.character(df$category) else if (term_col != "term") df$term <- df[[term_col]]
  if (!is.null(ont_col) && ont_col != "ontology") df$ontology <- df[[ont_col]]
  if (!is.null(nde_col) && nde_col != "numDEInCat") df$numDEInCat <- df[[nde_col]]
  if (!is.null(nin_col) && nin_col != "numInCat")  df$numInCat  <- df[[nin_col]]
  if (!is.null(fdr_col) && fdr_col != "over_represented_FDR") df$over_represented_FDR <- df[[fdr_col]]
  if (!is.null(gid_col) && gid_col != "gene_ids") df$gene_ids <- df[[gid_col]]

  names_lc <- tolower(names(df))

  # Now validate using canonical names
  if (!"category" %in% names(df)) {
    stop("GOseq file missing required column: 'category'. Columns found: ",
         paste(names(df), collapse = ", "))
  }
  if (!"ontology" %in% names(df)) {
    stop("GOseq file missing required column: 'ontology'. Columns found: ",
         paste(names(df), collapse = ", "))
  }
  if (!"numDEInCat" %in% names(df) || !"numInCat" %in% names(df)) {
    stop("GOseq file must include 'numDEInCat' and 'numInCat' (or synonyms). Columns found: ",
         paste(names(df), collapse = ", "))
  }
  if (!"over_represented_FDR" %in% names(df)) {
    stop("GOseq file must include 'over_represented_FDR' (or synonyms like FDR/padj). Columns found: ",
         paste(names(df), collapse = ", "))
  }
  if (!"gene_ids" %in% names(df)) {
    warning("GOseq file has no 'gene_ids' column; gene name mapping will be empty.")
  }

  df$clean_go_term <- trimws(df$category)

  # ---- Load Trinotate and derive primary names (ROBUST) ----
  tri <- tryCatch(
    {
      .read_delim_robust(
        trinotate_file,
        candidates = c("\t", ";", ",")
      )
    },
    error = function(e) {
      message("⚠️  Trinotate/eggNOG read failed even in robust mode: ", conditionMessage(e))
      data.frame()
    }
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

   # ---- Define EggNOG / Trinotate-based background universe (Mode A) ----
  bg_universe <- unique(stats::na.omit(transcript_map$primary_name))
  if (!length(bg_universe)) {
    message("⚠️ No primary_name universe derived from Trinotate; ",
            "g:Profiler will fall back to GOseq gene_names for background.")
  }
  # Attach as attribute so the pipeline can recover it later
  attr(df, "echogo_bg_universe") <- bg_universe

  # ---- Map gene IDs to names (robust against missing/NA 'gene_ids') ----
  if ("gene_ids" %in% names(df)) {
    df$gene_names <- vapply(df$gene_ids, function(glist) {
      if (is.na(glist) || !nzchar(glist)) return("")

      ids <- trimws(unlist(strsplit(as.character(glist), ",", fixed = TRUE)))
            ids <- ids[nzchar(ids)]
      if (!length(ids)) return("")

      # If transcript_map is empty OR none of the ids match transcript_id,
      # treat gene_ids as already-usable names (reference-based mode).
      if (!nrow(transcript_map) || !any(ids %in% transcript_map$transcript_id)) {
        return(paste(unique(ids), collapse = ", "))
      }

      mapped <- transcript_map[transcript_map$transcript_id %in% ids, , drop = FALSE]
      out <- unique(stats::na.omit(mapped$primary_name))

      # fallback if mapping produces nothing
      if (!length(out)) paste(unique(ids), collapse = ", ") else paste(out, collapse = ", ")
    }, character(1))

  } else {
    df$gene_names <- ""
  }

  # ---- Load DE results (ROBUST) ----
  de_results <- .read_delim_robust(
    de_file,
    candidates = c("\t",";",",")
  )

  # Count matrix (ROBUST): we only need background size (n rows)
  count_matrix <- tryCatch(
    {
      # detect delimiter and read just the first column
      sep_cm <- .detect_delim(count_matrix_file, candidates = c("\t",";",","))
      if (requireNamespace("data.table", quietly = TRUE)) {
        data.table::fread(
          count_matrix_file,
          sep = sep_cm,
          data.table = FALSE,
          fill = TRUE,
          quote = "",
          showProgress = FALSE,
          select = 1
        )
      } else {
        utils::read.delim(
          count_matrix_file,
          sep = sep_cm,
          stringsAsFactors = FALSE,
          check.names = FALSE,
          quote = "",
          fill = TRUE,
          comment.char = ""
        )[, 1, drop = FALSE]
      }
    },
    error = function(e) {
      message(
        "⚠️  Primary count matrix read failed (", conditionMessage(e),
        "). Falling back to line counting..."
      )
      ln <- readLines(count_matrix_file, warn = FALSE)
      ln <- ln[nzchar(trimws(ln))]
      if (!length(ln)) return(data.frame(dummy = integer(0)))
      header_guess <- grepl("[A-Za-z]", ln[1])
      n <- length(ln) - if (header_guess) 1L else 0L
      if (n < 0L) n <- 0L
      data.frame(dummy = seq_len(n))
    }
  )

  totalDE <- nrow(de_results)
  totalBG <- nrow(count_matrix)

  # ---- Coerce GOseq counts + FDR to numeric (MUST happen before foldEnrichment) ----
  .to_num <- function(x) {
    if (is.numeric(x)) return(x)
    x <- as.character(x)
    x <- trimws(x)
    x <- gsub(",", ".", x, fixed = TRUE)          # tolerate comma decimals
    x <- gsub("[^0-9eE+\\-\\.]", "", x)           # strip stray text
    suppressWarnings(as.numeric(x))
  }

  df$numDEInCat <- .to_num(df$numDEInCat)
  df$numInCat   <- .to_num(df$numInCat)

  # FDR can be NA for some rows; still coerce safely
  df$over_represented_FDR <- .to_num(df$over_represented_FDR)

  if (all(is.na(df$numDEInCat)) || all(is.na(df$numInCat))) {
    stop("numDEInCat/numInCat could not be parsed as numeric. Columns present: ",
         paste(names(df), collapse = ", "))
  }


  if (totalDE <= 0 || totalBG <= 0) {
    warning("DE or background size is zero; foldEnrichment will be NA/Inf.")
  }

  df$foldEnrichment <- with(
    df,
    (numDEInCat / max(totalDE, 1)) / (numInCat / max(totalBG, 1))
  )


  # ---- Remove empty-name rows (unchanged behavior) ----
  df$gene_names <- ifelse(is.na(df$gene_names), "", df$gene_names)
  df$gene_names <- trimws(df$gene_names)

  # only drop if still empty after all fallbacks
  df <- dplyr::filter(df, nzchar(.data$gene_names))


  # ---- Add GO term depth (length-safe) ----
  normalize_goid <- function(x) {
    x <- as.character(x)
    m <- regmatches(x, regexpr("GO:\\d{7}", x))
    ifelse(nzchar(m), m, NA_character_)
  }

  goids <- normalize_goid(df$clean_go_term)

  # CRITICAL: depth_vec must match nrow(df), not length(goids) (even if something goes odd)
  depth_vec <- rep(NA_integer_, nrow(df))
  used <- "none"

  if (requireNamespace("GO.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)) {

    valid <- !is.na(goids)

    if (any(valid)) {
      keys <- unique(goids[valid])

      safe_mget <- function(keys, map) {
        if (!length(keys)) return(setNames(vector("list", 0), character(0)))
        AnnotationDbi::mget(keys, map, ifnotfound = NA)
      }

      anc_bp_lu <- safe_mget(keys, GO.db::GOBPANCESTOR)
      anc_mf_lu <- safe_mget(keys, GO.db::GOMFANCESTOR)
      anc_cc_lu <- safe_mget(keys, GO.db::GOCCANCESTOR)

      depth_vec[valid] <- vapply(goids[valid], function(gi) {
        ai_raw <- c(
          unlist(anc_bp_lu[[gi]], use.names = FALSE),
          unlist(anc_mf_lu[[gi]], use.names = FALSE),
          unlist(anc_cc_lu[[gi]], use.names = FALSE)
        )
        ai <- ai_raw[!is.na(ai_raw) & ai_raw != "all"]
        if (!length(ai)) 0L else length(unique(ai))
      }, integer(1))
    }

    used <- "GO.db"

  } else {
    # Fallback: ontologyIndex with bundled/cached OBO
    go_obo <- system.file("extdata", "go-basic.obo", package = "EchoGO")
    if (go_obo == "") {
      go_obo <- file.path(tempdir(), "go-basic.obo")
      if (!file.exists(go_obo)) {
        utils::download.file("http://purl.obolibrary.org/obo/go.obo",
                             destfile = go_obo, mode = "wb", quiet = TRUE)
      }
    }

    if (requireNamespace("ontologyIndex", quietly = TRUE) && file.exists(go_obo)) {
      go_ont <- ontologyIndex::get_ontology(go_obo, extract_tags = "minimal")
      valid <- !is.na(goids)

      if (any(valid)) {
        depth_vec[valid] <- vapply(goids[valid], function(term_id) {
          if (!is.na(term_id) && term_id %in% go_ont$id) {
            length(ontologyIndex::get_ancestors(go_ont, term_id))
          } else NA_integer_
        }, integer(1))
      }

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

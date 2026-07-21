.echogo_compute_goseq_fold <- function(
    df,
    total_significant_genes,
    total_tested_genes
) {
  scalar_positive_integer <- function(x, label) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
        x <= 0 || x != as.integer(x)) {
      stop(label, " must be one positive integer.", call. = FALSE)
    }
    as.integer(x)
  }
  total_significant_genes <- scalar_positive_integer(
    total_significant_genes, "total_significant_genes"
  )
  total_tested_genes <- scalar_positive_integer(total_tested_genes, "total_tested_genes")
  if (total_significant_genes > total_tested_genes) {
    stop("total_significant_genes cannot exceed total_tested_genes.", call. = FALSE)
  }
  if (!all(c("numDEInCat", "numInCat") %in% names(df))) {
    stop("GOseq data require numDEInCat and numInCat columns.", call. = FALSE)
  }
  num_de <- suppressWarnings(as.numeric(df$numDEInCat))
  num_in <- suppressWarnings(as.numeric(df$numInCat))
  if (any(is.na(num_de)) || any(is.na(num_in)) || any(num_de < 0) || any(num_in < 0)) {
    stop("GOseq category counts must be finite non-negative numbers.", call. = FALSE)
  }
  if (any(num_de > total_significant_genes)) {
    stop("GOseq numDEInCat exceeds total_significant_genes.", call. = FALSE)
  }
  if (any(num_in > total_tested_genes)) {
    stop("GOseq numInCat exceeds total_tested_genes.", call. = FALSE)
  }
  if (any(num_de > num_in)) {
    stop("GOseq numDEInCat exceeds numInCat; the universes are incompatible.", call. = FALSE)
  }
  df$numDEInCat <- num_de
  df$numInCat <- num_in
  df$total_significant_genes <- rep(total_significant_genes, nrow(df))
  df$total_tested_genes <- rep(total_tested_genes, nrow(df))
  df$foldEnrichment <- (num_de / total_significant_genes) /
    (num_in / total_tested_genes)
  df$foldEnrichment[num_in == 0 & num_de == 0] <- NA_real_
  df
}

#' @name load_and_annotate_goseq
#' @title Load and annotate GOseq enrichment results
#' @description
#' Loads GOseq results, adds display names and annotation provenance, computes
#' fold enrichment and GO term depth, and exports cleaned results. In reference-based
#' workflows, existing GOseq gene symbols are preserved for display. g:Profiler
#' submission remains governed by the portable canonical-name resolver.
#'
#' @param goseq_file Path to the GOseq enrichment file (TSV/CSV with columns including
#'   \code{category}, \code{term}, \code{ontology}, \code{numDEInCat}, \code{numInCat},
#'   \code{over_represented_FDR}, and \code{gene_ids}).
#' @param trinotate_file Path to the Trinotate report (TSV/XLS; minimally needs
#'   \code{transcript_id}, and ideally \code{sprot_Top_BLASTX_hit}, \code{EggNM.Preferred_name},
#'   \code{EggNM.max_annot_lvl}).
#' @param de_file Path to the contrast-specific DE results table.
#' @param count_matrix_file Path to the matching tested-gene count matrix.
#' @param output_dir Output folder to save results (default: \code{"goseq"}). If
#'   \code{options(EchoGO.legacy_aliases)=TRUE}, a compatibility mirror is written.
#' @param gene_sets Optional object returned by `prepare_gprofiler_gene_sets()`.
#' @param total_significant_genes,total_tested_genes Optional explicit GOseq denominators.
#' @param de_id_column,de_significant_column,de_padj_column,de_lfc_column Optional DE columns.
#' @param count_id_column,annotation_id_column Optional count/annotation ID columns.
#' @param padj_threshold,log2fc_threshold Default DE significance thresholds.
#' @param de_table_significant_only Declare that the DE table contains significant rows only.
#' @param use_trinotate_universe Deprecated compatibility flag. It does not
#'   intersect the tested universe or permit raw-ID fallback.
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
    output_dir = "goseq",
    gene_sets = NULL,
    total_significant_genes = NULL,
    total_tested_genes = NULL,
    de_id_column = NULL,
    de_significant_column = NULL,
    de_padj_column = NULL,
    de_lfc_column = NULL,
    count_id_column = NULL,
    annotation_id_column = NULL,
    padj_threshold = 0.05,
    log2fc_threshold = 1,
    de_table_significant_only = FALSE,
    use_trinotate_universe = FALSE
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

  # Package-level readers keep parsing behavior consistent and testable.
  .detect_delim <- .echogo_detect_delim
  .read_delim_robust <- .echogo_read_delim_robust

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

  # ---- Sanitize GOseq column names (BOM/whitespace/duplicates) ----
  names(df) <- gsub("^\ufeff", "", names(df))  # remove UTF-8 BOM if present
  names(df) <- trimws(names(df))              # remove leading/trailing spaces
  names(df) <- make.unique(names(df))         # avoid duplicate names (e.g. 'term' + ' term')
  .echogo_check_goseq_parse(df)

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
          quote = "\"",
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
  numeric_columns <- intersect(
    c("numDEInCat", "numInCat", "over_represented_FDR"),
    names(df)
  )
  if (length(numeric_columns)) {
    message(
      "Classes: ",
      paste(vapply(df[numeric_columns], function(x) class(x)[1], character(1)), collapse = ", ")
    )
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

  # ---- Load Trinotate and derive display names ----
  tri <- tryCatch(
    {
      .read_delim_robust(
        trinotate_file,
        candidates = c("\t", ";", ",")
      )
    },
    error = function(e) {
      message("Warning: Trinotate/eggNOG input could not be read: ", conditionMessage(e))
      data.frame()
    }
  )

  if (!nrow(tri)) {
    tri <- tri[0, , drop = FALSE]
  }
  tri[tri == "."] <- NA

  if (is.null(gene_sets)) {
    gene_sets <- prepare_gprofiler_gene_sets(
      de_results = de_file,
      annotation = tri,
      count_matrix = count_matrix_file,
      de_id_column = de_id_column,
      significant_column = de_significant_column,
      padj_column = de_padj_column,
      lfc_column = de_lfc_column,
      count_id_column = count_id_column,
      annotation_id_column = annotation_id_column,
      padj_threshold = padj_threshold,
      log2fc_threshold = log2fc_threshold,
      de_table_significant_only = de_table_significant_only,
      use_trinotate_universe = use_trinotate_universe
    )
  }

  # blast taxonomy (if present)
  if ("sprot_Top_BLASTX_hit" %in% names(tri)) {
    tri$blast_taxonomy <- vapply(
      strsplit(as.character(tri$sprot_Top_BLASTX_hit), "\\^", fixed = FALSE),
      function(x) if (length(x)) tail(x, 1) else NA_character_, character(1)
    )
  } else {
    tri$blast_taxonomy <- rep(NA_character_, nrow(tri))
  }

  # Taxonomy is descriptive annotation provenance. Eligibility is based on
  # explicit lineage labels, never on numeric taxon ordering (taxon IDs are not
  # ordinal ranks).
  egg_col  <- "EggNM.max_annot_lvl"
  pref_col <- "EggNM.Preferred_name"
  is_animal_eggnog <- egg_col %in% names(tri) &
    stringr::str_detect(as.character(tri[[egg_col]]) %||% "", "(?i)Metazoa")
  is_animal_blastx <- stringr::str_detect(tri$blast_taxonomy %||% "", "(?i)Metazoa")
  is_animal_eggnog[is.na(is_animal_eggnog)] <- FALSE
  is_animal_blastx[is.na(is_animal_blastx)] <- FALSE
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

  contract_map <- gene_sets$mapping_table
  contract_map <- contract_map[
    !is.na(contract_map$resolved_name) & nzchar(contract_map$resolved_name),
    c("original_id", "resolved_name"),
    drop = FALSE
  ]
  if (nrow(contract_map)) {
    transcript_map <- data.frame(
      transcript_id = contract_map$original_id,
      primary_name = contract_map$resolved_name,
      stringsAsFactors = FALSE
    )
  }

  # Preserve resolved annotation names for diagnostics only. They never replace
  # the tested experiment universe used for a custom background.
  attr(df, "echogo_annotation_names") <- unique(stats::na.omit(transcript_map$primary_name))

  # ---- Map gene IDs to display names ----
  if ("gene_ids" %in% names(df)) {
    df$gene_names <- vapply(df$gene_ids, function(glist) {
      if (is.na(glist) || !nzchar(glist)) return("")

      ids <- trimws(unlist(strsplit(as.character(glist), ",", fixed = TRUE)))
            ids <- ids[nzchar(ids)]
      if (!length(ids)) return("")

      # Preserve existing GOseq display names in reference-based mode when no
      # transcript mapping applies.
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

  if ("gene_ids" %in% names(df)) {
    goseq_ids <- unique(unlist(strsplit(
      paste(stats::na.omit(as.character(df$gene_ids)), collapse = ","),
      ",",
      fixed = TRUE
    )))
    goseq_ids <- trimws(goseq_ids)
    goseq_ids <- goseq_ids[nzchar(goseq_ids)]
    incompatible_ids <- setdiff(goseq_ids, gene_sets$foreground_original)
    if (length(incompatible_ids)) {
      stop(
        "The GOseq table contains gene IDs outside the significant-gene universe: ",
        paste(utils::head(incompatible_ids, 8L), collapse = ", "),
        ". Confirm that GOseq and the DE contrast use the same identifiers and filtering universe.",
        call. = FALSE
      )
    }
  }

  # Denominators come from explicit GOseq metadata or the matched DE gene sets.
  constant_integer <- function(column, label) {
    if (!column %in% names(df)) return(NULL)
    values <- unique(stats::na.omit(suppressWarnings(as.numeric(df[[column]]))))
    if (length(values) != 1L) {
      stop("GOseq ", label, " metadata must contain one constant value.", call. = FALSE)
    }
    as.integer(values[[1]])
  }
  totalDE <- total_significant_genes %||%
    constant_integer("total_significant_genes", "significant-gene denominator") %||%
    length(gene_sets$foreground_original)
  totalBG <- total_tested_genes %||%
    constant_integer("total_tested_genes", "tested-gene denominator") %||%
    length(gene_sets$background_original)

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


  df <- .echogo_compute_goseq_fold(
    df,
    total_significant_genes = totalDE,
    total_tested_genes = totalBG
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
                  total_significant_genes, total_tested_genes,
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

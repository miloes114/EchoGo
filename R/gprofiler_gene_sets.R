# Matched portable gene-set preparation --------------------------------------

.echogo_gene_set_table <- function(x, label) {
  if (is.data.frame(x)) return(as.data.frame(x, stringsAsFactors = FALSE))
  if (is.character(x) && length(x) == 1L && file.exists(x)) {
    return(.echogo_read_delim_robust(
      x,
      candidates = c("\t", ",", ";"),
      expected = c(
        "gene_id", "id", "gene", "transcript_id", "padj", "fdr",
        "log2foldchange", "significant", "sprot_top_blastx_hit",
        "eggnm.preferred_name", "preferred_name"
      )
    ))
  }
  stop(label, " must be a data frame or an existing table path.", call. = FALSE)
}

.echogo_clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  missing <- is.na(x) | tolower(x) %in% c("", ".", "-", "na", "n/a", "nan", "null", "none")
  x[missing] <- NA_character_
  x
}

.echogo_is_portable_name <- function(x) {
  x <- .echogo_clean_gene_id(x)
  ok <- !is.na(x) & nchar(x) <= 64L &
    grepl("^[A-Za-z][A-Za-z0-9.-]*$", x) &
    grepl("[A-Za-z]", x)
  raw_id <- grepl(
    paste(
      c(
        "^TRINITY([_.-]|$)", "^(tri|tx|transcript|contig|scaffold|NODE)[_.-]?[0-9]",
        "^[0-9]+\\.", "^ENS[A-Z]*[GPT][0-9]", "^[NXWY][MRP]_[0-9]",
        "^LOC[0-9]+$", "^[A-Z]{1,4}[0-9]{5,}\\.[0-9]+$"
      ),
      collapse = "|"
    ),
    x,
    ignore.case = TRUE
  )
  ok & !raw_id
}

.echogo_is_uniprot_accession <- function(x) {
  x <- toupper(.echogo_clean_gene_id(x))
  !is.na(x) & (
    grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$", x) |
      grepl("^A0A[A-Z0-9]{7}$", x)
  )
}

.echogo_extract_swissprot_symbol <- function(x) {
  raw <- .echogo_clean_gene_id(x)
  vapply(raw, function(value) {
    if (is.na(value)) return(NA_character_)
    first <- sub("\\^.*$", "", value)
    parts <- strsplit(first, "|", fixed = TRUE)[[1]]
    if (length(parts) >= 3L && identical(tolower(parts[[1]]), "tr")) {
      return(NA_character_)
    }
    entry <- if (length(parts) >= 3L && identical(tolower(parts[[1]]), "sp")) {
      parts[[3]]
    } else {
      parts[[1]]
    }
    has_species_suffix <- grepl("_[A-Za-z0-9]{2,15}$", entry)
    if (!has_species_suffix && .echogo_is_uniprot_accession(entry)) {
      return(NA_character_)
    }
    symbol <- sub("_[A-Za-z0-9]{2,15}$", "", entry)
    if (.echogo_is_portable_name(symbol)) symbol else NA_character_
  }, character(1), USE.NAMES = FALSE)
}

.echogo_extract_swissprot_accession <- function(x) {
  raw <- .echogo_clean_gene_id(x)
  parts <- strsplit(sub("\\^.*$", "", raw), "|", fixed = TRUE)
  out <- vapply(parts, function(value) {
    if (length(value) >= 3L && tolower(value[[1]]) %in% c("sp", "tr")) value[[2]] else NA_character_
  }, character(1))
  .echogo_clean_gene_id(out)
}

.echogo_pick_gene_column <- function(df, explicit, candidates, label) {
  if (!is.null(explicit)) {
    if (!explicit %in% names(df)) {
      stop("Column '", explicit, "' was not found for ", label, ".", call. = FALSE)
    }
    return(explicit)
  }
  normalize <- function(x) tolower(gsub("[^a-z0-9]", "", x))
  hit <- match(normalize(candidates), normalize(names(df)))
  hit <- hit[!is.na(hit)]
  if (!length(hit)) {
    stop(
      "Could not identify the ", label, " column. Available columns: ",
      paste(names(df), collapse = ", "),
      call. = FALSE
    )
  }
  names(df)[hit[[1]]]
}

.echogo_parse_significance <- function(x, column) {
  if (is.logical(x)) return(dplyr::coalesce(x, FALSE))
  values <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(values))
  out[is.na(values) | values %in% c("", "na", "nan")] <- FALSE
  out[values %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[values %in% c("false", "f", "0", "no", "n")] <- FALSE
  if (anyNA(out)) {
    bad <- unique(values[is.na(out)])
    stop(
      "Significance column '", column, "' contains an ambiguous value: ",
      paste(utils::head(bad, 5L), collapse = ", "),
      ". Use logical TRUE/FALSE or documented 0/1 values.",
      call. = FALSE
    )
  }
  out
}

.echogo_read_tested_ids <- function(x, id_column = NULL) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) {
    column <- if (is.null(id_column)) names(x)[[1]] else id_column
    if (!column %in% names(x)) stop("Tested-gene ID column not found: ", column, call. = FALSE)
    return(.echogo_clean_gene_id(x[[column]]))
  }
  if (is.character(x) && length(x) == 1L && file.exists(x)) {
    table <- tryCatch(
      .echogo_read_delim_robust(x, candidates = c("\t", ",", ";")),
      error = function(e) NULL
    )
    if (!is.null(table) && ncol(table)) {
      column <- if (is.null(id_column)) names(table)[[1]] else id_column
      if (!column %in% names(table)) stop("Tested-gene ID column not found: ", column, call. = FALSE)
      return(.echogo_clean_gene_id(table[[column]]))
    }
    return(.echogo_clean_gene_id(readLines(x, warn = FALSE)))
  }
  .echogo_clean_gene_id(x)
}

#' Prepare matched portable gene sets for g:Profiler
#'
#' Derives a significant foreground and a tested experiment background, then
#' resolves both through one deterministic canonical-name mapping. The returned
#' vectors are shared across requested organism contexts; this function does not
#' construct species-specific ortholog vectors.
#'
#' @param de_results DE table or path.
#' @param tested_gene_ids Optional tested-gene vector, table, or file. This has
#'   precedence over a full DE table and count matrix.
#' @param annotation Annotation table or path.
#' @param count_matrix Optional matching count matrix or path.
#' @param de_id_column,significant_column,padj_column,lfc_column Optional column names.
#' @param tested_id_column,count_id_column,annotation_id_column Optional ID columns.
#' @param preferred_name_columns Legacy annotation-name aliases. Resolution is
#'   always semantic: SwissProt gene symbol, eggNOG preferred name, then a
#'   portable native symbol. Raw transcript/contig/accession IDs are excluded.
#' @param padj_threshold,log2fc_threshold Default DE significance thresholds.
#' @param de_table_significant_only Declare that every DE row is significant.
#' @param use_trinotate_universe Deprecated compatibility flag. It never permits
#'   raw-ID fallback or changes the tested universe. Portable original IDs may be
#'   used as native symbols; non-portable IDs remain provenance-only exclusions.
#' @return A structured list containing mapping audit tables and submitted vectors.
#' @keywords internal
prepare_gprofiler_gene_sets <- function(
    de_results,
    tested_gene_ids = NULL,
    annotation,
    count_matrix = NULL,
    de_id_column = NULL,
    significant_column = NULL,
    padj_column = NULL,
    lfc_column = NULL,
    tested_id_column = NULL,
    count_id_column = NULL,
    annotation_id_column = NULL,
    preferred_name_columns = c(
      "sprot_Top_BLASTX_hit", "EggNM.Preferred_name", "Preferred_name", "primary_name"
    ),
    padj_threshold = 0.05,
    log2fc_threshold = 1,
    de_table_significant_only = FALSE,
    use_trinotate_universe = FALSE
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  if (!is.numeric(padj_threshold) || length(padj_threshold) != 1L ||
      is.na(padj_threshold) || padj_threshold <= 0 || padj_threshold > 1) {
    stop("padj_threshold must be a number in (0, 1].", call. = FALSE)
  }
  if (!is.numeric(log2fc_threshold) || length(log2fc_threshold) != 1L ||
      is.na(log2fc_threshold) || log2fc_threshold < 0) {
    stop("log2fc_threshold must be a non-negative number.", call. = FALSE)
  }

  de <- .echogo_gene_set_table(de_results, "de_results")
  ann <- .echogo_gene_set_table(annotation, "annotation")
  de_id_column <- .echogo_pick_gene_column(
    de, de_id_column,
    c("gene_id", "ID", "gene", "GeneID", "transcript_id", "symbol"),
    "DE gene ID"
  )
  de_ids <- .echogo_clean_gene_id(de[[de_id_column]])

  detected_sig <- significant_column
  if (is.null(detected_sig)) {
    normalized <- tolower(gsub("[^a-z0-9]", "", names(de)))
    position <- match(c("significant", "isde", "issignificant"), normalized)
    position <- position[!is.na(position)]
    if (length(position)) detected_sig <- names(de)[position[[1]]]
  }

  if (!is.null(detected_sig)) {
    if (!detected_sig %in% names(de)) {
      stop("Significance column not found: ", detected_sig, call. = FALSE)
    }
    significant <- .echogo_parse_significance(de[[detected_sig]], detected_sig)
    significance_rule <- list(
      type = "logical_column",
      column = detected_sig,
      padj_threshold = NULL,
      log2fc_threshold = NULL,
      missing_adjusted_p = "not_applicable"
    )
  } else if (isTRUE(de_table_significant_only)) {
    significant <- rep(TRUE, nrow(de))
    significance_rule <- list(
      type = "declared_significant_only_table",
      column = NULL,
      padj_threshold = NULL,
      log2fc_threshold = NULL,
      missing_adjusted_p = "not_applicable"
    )
  } else {
    padj_column <- .echogo_pick_gene_column(
      de, padj_column,
      c("padj", "FDR", "adj.P.Val", "adjusted_pvalue", "qvalue"),
      "adjusted p-value"
    )
    lfc_column <- .echogo_pick_gene_column(
      de, lfc_column,
      c("log2FoldChange", "log2FC", "LFC", "log2FoldChg"),
      "log2 fold-change"
    )
    padj <- suppressWarnings(as.numeric(de[[padj_column]]))
    lfc <- suppressWarnings(as.numeric(de[[lfc_column]]))
    significant <- !is.na(padj) & padj <= padj_threshold &
      !is.na(lfc) & abs(lfc) >= log2fc_threshold
    significance_rule <- list(
      type = "threshold",
      padj_column = padj_column,
      lfc_column = lfc_column,
      padj_threshold = padj_threshold,
      log2fc_threshold = log2fc_threshold,
      missing_adjusted_p = "excluded"
    )
  }

  valid_de <- !is.na(de_ids)
  de_status <- data.frame(
    original_id = de_ids[valid_de],
    significant = significant[valid_de],
    row_order = which(valid_de),
    stringsAsFactors = FALSE
  )
  de_status <- dplyr::summarise(
    dplyr::group_by(de_status, .data$original_id),
    significant = any(.data$significant, na.rm = TRUE),
    row_order = min(.data$row_order),
    .groups = "drop"
  )
  de_status <- dplyr::arrange(de_status, .data$row_order)
  foreground_original <- de_status$original_id[de_status$significant]
  if (!length(foreground_original)) {
    stop("The DE significance rule selected zero genes; no g:Profiler foreground can be built.", call. = FALSE)
  }

  tested_source <- NULL
  background_original <- .echogo_read_tested_ids(tested_gene_ids, tested_id_column)
  if (!is.null(background_original)) {
    tested_source <- "explicit_tested_gene_universe"
  } else if (!isTRUE(de_table_significant_only)) {
    background_original <- de_status$original_id
    tested_source <- "full_de_table"
  } else if (!is.null(count_matrix)) {
    counts <- .echogo_gene_set_table(count_matrix, "count_matrix")
    count_id_column <- count_id_column %||% names(counts)[[1]]
    if (!count_id_column %in% names(counts)) {
      stop("Count-matrix ID column not found: ", count_id_column, call. = FALSE)
    }
    background_original <- .echogo_clean_gene_id(counts[[count_id_column]])
    tested_source <- "count_matrix"
  } else {
    stop(
      "No tested-gene universe is available. Supply tested_gene_ids, use a full DE table, ",
      "or provide count_matrix when de_table_significant_only = TRUE.",
      call. = FALSE
    )
  }
  background_original <- unique(background_original[!is.na(background_original)])
  if (!length(background_original)) {
    stop("The tested-gene universe contains no usable identifiers.", call. = FALSE)
  }

  missing_tested <- setdiff(foreground_original, background_original)
  if (length(missing_tested)) {
    stop(
      length(missing_tested), " significant gene(s) are not present in the tested universe: ",
      paste(utils::head(missing_tested, 8L), collapse = ", "),
      ". Check that the DE contrast and tested/background input describe the same experiment.",
      call. = FALSE
    )
  }

  annotation_id_column <- .echogo_pick_gene_column(
    ann, annotation_id_column,
    c("transcript_id", "gene_id", "#query", "query"),
    "annotation ID"
  )
  ann$.echogo_id <- .echogo_clean_gene_id(ann[[annotation_id_column]])
  all_original <- unique(c(background_original, foreground_original))
  resolved <- rep(NA_character_, length(all_original))
  resolution_source <- rep(NA_character_, length(all_original))
  valid_annotation_ids <- ann$.echogo_id[!is.na(ann$.echogo_id)]
  annotation_match <- all_original %in% valid_annotation_ids

  stable_lookup <- function(values) {
    values <- .echogo_clean_gene_id(values)
    usable <- !is.na(ann$.echogo_id) & !is.na(values)
    if (!any(usable)) return(stats::setNames(character(), character()))
    ids <- ann$.echogo_id[usable]
    values <- values[usable]
    keep <- !duplicated(ids)
    stats::setNames(values[keep], ids[keep])
  }
  lookup_for <- function(values) {
    lookup <- stable_lookup(values)
    unname(lookup[match(all_original, names(lookup))])
  }
  first_from_columns <- function(columns, transform = identity, portable_only = FALSE) {
    out <- rep(NA_character_, length(all_original))
    for (column in intersect(columns, names(ann))) {
      candidate <- lookup_for(transform(ann[[column]]))
      if (portable_only) candidate[!.echogo_is_portable_name(candidate)] <- NA_character_
      take <- is.na(out) & !is.na(candidate)
      out[take] <- candidate[take]
    }
    out
  }
  first_raw_from_columns <- function(columns) {
    out <- rep(NA_character_, length(all_original))
    for (column in intersect(columns, names(ann))) {
      values <- as.character(ann[[column]])
      usable <- !is.na(ann$.echogo_id) & !is.na(values) & nzchar(values)
      if (!any(usable)) next
      ids <- ann$.echogo_id[usable]
      values <- values[usable]
      keep <- !duplicated(ids)
      lookup <- stats::setNames(values[keep], ids[keep])
      candidate <- unname(lookup[match(all_original, names(lookup))])
      take <- is.na(out) & !is.na(candidate)
      out[take] <- candidate[take]
    }
    out
  }

  swiss_columns <- unique(c(
    "sprot_Top_BLASTX_hit", "sprot_Top_BLASTP_hit",
    preferred_name_columns[grepl("sprot|swiss", preferred_name_columns, ignore.case = TRUE)]
  ))
  eggnog_columns <- unique(c(
    "EggNM.Preferred_name", "Preferred_name", "eggnog_preferred_name",
    preferred_name_columns[grepl("preferred", preferred_name_columns, ignore.case = TRUE)]
  ))
  native_columns <- unique(c(
    "native_symbol", "gene_symbol", "symbol", "primary_name",
    preferred_name_columns[grepl("native|symbol|primary", preferred_name_columns, ignore.case = TRUE)]
  ))

  swissprot_gene_symbol <- first_from_columns(
    swiss_columns, .echogo_extract_swissprot_symbol, portable_only = TRUE
  )
  swissprot_accession <- first_from_columns(
    swiss_columns, .echogo_extract_swissprot_accession
  )
  eggnog_preferred_name <- first_from_columns(
    eggnog_columns, .echogo_clean_gene_id, portable_only = TRUE
  )
  native_symbol <- first_from_columns(
    native_columns, .echogo_clean_gene_id, portable_only = TRUE
  )

  # A portable annotation or DE identifier may serve as a native symbol. The
  # same validation applies to annotation-matched and unmatched identifiers.
  annotation_native <- lookup_for(.echogo_clean_gene_id(ann$.echogo_id))
  annotation_native[!.echogo_is_portable_name(annotation_native)] <- NA_character_
  take_native <- is.na(native_symbol) & !is.na(annotation_native)
  native_symbol[take_native] <- annotation_native[take_native]
  original_native <- all_original
  original_native[!.echogo_is_portable_name(original_native)] <- NA_character_
  take_original <- is.na(native_symbol) & !is.na(original_native)
  native_symbol[take_original] <- original_native[take_original]

  resolution_candidates <- list(
    swissprot_gene_symbol = swissprot_gene_symbol,
    eggnog_preferred_name = eggnog_preferred_name,
    native_symbol = native_symbol
  )
  for (source in names(resolution_candidates)) {
    candidate <- resolution_candidates[[source]]
    take <- is.na(resolved) & !is.na(candidate)
    resolved[take] <- candidate[take]
    resolution_source[take] <- source
  }

  taxonomy_columns <- c(
    "annotation_taxonomy", "EggNM.max_annot_lvl", "max_annot_lvl",
    "blast_taxonomy"
  )
  annotation_taxonomy <- first_from_columns(taxonomy_columns, .echogo_clean_gene_id)
  if (all(is.na(annotation_taxonomy)) && "sprot_Top_BLASTX_hit" %in% names(ann)) {
    blast_taxonomy <- sub("^[^^]*\\^", "", .echogo_clean_gene_id(ann$sprot_Top_BLASTX_hit))
    blast_taxonomy[blast_taxonomy == ann$sprot_Top_BLASTX_hit] <- NA_character_
    annotation_taxonomy <- lookup_for(blast_taxonomy)
  }

  mapping <- data.frame(
    original_gene_id = all_original,
    original_id = all_original,
    tested = all_original %in% background_original,
    significant = all_original %in% foreground_original,
    portable_name = resolved,
    resolved_name = resolved,
    portable_name_source = resolution_source,
    resolution_source = resolution_source,
    native_symbol = native_symbol,
    swissprot_accession = swissprot_accession,
    swissprot_gene_symbol = swissprot_gene_symbol,
    eggnog_preferred_name = eggnog_preferred_name,
    annotation_taxonomy = annotation_taxonomy,
    annotation_match = annotation_match,
    stringsAsFactors = FALSE
  )
  # Preserve richer annotation fields when supplied.  These are descriptive
  # provenance only: none participates in resolution, enrichment, recurrence,
  # or a confidence/ranking calculation.
  optional_provenance <- list(
    eggnog_seed_ortholog = c("EggNM.seed_ortholog", "seed_ortholog", "seed_ortholog_id"),
    eggnog_seed_evalue = c("EggNM.seed_evalue", "seed_evalue", "seed_evalue_score"),
    eggnog_seed_score = c("EggNM.seed_score", "EggNM.seed_ortholog_score", "seed_score"),
    eggnog_ogs = c("EggNM.OGs", "eggnog_ogs", "eggNOG_OGs", "OGs"),
    annotation_go_field = c("EggNM.GOs", "GO", "go_terms", "annotation_go_field"),
    eggnog_mapper_version = c("eggnog_mapper_version", "eggNOG_mapper_version"),
    go_evidence_filter_policy = c("go_evidence_filter_policy", "GO_evidence_filter_policy")
  )
  for (field in names(optional_provenance)) {
    source_columns <- intersect(optional_provenance[[field]], names(ann))
    if (length(source_columns)) mapping[[field]] <- first_raw_from_columns(source_columns)
  }

  bg_candidates <- which(mapping$tested & !is.na(mapping$resolved_name))
  fg_candidates <- which(mapping$tested & mapping$significant & !is.na(mapping$resolved_name))
  bg_keep <- bg_candidates[!duplicated(mapping$resolved_name[bg_candidates])]
  fg_keep <- fg_candidates[!duplicated(mapping$resolved_name[fg_candidates])]
  mapping$included_background <- seq_len(nrow(mapping)) %in% bg_keep
  mapping$included_foreground <- seq_len(nrow(mapping)) %in% fg_keep
  mapping$included_in_gprofiler <- mapping$included_background
  duplicated_names <- unique(mapping$resolved_name[
    !is.na(mapping$resolved_name) & duplicated(mapping$resolved_name)
  ])
  mapping$duplicate_group <- ifelse(
    !is.na(mapping$resolved_name) & mapping$resolved_name %in% duplicated_names,
    mapping$resolved_name,
    NA_character_
  )
  mapping$exclusion_reason <- NA_character_
  mapping$exclusion_reason[is.na(mapping$resolved_name)] <- "no portable canonical name"
  mapping$exclusion_reason[!mapping$tested & mapping$significant] <- "not_in_tested_universe"
  duplicate_only <- !mapping$included_background & mapping$tested &
    !is.na(mapping$resolved_name) & mapping$resolved_name %in% duplicated_names
  mapping$exclusion_reason[duplicate_only] <- "duplicate_resolved_name"

  foreground_resolved <- mapping$resolved_name[mapping$included_foreground]
  background_resolved <- mapping$resolved_name[mapping$included_background]
  outside <- setdiff(foreground_resolved, background_resolved)
  if (length(outside)) {
    stop(
      "Resolved foreground is not a subset of the resolved background; ",
      length(outside), " offending identifier(s): ",
      paste(utils::head(outside, 8L), collapse = ", "),
      ". Check resolver inputs and the tested universe.",
      call. = FALSE
    )
  }
  if (!length(background_resolved)) {
    stop("No tested genes could be resolved for a custom g:Profiler background.", call. = FALSE)
  }
  if (!length(foreground_resolved)) {
    stop("No significant genes could be resolved for the g:Profiler foreground.", call. = FALSE)
  }
  if (identical(foreground_resolved, background_resolved)) {
    stop(
      "The foreground and background resolve to the same vector. ",
      "Provide the complete tested universe; EchoGO will not submit a degenerate custom background.",
      call. = FALSE
    )
  }

  list(
    mapping_table = tibble::as_tibble(mapping),
    foreground_original = foreground_original,
    background_original = background_original,
    foreground_resolved = foreground_resolved,
    background_resolved = background_resolved,
    foreground_unmapped = mapping$original_id[mapping$significant & is.na(mapping$resolved_name)],
    background_unmapped = mapping$original_id[mapping$tested & is.na(mapping$resolved_name)],
    duplicates_collapsed = list(
      foreground = as.integer(length(fg_candidates) - length(foreground_resolved)),
      background = as.integer(length(bg_candidates) - length(background_resolved))
    ),
    significance_rule = significance_rule,
    tested_universe_source = tested_source,
    resolver_definition = list(
      route = "shared_portable_canonical_organism_context",
      annotation_id_column = annotation_id_column,
      priority = c("SwissProt gene symbol", "eggNOG preferred name", "portable native symbol"),
      swissprot_columns = intersect(swiss_columns, names(ann)),
      eggnog_preferred_columns = intersect(eggnog_columns, names(ann)),
      native_symbol_columns = intersect(native_columns, names(ann)),
      swissprot_transform = "Trinotate SwissProt hit to gene symbol; accessions and species suffixes are not submitted",
      fallback = "none; raw target IDs are retained only in mapping provenance",
      portability_filter = "symbol-like names only; transcript, contig, seed-ortholog, Ensembl protein and arbitrary accessions excluded",
      taxonomy_policy = "taxonomy-independent resolver; annotation taxonomy retained as provenance only",
      duplicate_policy = "stable first portable canonical name",
      use_trinotate_universe = isTRUE(use_trinotate_universe),
      use_trinotate_universe_effect = "compatibility metadata only; does not redefine the tested universe",
      explicit_species_specific_ortholog_mapping = FALSE
    )
  )
}

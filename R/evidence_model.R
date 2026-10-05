# Scoreless evidence-model helpers -------------------------------------------

.echogo_normalize_context_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- NA_character_
  x <- trimws(gsub("[[:space:]]+", " ", x))
  x[!nzchar(x)] <- NA_character_
  x
}

.echogo_context_configuration <- function(species_map, target_context = NULL,
                                          context_metadata = NULL) {
  if (is.null(names(species_map)) || !any(nzchar(names(species_map)))) {
    species_map <- stats::setNames(as.character(species_map), as.character(species_map))
  }
  codes <- as.character(names(species_map))
  labels <- as.character(unname(species_map))
  if (anyDuplicated(labels)) {
    stop("Each researcher-selected context must have a unique label.", call. = FALSE)
  }
  target_label <- NULL
  if (!is.null(target_context)) {
    target_context <- as.character(target_context)
    if (length(target_context) != 1L || is.na(target_context) || !nzchar(target_context)) {
      stop("target_context must be NULL or one queried g:Profiler context.", call. = FALSE)
    }
    index <- match(target_context, codes)
    if (is.na(index)) index <- match(target_context, labels)
    if (is.na(index)) {
      stop(
        "target_context ('", target_context,
        "') is not among the queried contexts: ", paste(codes, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    target_label <- labels[[index]]
  }
  contexts <- tibble::tibble(
    context_code = codes,
    context_label = labels,
    context_role = if (is.null(target_label)) {
      rep("ALTERNATIVE", length(labels))
    } else {
      ifelse(labels == target_label, "TARGET", "ALTERNATIVE")
    },
    context_rationale = NA_character_,
    selection_basis = NA_character_
  )

  if (is.null(context_metadata)) return(contexts)
  if (!is.data.frame(context_metadata)) {
    stop("context_metadata must be NULL or a data frame.", call. = FALSE)
  }
  allowed <- c("context_code", "context_label", "context_rationale", "selection_basis")
  supplied <- intersect(allowed, names(context_metadata))
  if (!length(intersect(c("context_code", "context_label"), supplied))) {
    stop(
      "context_metadata must identify each row with context_code and/or context_label.",
      call. = FALSE
    )
  }
  if (!nrow(context_metadata)) return(contexts)

  metadata <- tibble::as_tibble(context_metadata)
  for (field in allowed) {
    if (!field %in% names(metadata)) metadata[[field]] <- NA_character_
    metadata[[field]] <- .echogo_normalize_context_text(metadata[[field]])
  }
  if (any(is.na(metadata$context_code) & is.na(metadata$context_label))) {
    stop(
      "Every context_metadata row must identify a queried context with context_code and/or context_label.",
      call. = FALSE
    )
  }

  metadata$context_index <- vapply(seq_len(nrow(metadata)), function(i) {
    code_index <- if (!is.na(metadata$context_code[[i]])) {
      match(metadata$context_code[[i]], contexts$context_code)
    } else {
      NA_integer_
    }
    label_index <- if (!is.na(metadata$context_label[[i]])) {
      match(metadata$context_label[[i]], contexts$context_label)
    } else {
      NA_integer_
    }
    if (!is.na(code_index) && !is.na(label_index) && code_index != label_index) {
      stop(
        "context_metadata context_code and context_label must identify the same queried context.",
        call. = FALSE
      )
    }
    index <- if (!is.na(code_index)) code_index else label_index
    if (is.na(index)) {
      supplied_id <- if (!is.na(metadata$context_code[[i]])) {
        paste0("context_code '", metadata$context_code[[i]], "'")
      } else {
        paste0("context_label '", metadata$context_label[[i]], "'")
      }
      stop(
        "context_metadata references an unqueried context (", supplied_id,
        "). Queried context codes are: ", paste(contexts$context_code, collapse = ", "), ".",
        call. = FALSE
      )
    }
    index
  }, integer(1))
  if (anyDuplicated(metadata$context_index)) {
    stop("context_metadata must contain at most one row per queried context.", call. = FALSE)
  }
  for (field in c("context_rationale", "selection_basis")) {
    contexts[[field]][metadata$context_index] <- metadata[[field]]
  }
  contexts
}

.echogo_gprofiler_mode_paths <- function(gprofiler_dir, mode) {
  if (identical(mode, "custom_experimental_background")) {
    c(
      file.path(gprofiler_dir, "custom_experimental_background"),
      file.path(gprofiler_dir, "bg"),
      file.path(gprofiler_dir, "with_custom_background")
    )
  } else {
    c(
      file.path(gprofiler_dir, "default_domain_exploratory"),
      file.path(gprofiler_dir, "nobg"),
      file.path(gprofiler_dir, "no_background_genome_wide")
    )
  }
}

.echogo_read_json_scalar <- function(path, field) {
  if (!file.exists(path)) return(NA)
  value <- tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(value) || is.null(value[[field]])) return(NA)
  value[[field]]
}

.echogo_read_gprofiler_long <- function(gprofiler_dir, species_map,
                                        mode = c("custom_experimental_background", "default_domain_exploratory"),
                                        target_context = NULL,
                                        context_metadata = NULL) {
  mode <- match.arg(mode)
  contexts <- .echogo_context_configuration(species_map, target_context, context_metadata)
  source_dir <- .echogo_gprofiler_mode_paths(gprofiler_dir, mode)
  source_dir <- source_dir[dir.exists(source_dir)][1]
  if (!length(source_dir) || is.na(source_dir)) return(tibble::tibble())

  suffix <- if (identical(mode, "custom_experimental_background")) "with_bg" else "nobg"
  rows <- vector("list", nrow(contexts))
  for (i in seq_len(nrow(contexts))) {
    label <- contexts$context_label[[i]]
    stem <- file.path(source_dir, paste0("gprofiler_", label, "_", suffix))
    csv <- paste0(stem, ".csv")
    if (!file.exists(csv)) next
    tab <- readr::read_csv(csv, show_col_types = FALSE)
    if (!all(c("term_id", "p_value") %in% names(tab))) {
      stop("g:Profiler CSV lacks term_id or p_value: ", csv, call. = FALSE)
    }
    metadata <- paste0(stem, "_metadata.json")
    scalar <- function(field, fallback = NA) {
      value <- .echogo_read_json_scalar(metadata, field)
      if (length(value) == 0L || all(is.na(value))) fallback else value[[1]]
    }
    rows[[i]] <- tibble::tibble(
      term_id = as.character(tab$term_id),
      term_name = if ("term_name" %in% names(tab)) as.character(tab$term_name) else NA_character_,
      ontology = if ("source" %in% names(tab)) as.character(tab$source) else NA_character_,
      source_type = "gprofiler",
      context_code = contexts$context_code[[i]],
      context_label = label,
      context_role = contexts$context_role[[i]],
      context_rationale = contexts$context_rationale[[i]],
      selection_basis = contexts$selection_basis[[i]],
      background_mode = mode,
      domain_scope = if (identical(mode, "custom_experimental_background")) "custom" else "annotated",
      source_adjusted_p = suppressWarnings(as.numeric(tab$p_value)),
      fold_enrichment = if ("fold_enrichment" %in% names(tab)) suppressWarnings(as.numeric(tab$fold_enrichment)) else NA_real_,
      contributing_genes = if ("intersection" %in% names(tab)) as.character(tab$intersection) else NA_character_,
      intersection_size = if ("intersection_size" %in% names(tab)) suppressWarnings(as.numeric(tab$intersection_size)) else NA_real_,
      query_size = if ("query_size" %in% names(tab)) suppressWarnings(as.numeric(tab$query_size)) else NA_real_,
      effective_domain_size = if ("effective_domain_size" %in% names(tab)) suppressWarnings(as.numeric(tab$effective_domain_size)) else NA_real_,
      submitted_foreground_n = suppressWarnings(as.numeric(scalar("submitted_foreground_count"))),
      submitted_background_n = suppressWarnings(as.numeric(scalar("submitted_background_count"))),
      effective_query_n = suppressWarnings(as.numeric(scalar("effective_query_size"))),
      effective_background_or_domain_n = suppressWarnings(as.numeric(scalar("effective_domain_size"))),
      provenance_file = normalizePath(csv, winslash = "/", mustWork = FALSE),
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

.echogo_read_goseq_long <- function(goseq_file) {
  tab <- .echogo_read_delim_robust(
    goseq_file,
    candidates = c("\t", ",", ";"),
    expected = c("clean_go_term", "category", "term_id", "GO.ID", "go_id")
  )
  id_col <- intersect(c("clean_go_term", "category", "term_id"), names(tab))[1]
  p_col <- intersect(c("over_represented_FDR", "FDR", "padj", "qvalue"), names(tab))[1]
  if (is.na(id_col) || is.na(p_col)) {
    stop("GOseq output must contain a term ID and adjusted-p column.", call. = FALSE)
  }
  tibble::tibble(
    term_id = as.character(tab[[id_col]]),
    term_name = if ("term" %in% names(tab)) as.character(tab$term) else as.character(tab[[id_col]]),
    ontology = if ("ontology" %in% names(tab)) as.character(tab$ontology) else NA_character_,
    source_type = "goseq",
    context_code = NA_character_,
    context_label = "target_goseq",
    context_role = "TARGET_GOSEQ",
    background_mode = "target_annotation_goseq",
    domain_scope = "target_goseq_tested_universe",
    source_adjusted_p = suppressWarnings(as.numeric(tab[[p_col]])),
    fold_enrichment = if ("foldEnrichment" %in% names(tab)) suppressWarnings(as.numeric(tab$foldEnrichment)) else NA_real_,
    contributing_genes = if ("gene_names" %in% names(tab)) as.character(tab$gene_names) else if ("gene_ids" %in% names(tab)) as.character(tab$gene_ids) else NA_character_,
    intersection_size = if ("numDEInCat" %in% names(tab)) suppressWarnings(as.numeric(tab$numDEInCat)) else NA_real_,
    query_size = if ("total_significant_genes" %in% names(tab)) suppressWarnings(as.numeric(tab$total_significant_genes)) else NA_real_,
    effective_domain_size = if ("total_tested_genes" %in% names(tab)) suppressWarnings(as.numeric(tab$total_tested_genes)) else NA_real_,
    submitted_foreground_n = NA_real_,
    submitted_background_n = NA_real_,
    effective_query_n = NA_real_,
    effective_background_or_domain_n = NA_real_,
    provenance_file = normalizePath(goseq_file, winslash = "/", mustWork = FALSE),
    stringsAsFactors = FALSE
  )
}

.echogo_profile_order <- function(x) {
  match(x, c("TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT", "TARGET_ONLY", "NO_PRIMARY_SUPPORT"), nomatch = 99L)
}

.echogo_add_display_order <- function(evidence) {
  if (!nrow(evidence)) return(evidence)
  ordering <- order(
    .echogo_profile_order(evidence$evidence_profile),
    -evidence$alternative_context_support_n,
    evidence$term_id,
    na.last = TRUE
  )
  evidence$display_order <- NA_integer_
  evidence$display_order[ordering] <- seq_along(ordering)
  evidence$display_rank <- evidence$display_order
  evidence$representative_order <- evidence$display_order
  evidence
}

#' Assemble scoreless exact-term evidence and source provenance
#'
#' @param goseq_long Long-form GOseq provenance.
#' @param gprofiler_custom_long Long-form custom-background g:Profiler provenance.
#' @param gprofiler_exploratory_long Long-form default-domain provenance.
#' @param context_configuration Researcher-selected context roles.
#' @param adjusted_p_threshold Source-specific qualification threshold.
#' @return A list with `exact_terms` and `source_provenance` tables.
#' @keywords internal
assemble_echogo_evidence <- function(goseq_long,
                                     gprofiler_custom_long = tibble::tibble(),
                                     gprofiler_exploratory_long = tibble::tibble(),
                                     context_configuration,
                                     adjusted_p_threshold = 0.05) {
  if (!is.numeric(adjusted_p_threshold) || length(adjusted_p_threshold) != 1L ||
      is.na(adjusted_p_threshold) || adjusted_p_threshold <= 0 || adjusted_p_threshold > 1) {
    stop("adjusted_p_threshold must be in (0, 1].", call. = FALSE)
  }
  required <- c("term_id", "term_name", "ontology", "source_type", "context_role",
                "background_mode", "source_adjusted_p", "fold_enrichment", "contributing_genes")
  provenance <- dplyr::bind_rows(goseq_long, gprofiler_custom_long, gprofiler_exploratory_long)
  missing <- setdiff(required, names(provenance))
  if (length(missing)) stop("Source provenance lacks: ", paste(missing, collapse = ", "), call. = FALSE)
  provenance <- provenance %>%
    dplyr::filter(!is.na(.data$term_id), nzchar(.data$term_id)) %>%
    dplyr::mutate(
      source_adjusted_p = suppressWarnings(as.numeric(.data$source_adjusted_p)),
      source_qualifies = !is.na(.data$source_adjusted_p) & .data$source_adjusted_p <= adjusted_p_threshold,
      provenance_reference = paste0(.data$source_type, ":", dplyr::coalesce(.data$context_label, "target"), ":", .data$background_mode, ":", .data$term_id)
    )
  terms <- sort(unique(provenance$term_id))
  alternatives <- context_configuration$context_label[context_configuration$context_role == "ALTERNATIVE"]
  queried <- context_configuration$context_label

  exact <- lapply(terms, function(id) {
    d <- provenance[provenance$term_id == id, , drop = FALSE]
    gs <- d[d$source_type == "goseq", , drop = FALSE]
    custom <- d[d$background_mode == "custom_experimental_background", , drop = FALSE]
    exploratory <- d[d$background_mode == "default_domain_exploratory", , drop = FALSE]
    target_supported <- any(gs$source_qualifies, na.rm = TRUE)
    custom_supported <- unique(custom$context_label[custom$source_qualifies])
    alternative_supported <- intersect(custom_supported, alternatives)
    target_gprof_supported <- any(custom$context_role == "TARGET" & custom$source_qualifies, na.rm = TRUE)
    alternative_n <- length(alternative_supported)
    profile <- if (target_supported && alternative_n == 0L) {
      "TARGET_ONLY"
    } else if (target_supported && alternative_n > 0L) {
      "TARGET_PLUS_CONTEXT"
    } else if (!target_supported && alternative_n > 0L) {
      "ALTERNATIVE_CONTEXT"
    } else {
      "NO_PRIMARY_SUPPORT"
    }
    first_non_missing <- function(x) { x <- x[!is.na(x) & nzchar(x)]; if (length(x)) x[[1]] else NA_character_ }
    first_source_value <- function(x) {
      x <- x[!is.na(x)]
      if (length(x)) x[[1]] else NA_real_
    }
    tibble::tibble(
      term_id = id,
      term_name = first_non_missing(d$term_name),
      ontology = first_non_missing(d$ontology),
      target_goseq_supported = target_supported,
      # A GOseq adjusted p-value is retained from its own source row only.
      # No best/minimum value is selected across contexts or analyses.
      target_goseq_adjusted_p = first_source_value(gs$source_adjusted_p),
      target_gprof_context_supported = target_gprof_supported,
      gprof_custom_context_support_n = length(custom_supported),
      gprof_custom_context_support_fraction = if (length(queried)) length(custom_supported) / length(queried) else NA_real_,
      alternative_context_support_n = alternative_n,
      alternative_context_support_fraction = if (length(alternatives)) alternative_n / length(alternatives) else NA_real_,
      queried_context_n = length(queried),
      alternative_queried_context_n = length(alternatives),
      recovered_custom_contexts = paste(custom_supported, collapse = ";"),
      recovered_alternative_contexts = paste(alternative_supported, collapse = ";"),
      default_domain_exploratory_context_support_n = length(unique(exploratory$context_label[exploratory$source_qualifies])),
      evidence_profile = profile,
      contributing_genes = paste(unique(stats::na.omit(d$contributing_genes)), collapse = ";"),
      source_provenance_references = paste(unique(d$provenance_reference), collapse = ";"),
      source_provenance_n = nrow(d),
      stringsAsFactors = FALSE
    )
  })
  exact <- dplyr::bind_rows(exact)
  list(exact_terms = .echogo_add_display_order(exact), source_provenance = provenance)
}

#' Retrieve the traceable evidence for an alternative-context hypothesis
#'
#' @param evidence_table Exact-term scoreless evidence table.
#' @param source_provenance Long-form source provenance table.
#' @param term_id GO term identifier.
#' @param annotation_provenance Optional gene-level annotation provenance table.
#' @return A list containing the exact term evidence, relevant source rows, normalized
#'   context-selection provenance, and matching gene provenance.
#' @export
echogo_hypothesis_drilldown <- function(evidence_table, source_provenance, term_id,
                                        annotation_provenance = NULL) {
  term <- evidence_table[evidence_table$term_id == term_id, , drop = FALSE]
  if (nrow(term) != 1L) stop("term_id must identify exactly one evidence-table row.", call. = FALSE)
  if (!identical(term$evidence_profile[[1]], "ALTERNATIVE_CONTEXT")) {
    warning("Requested term is not an ALTERNATIVE_CONTEXT hypothesis; returning its full evidence trail.", call. = FALSE)
  }
  sources <- source_provenance[source_provenance$term_id == term_id, , drop = FALSE]
  genes <- unique(trimws(unlist(strsplit(
    paste(stats::na.omit(as.character(sources$contributing_genes)), collapse = ","), ",", fixed = TRUE
  ))))
  genes <- genes[nzchar(genes)]
  annotation <- NULL
  if (!is.null(annotation_provenance) && nrow(annotation_provenance)) {
    key <- intersect(c("resolved_name", "portable_name"), names(annotation_provenance))[1]
    if (!is.na(key)) annotation <- annotation_provenance[annotation_provenance[[key]] %in% genes, , drop = FALSE]
  }
  context_fields <- c("context_code", "context_label", "context_role",
                      "context_rationale", "selection_basis")
  for (field in setdiff(context_fields, names(sources))) sources[[field]] <- NA_character_
  context_selection <- sources[sources$source_type == "gprofiler", context_fields, drop = FALSE]
  context_selection <- unique(context_selection[order(
    context_selection$context_code, context_selection$context_label
  ), , drop = FALSE])
  list(
    term_evidence = term,
    source_provenance = sources,
    context_selection = context_selection,
    contributing_genes = genes,
    annotation_provenance = annotation
  )
}

.echogo_write_annotation_provenance <- function(mapping_table, evidence_table, output_dir) {
  if (is.null(mapping_table) || !nrow(mapping_table)) return(invisible(NULL))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  mapping <- as.data.frame(mapping_table, stringsAsFactors = FALSE)
  # Provenance is optional.  Keep all supplied resolver fields intact while
  # representing unsupported fields as explicitly unavailable rather than
  # requiring an eggNOG-annotated input.
  for (field in c("resolved_name", "swissprot_accession", "eggnog_preferred_name",
                  "annotation_taxonomy", "tested", "significant")) {
    if (!field %in% names(mapping)) mapping[[field]] <- NA
  }
  mapping$annotation_provenance_available <- !is.na(mapping$resolved_name) |
    !is.na(mapping$swissprot_accession) | !is.na(mapping$eggnog_preferred_name) |
    !is.na(mapping$annotation_taxonomy)
  readr::write_csv(mapping, file.path(output_dir, "annotation_provenance_gene_level.csv"))
  fields <- intersect(c("native_symbol", "swissprot_accession", "swissprot_gene_symbol",
                        "eggnog_preferred_name", "annotation_taxonomy", "annotation_match",
                        "resolution_source"), names(mapping))
  summary <- tibble::tibble(
    field = fields,
    available_n = vapply(fields, function(x) sum(!is.na(mapping[[x]]) & nzchar(as.character(mapping[[x]]))), numeric(1)),
    tested_n = sum(mapping$tested %in% TRUE),
    significant_n = sum(mapping$significant %in% TRUE)
  )
  readr::write_csv(summary, file.path(output_dir, "annotation_provenance_summary.csv"))
  term_summary <- evidence_table %>%
    dplyr::transmute(
      term_id, term_name, evidence_profile,
      target_goseq_supported, alternative_context_support_n,
      annotation_provenance_note = "Gene-level provenance is retained separately; unavailable fields remain missing."
    )
  readr::write_csv(term_summary, file.path(output_dir, "term_annotation_provenance_summary.csv"))
  invisible(list(mapping = mapping, summary = summary, term_summary = term_summary))
}

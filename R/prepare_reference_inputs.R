#' Prepare EchoGO inputs from a reference-based RNA-seq project
#'
#' Builds the reference-based EchoGO input bundle in two resumable stages.
#' The first stage validates counts and DESeq2 tables, parses an NCBI-style
#' GFF3, creates gene/protein and gene-length maps, and writes the background
#' protein FASTA for eggNOG-mapper. If the eggNOG-mapper annotations file is
#' present, the second stage creates the Trinotate-like annotation tables and
#' runs GOseq for every detected contrast.
#'
#' @param root Project directory containing the source files.
#' @param gff_file Genome annotation in GFF3 or gzip-compressed GFF3 format.
#' @param protein_fasta Protein FASTA matching the GFF3; gzip is supported.
#' @param counts_file Count matrix path or filename relative to `root`.
#' @param dge_pattern Regular expression selecting DESeq2 CSV files in `root`.
#' @param emapper_file eggNOG-mapper annotations path. The default is inside
#'   `work_dir/04_emapper`.
#' @param reference_label Short label used for the species-labelled annotation
#'   filename.
#' @param output_dir Final EchoGO input directory.
#' @param work_dir Directory for intermediate files.
#' @param padj_threshold Adjusted p-value threshold for DE genes.
#' @param log2fc_threshold Absolute log2 fold-change threshold for DE genes.
#' @param strip_gene_prefix Remove an NCBI-style `gene-` prefix from IDs.
#' @param counts_id_column Optional count-matrix gene ID column. Defaults to the
#'   first column.
#' @param dge_id_column,dge_lfc_column,dge_padj_column Optional DESeq2 column
#'   names. Common synonyms are detected when omitted.
#' @param dge_significant_column Optional logical/0-1 significance column.
#'   When omitted, a column named `significant` or `is_de` is used when found;
#'   otherwise `padj_threshold` and `log2fc_threshold` are applied.
#' @param annotation_taxid Deprecated compatibility fallback used only when the
#'   eggNOG output has no annotation-level field. It is never used as a numeric
#'   eligibility cutoff.
#' @param blast_taxonomy Taxonomy label paired with `annotation_taxid` for that
#'   provenance fallback. It is not written as a synthetic SwissProt hit.
#' @param orgdb Optional installed OrgDb package name, such as
#'   `"org.Dr.eg.db"`. When supplied, GOseq categories are taken from the
#'   species-native OrgDb using SYMBOL/ALIAS/ENTREZID matching; eggNOG remains
#'   the protein-name annotation bridge. When `NULL`, eggNOG GO terms are used.
#' @param overwrite Overwrite generated files when rerunning the helper.
#' @param reuse_reference Reuse existing gene maps and background FASTA in
#'   `work_dir`. Set to `FALSE` after changing the reference files or ID rules.
#' @param verbose Print stage and validation messages.
#'
#' @return Invisibly returns a list with `status` equal to
#'   `"awaiting_emapper"` or `"complete"`, plus generated paths and summaries.
#' @export
echogo_prepare_reference_inputs <- function(
    root,
    gff_file,
    protein_fasta,
    counts_file = "allcounts_table.txt",
    dge_pattern = "^dge_.*\\.csv$",
    emapper_file = NULL,
    reference_label = "reference",
    output_dir = file.path(root, "echogo_input"),
    work_dir = file.path(root, "echogo_work"),
    padj_threshold = 0.05,
    log2fc_threshold = 1,
    strip_gene_prefix = TRUE,
    counts_id_column = NULL,
    dge_id_column = NULL,
    dge_lfc_column = NULL,
    dge_padj_column = NULL,
    dge_significant_column = NULL,
    annotation_taxid = 33208L,
    blast_taxonomy = "Metazoa",
    orgdb = NULL,
    overwrite = TRUE,
    reuse_reference = TRUE,
    verbose = TRUE
) {
  resolve_input <- function(path) {
    if (length(path) != 1L || is.na(path) || !nzchar(path)) {
      stop("Input paths must be single non-empty strings.", call. = FALSE)
    }
    if (!grepl("^([A-Za-z]:)?[\\/]", path)) path <- file.path(root, path)
    normalizePath(path, winslash = "/", mustWork = FALSE)
  }
  say <- function(...) if (isTRUE(verbose)) message(...)
  normalize_gene <- function(x) {
    x <- trimws(as.character(x))
    if (isTRUE(strip_gene_prefix)) x <- sub("^gene-", "", x)
    x
  }
  first_nonempty <- function(x, fallback = NA_character_) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x) & x != "-"]
    if (length(x)) x[[1]] else fallback
  }
  clean_name <- function(x) tolower(gsub("[^a-z0-9]", "", x))
  pick_column <- function(df, explicit, candidates, label) {
    if (!is.null(explicit)) {
      if (!explicit %in% names(df)) {
        stop("Column '", explicit, "' was not found for ", label, ".", call. = FALSE)
      }
      return(explicit)
    }
    index <- match(clean_name(candidates), clean_name(names(df)))
    index <- index[!is.na(index)]
    if (!length(index)) {
      stop(
        "Could not identify the ", label, " column. Available columns: ",
        paste(names(df), collapse = ", "),
        call. = FALSE
      )
    }
    names(df)[index[[1]]]
  }
  optional_column <- function(df, candidates) {
    index <- match(clean_name(candidates), clean_name(names(df)))
    index <- index[!is.na(index)]
    if (length(index)) names(df)[index[[1]]] else NULL
  }
  require_optional <- function(packages) {
    missing <- packages[
      !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing)) {
      stop(
        "Reference-input preparation requires: ",
        paste(missing, collapse = ", "),
        ". Install with BiocManager::install(c(",
        paste(sprintf("'%s'", missing), collapse = ", "),
        ")).",
        call. = FALSE
      )
    }
  }

  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  gff_file <- resolve_input(gff_file)
  protein_fasta <- resolve_input(protein_fasta)
  counts_file <- resolve_input(counts_file)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  work_dir <- normalizePath(work_dir, winslash = "/", mustWork = FALSE)
  emapper_work_dir <- file.path(work_dir, "04_emapper")
  if (is.null(emapper_file)) {
    emapper_file <- file.path(
      emapper_work_dir,
      "BG_universe.emapper.annotations"
    )
  }
  emapper_file <- normalizePath(emapper_file, winslash = "/", mustWork = FALSE)

  required <- c(gff = gff_file, proteins = protein_fasta, counts = counts_file)
  missing_files <- names(required)[!file.exists(required)]
  if (length(missing_files)) {
    stop(
      "Missing required source files: ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.numeric(padj_threshold) || length(padj_threshold) != 1L ||
      is.na(padj_threshold) || padj_threshold <= 0 || padj_threshold > 1) {
    stop("'padj_threshold' must be in (0, 1].", call. = FALSE)
  }
  if (!is.numeric(log2fc_threshold) || length(log2fc_threshold) != 1L ||
      is.na(log2fc_threshold) || log2fc_threshold < 0) {
    stop("'log2fc_threshold' must be a non-negative number.", call. = FALSE)
  }

  require_optional(c(
    "Biostrings", "rtracklayer", "GenomicRanges", "IRanges", "S4Vectors"
  ))

  lists_dir <- file.path(work_dir, "01_lists")
  maps_dir <- file.path(work_dir, "02_maps")
  proteins_dir <- file.path(work_dir, "03_proteins")
  goseq_dir <- file.path(work_dir, "05_goseq")
  invisible(lapply(
    c(output_dir, lists_dir, maps_dir, proteins_dir, emapper_work_dir, goseq_dir),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))

  say("[1/6] Reading counts and DESeq2 contrasts...")
  counts <- .echogo_read_delim_robust(
    counts_file,
    candidates = c("\t", ",", ";")
  )
  if (!nrow(counts)) stop("The count matrix has no data rows.", call. = FALSE)
  if (is.null(counts_id_column)) counts_id_column <- names(counts)[[1]]
  if (!counts_id_column %in% names(counts)) {
    stop("Count ID column not found: ", counts_id_column, call. = FALSE)
  }
  background_genes <- unique(normalize_gene(counts[[counts_id_column]]))
  background_genes <- background_genes[
    !is.na(background_genes) & nzchar(background_genes)
  ]
  if (!length(background_genes)) {
    stop("No background gene IDs were found in the count matrix.", call. = FALSE)
  }

  dge_files <- sort(list.files(root, pattern = dge_pattern, full.names = TRUE))
  if (!length(dge_files)) {
    stop("No DESeq2 files matched dge_pattern in: ", root, call. = FALSE)
  }
  if (!file.copy(
    counts_file,
    file.path(output_dir, "allcounts_table.txt"),
    overwrite = overwrite
  )) {
    stop("Could not copy the count matrix into the EchoGO input directory.", call. = FALSE)
  }

  dge_info <- lapply(dge_files, function(path) {
    dge <- .echogo_read_delim_robust(
      path,
      candidates = c(",", "\t", ";"),
      expected = c("id", "gene", "gene_id", "log2foldchange", "padj")
    )
    id_col <- pick_column(
      dge,
      dge_id_column,
      c("ID", "gene", "gene_id", "GeneID", "symbol"),
      "DE gene ID"
    )
    ids <- normalize_gene(dge[[id_col]])
    significant_col <- dge_significant_column
    if (is.null(significant_col)) {
      significant_index <- match(
        clean_name(c("significant", "is_de", "is_significant")),
        clean_name(names(dge))
      )
      significant_index <- significant_index[!is.na(significant_index)]
      if (length(significant_index)) significant_col <- names(dge)[significant_index[[1]]]
    }
    if (!is.null(significant_col)) {
      if (!significant_col %in% names(dge)) {
        stop("Column '", significant_col, "' was not found for DE significance.", call. = FALSE)
      }
      significant <- tolower(trimws(as.character(dge[[significant_col]]))) %in%
        c("true", "t", "1", "yes", "y")
      de <- unique(ids[!is.na(ids) & nzchar(ids) & significant])
    } else {
      lfc_col <- pick_column(
        dge,
        dge_lfc_column,
        c("log2FoldChange", "log2FC", "LFC", "log2FoldChg"),
        "log2 fold-change"
      )
      padj_col <- pick_column(
        dge,
        dge_padj_column,
        c("padj", "FDR", "adj.P.Val", "adjusted_pvalue"),
        "adjusted p-value"
      )
      lfc <- suppressWarnings(as.numeric(dge[[lfc_col]]))
      padj <- suppressWarnings(as.numeric(dge[[padj_col]]))
      de <- unique(ids[
        !is.na(ids) & nzchar(ids) &
          !is.na(lfc) & abs(lfc) >= log2fc_threshold &
          !is.na(padj) & padj <= padj_threshold
      ])
    }
    overlap <- intersect(de, background_genes)
    if (!length(overlap)) {
      stop(
        "No DE genes from ", basename(path),
        " overlap the count-matrix background IDs.",
        call. = FALSE
      )
    }
    contrast <- tools::file_path_sans_ext(basename(path))
    writeLines(overlap, file.path(lists_dir, paste0(contrast, "_DE_genes.txt")))
    writeLines(background_genes, file.path(lists_dir, paste0(contrast, "_BG_genes.txt")))
    if (!file.copy(path, file.path(output_dir, basename(path)), overwrite = overwrite)) {
      stop("Could not copy DE table: ", basename(path), call. = FALSE)
    }
    list(path = path, contrast = contrast, de_genes = overlap)
  })

  gene_protein_map_file <- file.path(maps_dir, "gene_protein_map.tsv")
  gene_lengths_file <- file.path(maps_dir, "gene_lengths.tsv")
  background_fasta <- file.path(proteins_dir, "BG_universe_proteins.fa")
  reference_cache_ready <- isTRUE(reuse_reference) && all(file.exists(c(
    gene_protein_map_file,
    gene_lengths_file,
    background_fasta
  )))

  if (reference_cache_ready) {
    say("[2/6] Reusing cached reference maps and background protein FASTA...")
    gene_protein_map <- .echogo_read_delim_robust(gene_protein_map_file)
    gene_lengths <- .echogo_read_delim_robust(gene_lengths_file)
  } else {
  say("[2/6] Importing GFF3/GTF gene, exon, CDS, and protein relationships...")
  gr <- rtracklayer::import(
    gff_file,
    feature.type = c("gene", "mRNA", "transcript", "exon", "CDS")
  )
  if (!length(gr)) {
    stop("No gene, transcript, exon, or CDS records were found.", call. = FALSE)
  }
  metadata <- S4Vectors::mcols(gr)
  metadata_value <- function(key) {
    if (!key %in% colnames(metadata)) return(rep(NA_character_, length(gr)))
    value <- as.character(metadata[[key]])
    sub(",.*", "", value)
  }
  first_metadata <- function(keys) {
    output <- rep(NA_character_, length(gr))
    for (key in keys) {
      value <- metadata_value(key)
      use <- (is.na(output) | !nzchar(output)) & !is.na(value) & nzchar(value)
      output[use] <- value[use]
    }
    output
  }
  feature_type <- metadata_value("type")
  feature_id <- first_metadata(c("ID", "transcript_id", "gene_id"))
  parent_id <- metadata_value("Parent")
  direct_gene <- first_metadata(c("gene", "gene_name", "gene_id"))

  gene_rows <- feature_type == "gene"
  gene_feature_id <- feature_id[gene_rows]
  gene_name <- first_metadata(c("gene", "gene_name", "Name", "gene_id", "ID"))[
    gene_rows
  ]
  gene_name <- normalize_gene(gene_name)
  gene_lookup <- stats::setNames(gene_name, gene_feature_id)

  transcript_rows <- feature_type %in% c("mRNA", "transcript")
  transcript_id <- feature_id[transcript_rows]
  transcript_gene <- direct_gene[transcript_rows]
  missing_transcript_gene <- is.na(transcript_gene) | !nzchar(transcript_gene)
  transcript_gene[missing_transcript_gene] <- unname(gene_lookup[
    parent_id[transcript_rows][missing_transcript_gene]
  ])
  transcript_gene <- normalize_gene(transcript_gene)
  transcript_lookup <- stats::setNames(transcript_gene, transcript_id)

  resolve_feature_gene <- function(rows) {
    gene <- direct_gene[rows]
    missing <- is.na(gene) | !nzchar(gene)
    gene[missing] <- unname(transcript_lookup[parent_id[rows][missing]])
    missing <- is.na(gene) | !nzchar(gene)
    gene[missing] <- unname(gene_lookup[parent_id[rows][missing]])
    normalize_gene(gene)
  }

  cds_rows <- which(feature_type == "CDS")
  cds_gene <- resolve_feature_gene(cds_rows)
  cds_protein <- metadata_value("protein_id")[cds_rows]
  valid_cds <- !is.na(cds_gene) & nzchar(cds_gene)
  if (!any(valid_cds)) {
    stop("CDS records could not be linked to annotation gene features.", call. = FALSE)
  }
  gene_protein_map <- unique(data.frame(
    gene_id = cds_gene[valid_cds],
    protein_id = cds_protein[valid_cds],
    stringsAsFactors = FALSE
  ))
  gene_protein_map <- gene_protein_map[
    !is.na(gene_protein_map$protein_id) & nzchar(gene_protein_map$protein_id),
    ,
    drop = FALSE
  ]

  length_rows <- which(feature_type == "exon")
  length_label <- "exon"
  if (!length(length_rows)) {
    length_rows <- cds_rows
    length_label <- "CDS"
  }
  length_gene <- resolve_feature_gene(length_rows)
  valid_length <- !is.na(length_gene) & nzchar(length_gene)
  length_gr <- gr[length_rows][valid_length]
  length_gene <- length_gene[valid_length]
  if (!length(length_gr)) {
    stop("No exon or CDS records could be linked to genes.", call. = FALSE)
  }
  reduced_by_gene <- GenomicRanges::reduce(
    S4Vectors::splitAsList(length_gr, length_gene)
  )
  gene_lengths <- data.frame(
    gene_id = names(reduced_by_gene),
    length = vapply(
      reduced_by_gene,
      function(x) sum(IRanges::width(x)),
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
  say("      Bias lengths use unioned ", length_label, " intervals per gene.")
  if (!nrow(gene_protein_map)) {
    stop("No protein_id attributes could be linked to genes.", call. = FALSE)
  }
  readr::write_tsv(gene_protein_map, gene_protein_map_file)
  readr::write_tsv(gene_lengths, gene_lengths_file)

  protein_fraction <- mean(background_genes %in% gene_protein_map$gene_id)
  length_fraction <- mean(background_genes %in% gene_lengths$gene_id)
  say("      Background genes with proteins: ", round(100 * protein_fraction, 1), "%")
  say("      Background genes with lengths: ", round(100 * length_fraction, 1), "%")
  if (protein_fraction < 0.5 || length_fraction < 0.5) {
    stop(
      "Less than 50% of background genes mapped to GFF3 proteins/lengths. ",
      "Check gene ID conventions before continuing.",
      call. = FALSE
    )
  }

  say("[3/6] Writing the expressed-background protein FASTA...")
  proteins <- Biostrings::readAAStringSet(protein_fasta)
  raw_names <- names(proteins)
  bracketed <- grepl("\\[protein_id=", raw_names)
  protein_names <- sub(" .*", "", raw_names)
  protein_names[bracketed] <- sub(
    ".*\\[protein_id=([^] ]+)\\].*",
    "\\1",
    raw_names[bracketed]
  )
  names(proteins) <- protein_names
  requested_proteins <- unique(gene_protein_map$protein_id[
    gene_protein_map$gene_id %in% background_genes
  ])
  selected <- proteins[names(proteins) %in% requested_proteins]
  if (!length(selected)) {
    stop(
      "No GFF3 protein IDs matched the protein FASTA headers.",
      call. = FALSE
    )
  }
  recovery <- length(unique(names(selected))) / length(requested_proteins)
  say("      Recovered proteins: ", length(selected), " / ", length(requested_proteins))
  if (recovery < 0.5) {
    stop(
      "Less than 50% of requested proteins were found in the FASTA. ",
      "Inspect FASTA header conventions.",
      call. = FALSE
    )
  }
  Biostrings::writeXStringSet(selected, background_fasta)
  }

  command_file <- file.path(emapper_work_dir, "RUN_EMAPPER.txt")
  command <- paste(
    "emapper.py",
    "-i", shQuote(background_fasta),
    "--itype proteins",
    "-m diamond",
    "--output BG_universe",
    "--output_dir", shQuote(dirname(emapper_file)),
    "--cpu 4",
    "--override"
  )
  writeLines(command, command_file)

  base_result <- list(
    root = root,
    output_dir = output_dir,
    work_dir = work_dir,
    background_fasta = background_fasta,
    emapper_file = emapper_file,
    emapper_command = command_file,
    gene_protein_map = gene_protein_map_file,
    gene_lengths = gene_lengths_file
  )
  if (!file.exists(emapper_file)) {
    say("[4/6] eggNOG-mapper output is not present yet.")
    say("      Run the command saved at: ", command_file)
    say("      Then rerun echogo_prepare_reference_inputs() with the same arguments.")
    base_result$status <- "awaiting_emapper"
    return(invisible(base_result))
  }

  say("[4/6] Reading eggNOG-mapper annotations...")
  require_optional("goseq")
  emapper_lines <- readLines(emapper_file, warn = FALSE)
  header_lines <- grep("^#query", emapper_lines, value = TRUE)
  if (!length(header_lines)) {
    stop("Could not find the '#query' header in the eggNOG annotations.", call. = FALSE)
  }
  emapper_columns <- strsplit(
    sub("^#", "", tail(header_lines, 1L)),
    "\t",
    fixed = TRUE
  )[[1]]
  emapper <- utils::read.delim(
    emapper_file,
    sep = "\t",
    header = FALSE,
    col.names = emapper_columns,
    comment.char = "#",
    quote = "\"",
    fill = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  query_col <- pick_column(emapper, NULL, c("query", "query_name"), "eggNOG query")
  go_col <- pick_column(emapper, NULL, c("GOs", "GO", "go_terms"), "eggNOG GO")
  preferred_col <- pick_column(
    emapper,
    NULL,
    c("Preferred_name", "preferred_name"),
    "eggNOG preferred name"
  )
  seed_col <- optional_column(emapper, c("seed_ortholog", "EggNM.seed_ortholog"))
  taxonomy_col <- optional_column(emapper, c("max_annot_lvl", "EggNM.max_annot_lvl"))
  emapper_small <- data.frame(
    protein_id = as.character(emapper[[query_col]]),
    preferred_name = as.character(emapper[[preferred_col]]),
    seed_ortholog = if (is.null(seed_col)) NA_character_ else as.character(emapper[[seed_col]]),
    annotation_taxonomy = if (is.null(taxonomy_col)) {
      paste0(annotation_taxid, "|", blast_taxonomy)
    } else {
      as.character(emapper[[taxonomy_col]])
    },
    GOs = as.character(emapper[[go_col]]),
    stringsAsFactors = FALSE
  )
  joined <- merge(gene_protein_map, emapper_small, by = "protein_id")
  if (!nrow(joined)) {
    stop("No eggNOG query IDs matched the GFF3 protein IDs.", call. = FALSE)
  }

  by_gene <- split(joined, joined$gene_id)
  gene_to_go <- lapply(by_gene, function(x) {
    values <- unlist(strsplit(x$GOs[!is.na(x$GOs)], ",", fixed = TRUE))
    values <- unique(trimws(values))
    values[grepl("^GO:[0-9]{7}$", values)]
  })
  gene_to_go <- gene_to_go[vapply(gene_to_go, length, integer(1)) > 0L]
  if (!length(gene_to_go)) {
    stop("No valid GO identifiers were recovered from eggNOG annotations.", call. = FALSE)
  }
  go_annotation_source <- "eggNOG-mapper"

  if (!is.null(orgdb)) {
    if (length(orgdb) != 1L || is.na(orgdb) || !nzchar(orgdb)) {
      stop("'orgdb' must be NULL or one installed OrgDb package name.", call. = FALSE)
    }
    require_optional(orgdb)
    odb <- getExportedValue(orgdb, orgdb)
    available_keytypes <- AnnotationDbi::keytypes(odb)
    candidate_keytypes <- intersect(
      c("SYMBOL", "ALIAS", "ENTREZID", "ENSEMBL"),
      available_keytypes
    )
    gene_entrez <- stats::setNames(
      rep(NA_character_, length(background_genes)),
      background_genes
    )
    for (keytype in candidate_keytypes) {
      unresolved <- names(gene_entrez)[is.na(gene_entrez)]
      valid <- intersect(unresolved, AnnotationDbi::keys(odb, keytype = keytype))
      if (!length(valid)) next
      mapped <- AnnotationDbi::mapIds(
        odb,
        keys = valid,
        keytype = keytype,
        column = "ENTREZID",
        multiVals = "first"
      )
      gene_entrez[names(mapped)] <- unname(mapped)
    }
    gene_entrez <- gene_entrez[!is.na(gene_entrez) & nzchar(gene_entrez)]
    if (!length(gene_entrez)) {
      stop(
        "None of the count-matrix gene IDs mapped through ", orgdb,
        ". Check the organism and gene ID type.",
        call. = FALSE
      )
    }
    go_raw <- AnnotationDbi::select(
      odb,
      keys = unique(unname(gene_entrez)),
      keytype = "ENTREZID",
      columns = "GOALL"
    )
    go_raw <- go_raw[!is.na(go_raw$GOALL) & nzchar(go_raw$GOALL), , drop = FALSE]
    entrez_to_gene <- data.frame(
      gene_id = names(gene_entrez),
      ENTREZID = unname(gene_entrez),
      stringsAsFactors = FALSE
    )
    gene_go <- merge(entrez_to_gene, go_raw, by = "ENTREZID")
    gene_to_go <- lapply(
      split(gene_go$GOALL, gene_go$gene_id),
      unique
    )
    gene_to_go <- gene_to_go[vapply(gene_to_go, length, integer(1)) > 0L]
    if (!length(gene_to_go)) {
      stop("No GOALL mappings were recovered from ", orgdb, ".", call. = FALSE)
    }
    go_annotation_source <- orgdb
    say("      GOseq categories: species-native ", orgdb)
  }

  annotation <- do.call(rbind, lapply(names(by_gene), function(gene_id) {
    preferred <- first_nonempty(by_gene[[gene_id]]$preferred_name)
    seed <- first_nonempty(by_gene[[gene_id]]$seed_ortholog)
    taxonomy <- first_nonempty(
      by_gene[[gene_id]]$annotation_taxonomy,
      paste0(annotation_taxid, "|", blast_taxonomy)
    )
    data.frame(
      transcript_id = gene_id,
      gene_id = gene_id,
      native_symbol = gene_id,
      sprot_Top_BLASTX_hit = NA_character_,
      swissprot_accession = NA_character_,
      swissprot_gene_symbol = NA_character_,
      EggNM.Preferred_name = preferred,
      EggNM.seed_ortholog = seed,
      EggNM.max_annot_lvl = taxonomy,
      annotation_taxonomy = taxonomy,
      stringsAsFactors = FALSE
    )
  }))
  trinotate_file <- file.path(output_dir, "Trinotate_for_EchoGO.tsv")
  labelled_file <- file.path(
    output_dir,
    paste0(reference_label, "_eggNOG_for_EchoGO.tsv")
  )
  readr::write_tsv(annotation, trinotate_file)
  readr::write_tsv(annotation, labelled_file)

  say("[5/6] Running GOseq for ", length(dge_info), " contrast(s)...")
  go_term_info <- function(go_id, field = c("term", "ontology")) {
    field <- match.arg(field)
    object <- GO.db::GOTERM[[go_id]]
    if (is.null(object)) return(NA_character_)
    if (field == "term") AnnotationDbi::Term(object) else AnnotationDbi::Ontology(object)
  }
  summaries <- lapply(dge_info, function(info) {
    analysis_genes <- Reduce(
      intersect,
      list(background_genes, gene_lengths$gene_id, names(gene_to_go))
    )
    if (length(analysis_genes) < 10L) {
      stop(
        "Fewer than 10 background genes have both lengths and GO annotations.",
        call. = FALSE
      )
    }
    de_genes <- intersect(info$de_genes, analysis_genes)
    if (!length(de_genes)) {
      stop("No DE genes remain after GO/length filtering for ", info$contrast, ".", call. = FALSE)
    }
    lengths <- stats::setNames(gene_lengths$length, gene_lengths$gene_id)
    lengths <- lengths[analysis_genes]
    indicator <- as.integer(analysis_genes %in% de_genes)
    names(indicator) <- analysis_genes
    pwf <- suppressWarnings(goseq::nullp(
      DEgenes = indicator,
      bias.data = lengths,
      plot.fit = FALSE
    ))
    enrichment <- suppressWarnings(goseq::goseq(
      pwf,
      gene2cat = gene_to_go[analysis_genes],
      use_genes_without_cat = FALSE
    ))
    enrichment$over_represented_FDR <- stats::p.adjust(
      enrichment$over_represented_pvalue,
      method = "BH"
    )
    keep <- !is.na(enrichment$over_represented_FDR) &
      enrichment$over_represented_FDR <= 0.05
    enrichment <- enrichment[keep, , drop = FALSE]
    output <- data.frame(
      category = as.character(enrichment$category),
      term = vapply(enrichment$category, go_term_info, character(1), field = "term"),
      ontology = vapply(enrichment$category, go_term_info, character(1), field = "ontology"),
      numDEInCat = enrichment$numDEInCat,
      numInCat = enrichment$numInCat,
      total_significant_genes = rep(length(de_genes), nrow(enrichment)),
      total_tested_genes = rep(length(analysis_genes), nrow(enrichment)),
      over_represented_FDR = enrichment$over_represented_FDR,
      gene_ids = vapply(enrichment$category, function(go_id) {
        genes <- names(gene_to_go)[
          vapply(gene_to_go, function(x) go_id %in% x, logical(1))
        ]
        paste(intersect(genes, de_genes), collapse = ", ")
      }, character(1)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    goseq_file <- file.path(
      output_dir,
      paste0(info$contrast, ".GOseq.enriched.tsv")
    )
    readr::write_tsv(output, goseq_file)
    if (!nrow(output)) {
      warning("GOseq returned no enriched rows for ", info$contrast, ".")
    }
    data.frame(
      contrast = info$contrast,
      background_genes = length(analysis_genes),
      de_genes = length(de_genes),
      enriched_terms = nrow(output),
      go_annotation_source = go_annotation_source,
      goseq_file = goseq_file,
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summaries)
  readr::write_csv(summary, file.path(goseq_dir, "summary_goseq_per_contrast.csv"))

  say("[6/6] Validating the final EchoGO input bundle...")
  expected_columns <- c(
    "category", "term", "ontology", "numDEInCat", "numInCat",
    "total_significant_genes", "total_tested_genes",
    "over_represented_FDR", "gene_ids"
  )
  for (path in summary$goseq_file) {
    table <- utils::read.delim(
      path,
      sep = "\t",
      quote = "\"",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    if (!identical(names(table), expected_columns)) {
      stop("Generated GOseq table has an unexpected schema: ", basename(path), call. = FALSE)
    }
  }
  final_files <- list.files(output_dir)
  if (!all(c("allcounts_table.txt", "Trinotate_for_EchoGO.tsv") %in% final_files)) {
    stop("The final EchoGO input bundle is incomplete.", call. = FALSE)
  }
  say("      EchoGO input bundle completed at: ", output_dir)

  base_result$status <- "complete"
  base_result$trinotate_file <- trinotate_file
  base_result$labelled_annotation_file <- labelled_file
  base_result$summary <- summary
  base_result$go_annotation_source <- go_annotation_source
  invisible(base_result)
}

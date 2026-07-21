# Build the deterministic EchoGO v0.1.3 matched-universe demonstration fixture.

suppressPackageStartupMessages({
  library(readr)
  library(jsonlite)
})

package_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
demo_dir <- file.path(package_root, "inst", "extdata", "echogo_demo")
archive_dir <- file.path(
  package_root, "audit", "v0.1.3_route1_implementation",
  "legacy_demo_inputs_pre_v013"
)
dir.create(demo_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)

legacy_files <- list.files(demo_dir, full.names = TRUE, recursive = FALSE)
legacy_files <- legacy_files[file.info(legacy_files)$isdir %in% FALSE]
for (path in legacy_files) {
  destination <- file.path(archive_dir, basename(path))
  if (!file.exists(destination)) file.copy(path, destination, overwrite = FALSE)
}

ids <- sprintf("tx%03d", seq_len(120L))
symbols <- c(
  "ISY1", "BUD13", "SF3B6", "U2AF2", "SF3B1", "DDX46", "SF3A2", "DDX42",
  "PUF60", "SF3B3", "SNIP1", "SLU7", "SRPK3", "SF3A1", "SF3A3", "PHF5A",
  "DDX23", "RBMX2", "NCBP1", "DDX1", "SF3B4", "SF3B2", "SF3B5", "SRSF7",
  "PTBP1", "LSM7", "RBM19", "TIA1", "RBM25", "RBM42", "DCPS", "LSM8",
  "CWC22", "DHX8", "LSM6", "SDE2", "RBM47", "RBM22", "RBM10", "PRP4K",
  "ZMAT2", "AQR", "SMU1", "DBR1", "RBM39", "LARP7", "PPIE", "RBM11",
  "SNW1", "LSM5", "NCBP2", "TRA2B", "NUP98", "CDC5L", "LSM4", "WBP4",
  "SRRM2", "LSM3", "PRPF3", "MFAP1", "SYF2", "NOVA1", "LSM2", "SART3",
  "SNRPA", "PQBP1", "RBM8A", "SRRM1", "NHP2", "RBM3", "SENP7", "U2AF1",
  "ZCRB1", "DDX17", "SRSF3", "PLRG1", "PPIL1", "DDX41", "UBL5", "BUD31",
  "YJU2", "DHX15", "CWC27", "MTREX", "TGS1", "ESS2", "MBNL1", "GLOD4",
  "PSMD1", "ADRM1", "PSMD2", "PSMD4", "PSMD7", "PSMD3", "PSMD6", "PSMD9",
  "PSMD8", "PSME3", "UBE3C", "UCHL5", "SEM1", "NDE1", "TTK", "CKAP5",
  "PLK1", "SIN3A", "KAT8", "CBX5", "PINX1", "KAT5", "DCTN1", "ORC2",
  "CDC20", "LRWD1", "CENPE", "BUB3", "XPO1", "CHMP3"
)
stopifnot(length(symbols) == 118L)
symbols[[20]] <- symbols[[19]]

counts <- data.frame(
  gene_id = ids,
  control_1 = 40L + (seq_along(ids) * 7L) %% 91L,
  control_2 = 45L + (seq_along(ids) * 11L) %% 89L,
  treated_1 = 50L + (seq_along(ids) * 13L) %% 97L,
  treated_2 = 55L + (seq_along(ids) * 17L) %% 101L,
  check.names = FALSE
)

significant <- seq_along(ids) <= 24L
de <- data.frame(
  gene_id = ids,
  baseMean = rowMeans(counts[-1]),
  log2FoldChange = ifelse(significant, ifelse(seq_along(ids) %% 2L, 2.2, -2.2), 0.25),
  lfcSE = 0.25,
  stat = ifelse(significant, 8.8, 1),
  pvalue = ifelse(significant, 1e-6, 0.3),
  padj = ifelse(significant, 0.001, 0.4),
  significant = significant,
  stringsAsFactors = FALSE
)
de$padj[[30]] <- NA_real_

annotation <- data.frame(
  gene_id = ids[seq_len(118L)],
  transcript_id = ids[seq_len(118L)],
  sprot_Top_BLASTX_hit = paste0(symbols, "^Metazoa"),
  EggNM.Preferred_name = symbols,
  EggNM.max_annot_lvl = "33208|Metazoa",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

term_ids <- c("GO:0003674", "GO:0005488", "GO:0008152", "GO:0009987", "GO:0006355")
term_names <- c(
  "molecular function", "binding", "metabolic process", "cellular process",
  "regulation of DNA-templated transcription"
)
ontologies <- c("MF", "MF", "BP", "BP", "BP")
num_de <- c(24L, 12L, 6L, 2L, 1L)
num_in <- c(24L, 30L, 60L, 20L, 40L)
goseq <- data.frame(
  category = term_ids,
  term = term_names,
  ontology = ontologies,
  numDEInCat = num_de,
  numInCat = num_in,
  total_significant_genes = 24L,
  total_tested_genes = 120L,
  over_represented_FDR = c(0.001, 0.004, 0.01, 0.02, 0.04),
  gene_ids = vapply(num_de, function(n) paste(ids[seq_len(n)], collapse = ", "), character(1)),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

write_tsv(counts, file.path(demo_dir, "counts_demo.tsv"))
write_tsv(de, file.path(demo_dir, "DE_results_demo.tsv"))
write_tsv(annotation, file.path(demo_dir, "Trinotate_demo.tsv"))
write_tsv(goseq, file.path(demo_dir, "GOseq_enrichment_demo.tsv"))

foreground <- unique(symbols[seq_len(24L)])
background <- unique(symbols)
cache_root <- file.path(demo_dir, "cached_gprofiler_v0.1.3")
bg_dir <- file.path(cache_root, "with_custom_background")
nobg_dir <- file.path(cache_root, "no_background_genome_wide")
dir.create(bg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(nobg_dir, recursive = TRUE, showWarnings = FALSE)

species <- c("hsapiens", "mmusculus", "drerio")
manifest_runs <- list()
relative_cache <- function(path) {
  root <- paste0(normalizePath(cache_root, winslash = "/", mustWork = FALSE), "/")
  sub(paste0("^", root), "", normalizePath(path, winslash = "/", mustWork = FALSE))
}
write_run <- function(organism, mode, index) {
  custom <- identical(mode, "custom_background")
  suffix <- if (custom) "with_bg" else "nobg"
  directory <- if (custom) bg_dir else nobg_dir
  stem <- file.path(directory, paste0("gprofiler_", organism, "_", suffix))
  query_file <- paste0(stem, "_query.txt")
  background_file <- if (custom) paste0(stem, "_background.txt") else NULL
  result_file <- paste0(stem, ".csv")
  metadata_file <- paste0(stem, "_metadata.json")
  writeLines(foreground, query_file)
  if (custom) writeLines(background, background_file)
  query_size <- length(foreground) - index + 1L
  domain_size <- if (custom) length(background) - index + 1L else 18000L + 1000L * index
  results <- data.frame(
    term_id = c(term_ids, "GO:0006397"),
    term_name = c(term_names, "mRNA processing"),
    source = "GO:BP",
    p_value = c(0.001, 0.004, 0.01, 0.02, 0.04, 0.03) * (1 + index / 10),
    intersection = vapply(c(10L, 8L, 6L, 4L, 2L, 5L), function(n) paste(foreground[seq_len(n)], collapse = ","), character(1)),
    intersection_size = c(10L, 8L, 6L, 4L, 2L, 5L),
    query_size = query_size,
    term_size = c(18L, 22L, 30L, 25L, 40L, 16L),
    effective_domain_size = domain_size,
    depth = c(1L, 2L, 3L, 4L, 7L, 8L),
    species_code = organism,
    species_label = organism,
    stringsAsFactors = FALSE
  )
  results$fold_enrichment <- with(
    results,
    (intersection_size / query_size) / (term_size / effective_domain_size)
  )
  write_csv(results, result_file)
  query_hash <- unname(tools::md5sum(query_file))
  background_hash <- if (custom) unname(tools::md5sum(background_file)) else NULL
  metadata <- list(
    echogo_version = "0.1.3",
    execution = "cached_demo_fixture",
    cache_version = "v0.1.3",
    organism_code = organism,
    organism_label = organism,
    requested_sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG"),
    correction_method = "fdr",
    user_threshold = 0.05,
    significant = FALSE,
    ordered_query = FALSE,
    multi_query = FALSE,
    background_mode = mode,
    submitted_foreground_count = length(foreground),
    submitted_background_count = if (custom) length(background) else NULL,
    effective_query_size = query_size,
    effective_domain_size = domain_size,
    de_significance_rule = list(type = "logical_column", column = "significant"),
    resolver_definition = list(
      route = "shared_portable_canonical_organism_context",
      explicit_species_specific_ortholog_mapping = FALSE
    ),
    query_hash = query_hash,
    background_hash = background_hash,
    response_metadata = list(
      source = "deterministic offline demonstration fixture",
      database_version = "not_a_live_gprofiler_release"
    ),
    status = "cached_fixture"
  )
  write_json(metadata, metadata_file, auto_unbox = TRUE, pretty = TRUE, na = "null")
  list(
    organism_code = organism,
    organism_label = organism,
    background_mode = mode,
    submitted_foreground_count = length(foreground),
    submitted_background_count = if (custom) length(background) else NULL,
    effective_query_size = query_size,
    effective_domain_size = domain_size,
    query_hash = query_hash,
    background_hash = background_hash,
    status = "cached_fixture",
    result_file = relative_cache(result_file),
    metadata_file = relative_cache(metadata_file),
    query_file = relative_cache(query_file),
    background_file = if (custom) relative_cache(background_file) else NULL
  )
}

for (i in seq_along(species)) {
  manifest_runs[[length(manifest_runs) + 1L]] <- write_run(species[[i]], "custom_background", i)
  manifest_runs[[length(manifest_runs) + 1L]] <- write_run(species[[i]], "no_background_genome_wide", i)
}
write_csv(
  expand.grid(species_code = species, mode = c("with_bg", "nobg"), stringsAsFactors = FALSE),
  file.path(cache_root, "gprofiler_enrichment_summary.csv")
)
write_json(
  list(
    schema_version = "1.0",
    echogo_version = "0.1.3",
    execution = "cached_demo_fixture",
    generated = "deterministic",
    vector_contract = "shared_portable_canonical_organism_context",
    explicit_species_specific_ortholog_mapping = FALSE,
    submitted_foreground_count = length(foreground),
    submitted_background_count = length(background),
    significance_rule = list(type = "logical_column", column = "significant"),
    resolver_definition = list(
      route = "shared_portable_canonical_organism_context",
      explicit_species_specific_ortholog_mapping = FALSE
    ),
    runs = manifest_runs
  ),
  file.path(cache_root, "run_manifest.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  na = "null"
)

readme <- c(
  "EchoGO v0.1.3 deterministic demonstration dataset",
  "=================================================",
  "",
  "This synthetic scientific fixture contains no research measurements or",
  "machine-specific paths.",
  "",
  "Contract checks represented by the fixture:",
  "- 120 unique tested genes in counts and the full DE table",
  "- 24 significant genes using the explicit logical 'significant' column",
  "- 118 annotation-matched genes, two intentionally unmapped tested genes",
  "- one intentional many-to-one canonical-name collision",
  "- GOseq rows with explicit 24/120 denominators; the first term has fold 5",
  "- a strict resolved foreground subset of the resolved background",
  "- deterministic cached g:Profiler responses for offline quickstart testing",
  "",
  "Offline quickstart (default):",
  "  echogo_quickstart(run_demo = TRUE)",
  "",
  "Optional live integration run:",
  "  echogo_quickstart(run_demo = TRUE, live_gprofiler = TRUE)",
  "",
  "Regeneration:",
  "  source('data-raw/build_demo_v013.R')"
)
writeLines(readme, file.path(demo_dir, "README_demo.txt"))

stopifnot(
  nrow(read_tsv(file.path(demo_dir, "counts_demo.tsv"), show_col_types = FALSE)) == 120L,
  nrow(read_tsv(file.path(demo_dir, "DE_results_demo.tsv"), show_col_types = FALSE)) == 120L,
  nrow(read_tsv(file.path(demo_dir, "Trinotate_demo.tsv"), show_col_types = FALSE)) == 118L,
  nrow(read_tsv(file.path(demo_dir, "GOseq_enrichment_demo.tsv"), show_col_types = FALSE)) == 5L
)
message("EchoGO v0.1.3 demo written to: ", normalizePath(demo_dir, winslash = "/"))

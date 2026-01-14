---
title: "Reference-based RNA-seq inputs for EchoGO"
output:
  rmarkdown::html_vignette:
    toc: true
    toc_depth: 2
    css: echogo-vignette.css
vignette: >
  %\VignetteIndexEntry{Reference-based RNA-seq inputs for EchoGO}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
editor_options:
  markdown:
    wrap: 72
---

\`\`\`{r setup, include=FALSE} knitr::opts_chunk\$set( collapse = TRUE,
comment = "#\>", error = TRUE, message = FALSE, warning = FALSE )

# Heavy steps are disabled during vignette build / CRAN.

# Run locally with: Sys.setenv(NOT_CRAN="true")

RUN_PIPELINE \<- identical(Sys.getenv("NOT_CRAN"), "true")

\`\`\`md

::: echogo-hero
<p class="echogo-badges">

<a href="https://doi.org/10.5281/zenodo.17660708" target="_blank" rel="noopener">
<img src="zenodo-badge.svg" alt="DOI badge"/> </a> [DOI:
<a href="https://doi.org/10.5281/zenodo.17660708" target="_blank" rel="noopener">
10.5281/zenodo.17660708 </a>]{.echogo-doi}

</p>

<h1 class="echogo-title">

🧬 EchoGO

</h1>

<p class="echogo-subtitle">

Reference-based RNA-seq inputs (HISAT2/STAR + featureCounts + DESeq2)

</p>
:::

::: {.alert .alert-info}
<strong>What this vignette does:</strong> prepares the exact
<code>input/</code> bundle EchoGO expects from a <strong>reference-based
RNA-seq</strong> workflow. EchoGO does not run aligners; it consumes
your outputs and standardizes the input folder.
:::

::: {.alert .alert-warning}
<strong>CRAN note:</strong> heavy steps are disabled during vignette
build. To run everything locally set
<code>Sys.setenv(NOT_CRAN="true")</code>.
:::

## Overview {.tabset}

### What you need before you run this

In one project folder (`root`), you must already have:

-   `allcounts_table.txt`\
    Gene-level count matrix (e.g., featureCounts output). First column
    must be the gene ID.
-   `dge_*.csv`\
    One DESeq2 result table per contrast (e.g.,
    `dge_Treated_vs_Control.csv`).
-   `GENOME_ANNOTATION.gff.gz`\
    Genome annotation in GFF3 format; ideally CDS features contain
    `protein_id=` attributes (common in NCBI GFF3).
-   `REFERENCE_PROTEINS.faa.gz`\
    Protein FASTA matching the GFF3 (e.g., NCBI
    `*_translated_cds.faa.gz` or equivalent).

### One manual external step (eggNOG-mapper)

This script will generate:

-   `echogo_work/03_proteins/BG_universe_proteins.fa`

Then you run eggNOG-mapper externally (Linux/WSL), producing:

-   `echogo_work/04_emapper/BG_universe.emapper.annotations`

After that, re-run the script (or continue from STEP 8 onward).

### What you should have at the end (EchoGO input bundle)

Inside `<root>/<echogo_input_folder>/` you should see:

-   `allcounts_table.txt`
-   `dge_<CONTRAST>.csv`
-   `dge_<CONTRAST>.GOseq.enriched.tsv`
-   `dge_<CONTRAST>.GOseq.depleted.tsv` (only if depleted terms exist)
-   `Trinotate_for_EchoGO.tsv`
-   `<reference_label>_eggNOG_for_EchoGO.tsv`

## Quickstart {.tabset}

### Minimal run (recommended)

1)  Edit the **USER SETTINGS** block\
2)  Run the script once → it will stop at STEP 7 with instructions\
3)  Run eggNOG-mapper externally\
4)  Re-run the script → it will finish and write the final
    `echogo_input/` folder.

### Run EchoGO after bundle exists

\`\`\`r library(EchoGO)

input_dir \<- "E:/your_reference_project/echogo_input" outdir \<-
"E:/your_reference_project/echogo_results"

res \<- run_full_echogo( input_dir = input_dir, outdir = outdir,
make_report = TRUE, verbose = TRUE )

suppressPackageStartupMessages({ library(readr) library(dplyr)
library(stringr) library(purrr) library(tibble) library(Biostrings)
library(goseq) library(GO.db) library(AnnotationDbi) library(qvalue) })

# ==============================================================================

# USER SETTINGS — EDIT THESE PATHS / VALUES

# ==============================================================================

# Project root: folder containing your reference-based RNA-seq outputs

# Must include: allcounts_table.txt and one or more dge\_\*.csv files

root \<- "E:/your_reference_project"

# A short label for your reference species (used to name the eggNOG-derived file)

# Examples: "Trout", "Zebrafish", "MySpecies"

reference_label \<- "MySpecies"

# Genome annotation (GFF3 with CDS features containing protein_id attributes)

gff_path \<- file.path(root, "GENOME_ANNOTATION.gff.gz")

# Reference protein FASTA (translated CDS / proteins)

# NCBI: \*\_translated_cds.faa.gz, Ensembl: \*.pep.all.fa.gz, etc.

ref_protein_fa \<- file.path(root, "REFERENCE_PROTEINS.faa.gz")

# Where eggNOG-mapper output will live (after you run emapper.py manually)

emapper_path \<- file.path( root, "echogo_work", "04_emapper",
"BG_universe.emapper.annotations" )

# Final EchoGO input folder name (created inside `root`)

echogo_input_folder \<- "echogo_input"

# Thresholds for defining DE genes from DESeq2 result tables

padj_thr \<- 0.05 l2fc_thr \<- 1

# If your gene IDs are prefixed like "gene-LOC123", set TRUE

strip_gene_prefix \<- TRUE

# Number of CPUs you plan to use for emapper (printed in instructions only)

emapper_cpus \<- 4

# ==============================================================================

# END USER SETTINGS

# ==============================================================================

# ==============================================================================

# Directory layout

# - intermediates in <root>/echogo_work/

# - final EchoGO input bundle in <root>/<echogo_input_folder>/

# ==============================================================================

work_dir \<- file.path(root, "echogo_work") lists_dir \<-
file.path(work_dir, "01_lists") \# DE/BG lists per contrast maps_dir \<-
file.path(work_dir, "02_maps") \# gene-protein map + gene lengths
proteins_dir \<- file.path(work_dir, "03_proteins") \# BG protein FASTA
emapper_dir \<- file.path(work_dir, "04_emapper") \# eggNOG output
goseq_dir \<- file.path(work_dir, "05_goseq") \# GOseq summaries
(optional) echogo_input \<- file.path(root, echogo_input_folder) \#
FINAL EchoGO input bundle

dirs \<- c(work_dir, lists_dir, maps_dir, proteins_dir, emapper_dir,
goseq_dir, echogo_input) invisible(lapply(dirs, dir.create, recursive =
TRUE, showWarnings = FALSE))

# ==============================================================================

# Helpers

# ==============================================================================

strip_gene \<- function(x) { if (!isTRUE(strip_gene_prefix)) return(x)
sub("\^gene-", "", x) }

nzc \<- function(x) x[!is.na(x) & nzchar(x)]

message("==============================================================================")
message("EchoGO — Reference-based RNA-seq input preparation")
message("==============================================================================\n")

# ==============================================================================

# STEP 1: Verify required files exist

# ==============================================================================

message("STEP 1: Verifying input files...")

counts_path \<- file.path(root, "allcounts_table.txt")
stopifnot(file.exists(counts_path)) stopifnot(file.exists(gff_path))
stopifnot(file.exists(ref_protein_fa))

dge_files \<- list.files(root, pattern = "\^dge\_.\*\\.csv\$",
full.names = TRUE) stopifnot(length(dge_files) \> 0)

message(" ✓ Counts file: ", basename(counts_path)) message(" ✓ GFF3
file: ", basename(gff_path)) message(" ✓ Protein FASTA: ",
basename(ref_protein_fa)) message(" ✓ DESeq2 DGE files: ",
length(dge_files), "\n")

# ==============================================================================

# STEP 2: Copy raw inputs into the final EchoGO input folder

# ==============================================================================

message("STEP 2: Copying raw inputs into EchoGO input folder...")

file.copy(counts_path, file.path(echogo_input, "allcounts_table.txt"),
overwrite = TRUE)

copied \<- file.copy( dge_files, file.path(echogo_input,
basename(dge_files)), overwrite = TRUE ) stopifnot(all(copied))

message(" ✓ Copied allcounts_table.txt + ", length(dge_files), "
dge\_\*.csv to:\n ", echogo_input, "\n")

# ==============================================================================

# STEP 3: Extract DE and background gene lists (written to echogo_work/01_lists)

# ==============================================================================

message("STEP 3: Extracting DE and background gene lists...")

count_mat \<- readr::read_tsv(counts_path, col_types = readr::cols())
gene_id_col \<- names(count_mat)[1]

bg_universe \<- unique(count_mat[[gene_id_col]]) bg_universe \<-
strip_gene(bg_universe) bg_universe \<- nzc(bg_universe)

message(" ✓ Background universe: ", length(bg_universe), " genes")

read_dge_table \<- function(path) { x \<- try(readr::read_csv(path,
guess_max = 100000, show_col_types = FALSE), silent = TRUE) if
(inherits(x, "try-error")) { message(" read_csv failed; trying
read_csv2: ", basename(path)) x \<- readr::read_csv2(path, guess_max =
100000, show_col_types = FALSE) return(x) } if (ncol(x) == 1L &&
grepl(";", names(x)[1])) { message(" Detected ';' separation; re-reading
with read_csv2: ", basename(path)) x \<- readr::read_csv2(path,
guess_max = 100000, show_col_types = FALSE) } x }

standardize_dge_columns \<- function(df, path) { orig_names \<-
names(df)

clean_names \<- orig_names \|\> stringr::str_squish() \|\>
stringr::str_replace_all("[[^1]\_]+", "") \|\> toupper()

[^1]: :alnum:

names(df) \<- clean_names

id_candidates \<- c("ID", "GENE", "GENEID", "SYMBOL") lfc_candidates \<-
c("LFC", "LOG2FC", "LOG2FOLDCHANGE", "LOG2FOLDCHG") padj_candidates \<-
c("PADJ", "QVALUE", "FDR", "ADJ_PVAL", "PVAL_ADJ")

pick_first \<- function(cands, available) { hit \<- intersect(cands,
available) if (!length(hit)) return(NA_character\_) hit[[1]] }

id_col \<- pick_first(id_candidates, names(df)) lfc_col \<-
pick_first(lfc_candidates, names(df)) padj_col \<-
pick_first(padj_candidates, names(df))

if (any(is.na(c(id_col, lfc_col, padj_col)))) { stop( "Could not find
ID/LFC/PADJ columns in: ", path, "\nOriginal names: ", paste(orig_names,
collapse = ", "), "\nCleaned names: ", paste(names(df), collapse = ",
"), call. = FALSE ) }

df$ID   <- df[[id_col]]
  df$LFC \<- df[[lfc_col]] df\$PADJ \<- df[[padj_col]] df }

process_dge_file \<- function(f) { df \<- read_dge_table(f) df \<-
standardize_dge_columns(df, f)

de \<- df \|\> dplyr::filter( !is.na(PADJ), PADJ \<= padj_thr,
!is.na(LFC), abs(LFC) \>= l2fc_thr ) \|\> dplyr::pull(ID) \|\> unique()
\|\> strip_gene() \|\> nzc()

base \<- tools::file_path_sans_ext(basename(f)) de_file \<-
file.path(lists_dir, paste0(base, "\_DE_genes.txt")) bg_file \<-
file.path(lists_dir, paste0(base, "\_BG_genes.txt"))

writeLines(de, de_file) writeLines(bg_universe, bg_file)

tibble::tibble( contrast = base, n_bg = length(bg_universe), n_de =
length(de) ) }

dge_summary \<- purrr::map_dfr(dge_files, process_dge_file)
readr::write_csv(dge_summary, file.path(lists_dir,
"summary_counts_per_contrast.csv"))

print(dge_summary) cat("\n")

# ==============================================================================

# STEP 4: Build gene → protein map from GFF3 (written to echogo_work/02_maps)

# ==============================================================================

message("STEP 4: Building gene → protein map from GFF3...")

out_map \<- file.path(maps_dir, "gene_protein_map.tsv")

con \<- gzfile(gff_path, open = "rt") gff_lines \<- readLines(con)
close(con)

gff_body \<- gff_lines[!startsWith(gff_lines, "\#")] tab_split \<-
strsplit(gff_body, "\t", fixed = TRUE) tab_split \<-
tab_split[vapply(tab_split, length, integer(1)) \>= 9]

gff \<- tibble::tibble( type = vapply(tab_split, `[`, character(1), 3),
start = vapply(tab_split, `[`, character(1), 4), end = vapply(tab_split,
`[`, character(1), 5), attrs = vapply(tab_split, `[`, character(1), 9) )

parse_attrs \<- function(attr_str) { fields \<- strsplit(attr_str, ";",
fixed = TRUE)[[1]] fields \<- trimws(fields) fields \<-
fields[nzchar(fields) & grepl("=", fields, fixed = TRUE)] kv \<-
strsplit(fields, "=", fixed = TRUE) keys \<- vapply(kv, `[`,
character(1), 1) vals \<- vapply(kv, `[`, character(1), 2)
stats::setNames(as.list(vals), keys) }

attr_list \<- lapply(gff\$attrs, parse_attrs) get_or_na \<- function(a,
key) if (!is.null(a[[key]])) a[[key]] else NA_character\_

# Genes

is_gene \<- gff\$type == "gene" gene_df \<- tibble::tibble(attrs =
attr_list[is_gene]) \|\> dplyr::mutate( gene_id_full = vapply(attrs,
get_or_na, character(1), "ID"), gene_name = vapply(attrs, get_or_na,
character(1), "Name") ) \|\> dplyr::mutate( gene_id =
strip_gene(gene_id_full), gene_id = ifelse(!is.na(gene_name) &
nzchar(gene_name), gene_name, gene_id), gene_id = strip_gene(gene_id) )
\|\> dplyr::select(gene_id, gene_id_full) \|\>
dplyr::filter(!is.na(gene_id), nzchar(gene_id)) \|\> dplyr::distinct()

# mRNA/transcripts

is_mrna \<- gff\$type %in% c("mRNA", "transcript") mrna_df \<-
tibble::tibble(attrs = attr_list[is_mrna]) \|\> dplyr::mutate( mrna_id =
vapply(attrs, get_or_na, character(1), "ID"), gene_id_full =
vapply(attrs, get_or_na, character(1), "Parent") ) \|\>
dplyr::filter(!is.na(mrna_id), !is.na(gene_id_full)) \|\>
dplyr::select(mrna_id, gene_id_full) \|\> dplyr::distinct()

# CDS → protein_id

is_cds \<- gff\$type == "CDS" cds_df \<- tibble::tibble(attrs =
attr_list[is_cds]) \|\> dplyr::mutate( mrna_id = vapply(attrs,
get_or_na, character(1), "Parent"), protein_id = vapply(attrs,
get_or_na, character(1), "protein_id") ) \|\>
dplyr::filter(!is.na(mrna_id), !is.na(protein_id), nzchar(protein_id))
\|\> dplyr::select(mrna_id, protein_id) \|\> dplyr::distinct()

# Link gene → protein via mRNA

gene_to_prot \<- mrna_df \|\> dplyr::inner_join(cds_df, by = "mrna_id")
\|\> dplyr::inner_join(gene_df, by = "gene_id_full") \|\>
dplyr::select(gene_id, protein_id) \|\> dplyr::distinct()

readr::write_tsv(gene_to_prot, out_map) message(" ✓ Wrote
gene_protein_map.tsv with ", nrow(gene_to_prot), " mappings\n")

# ==============================================================================

# STEP 5: Gene lengths from GFF3 (written to echogo_work/02_maps)

# ==============================================================================

message("STEP 5: Calculating gene lengths for GOseq...")

gene_len_file \<- file.path(maps_dir, "gene_lengths.tsv")

mrna_gene \<- mrna_df \|\> dplyr::mutate(gene_id_full =
strip_gene(gene_id_full)) \|\> dplyr::inner_join(gene_df, by =
"gene_id_full") \|\> dplyr::select(mrna_id, gene_id)

cds_coord \<- tibble::tibble( start =
suppressWarnings(as.integer(gff$start[is_cds])),
  end   = suppressWarnings(as.integer(gff$end[is_cds])), attrs =
attr_list[is_cds] ) \|\> dplyr::mutate( mrna_id = vapply(attrs,
get_or_na, character(1), "Parent"), cds_len = (end - start + 1L) ) \|\>
dplyr::filter(!is.na(mrna_id), !is.na(cds_len), cds_len \> 0) \|\>
dplyr::select(mrna_id, cds_len)

gene_lengths \<- cds_coord \|\> dplyr::inner_join(mrna_gene, by =
"mrna_id") \|\> dplyr::group_by(gene_id) \|\> dplyr::summarise(length =
sum(cds_len), .groups = "drop") \|\> dplyr::arrange(gene_id)

bg_with_len \<- intersect(bg_universe, gene_lengths\$gene_id) dropped_bg
\<- setdiff(bg_universe, bg_with_len)

if (length(dropped_bg) \> 0L) { message(" ⚠ Dropping ",
length(dropped_bg), " BG genes with no CDS length info") }

readr::write_tsv(gene_lengths, gene_len_file) message(" ✓ Wrote
gene_lengths.tsv with ", nrow(gene_lengths), " genes\n")

# ==============================================================================

# STEP 6: Build universal BG protein FASTA (written to echogo_work/03_proteins)

# ==============================================================================

message("STEP 6: Creating background protein FASTA for
eggNOG-mapper...")

bg_fasta \<- file.path(proteins_dir, "BG_universe_proteins.fa")

gene_prot_map \<- readr::read_tsv(out_map, show_col_types = FALSE) \|\>
dplyr::distinct(gene_id, protein_id)

bg_prot \<- tibble::tibble(gene_id = bg_universe) \|\>
dplyr::inner_join(gene_prot_map, by = "gene_id") \|\>
dplyr::distinct(protein_id) \|\> dplyr::arrange(protein_id)

all_prot \<- Biostrings::readAAStringSet(ref_protein_fa)

# Try to extract protein_id from NCBI-style headers: [protein_id=XP\_...]

prot_ids \<- sub(".*\\[protein_id=([\^\\\\]]+)\\].*", "\\1",
names(all_prot)) valid \<- nzchar(prot_ids) & !is.na(prot_ids)

if (mean(valid) \< 0.2) { message(" ⚠ Protein ID extraction matched
\<20% headers; using raw FASTA names as IDs") prot_ids \<-
names(all_prot) valid \<- nzchar(prot_ids) & !is.na(prot_ids) }

all_prot \<- all_prot[valid] prot_ids \<- prot_ids[valid]
names(all_prot) \<- prot_ids

keep \<- all_prot[names(all_prot) %in% bg_prot\$protein_id]
Biostrings::writeXStringSet(keep, bg_fasta, format = "fasta")

message(" ✓ BG proteins requested: ", length(bg_prot\$protein_id))
message(" ✓ BG proteins found: ", length(keep)) message(" ✓ Wrote: ",
bg_fasta, "\n")

# ==============================================================================

# STEP 7: Run eggNOG-mapper (manual external step)

# ==============================================================================

message("==============================================================================")
message("STEP 7: MANUAL STEP — Run eggNOG-mapper on
BG_universe_proteins.fa")
message("==============================================================================\n")

message("Run this in WSL/Linux (adjust conda env, CPUs, paths as
needed):\n") message(" conda activate eggnog") message(" cd ",
proteins_dir) message(" emapper.py \\") message(" -i
BG_universe_proteins.fa \\") message(" --itype proteins \\") message("
-o BG_universe \\") message(" --output_dir ", emapper_dir, " \\")
message(" --cpu ", emapper_cpus, "\n")

message("Expected output: ", emapper_path, "\n") message("When that file
exists, re-run this script (or continue from STEP 8).\n")

if (!file.exists(emapper_path)) { message("STOP: eggNOG output not found
yet at:\n ", emapper_path) message("Run emapper.py first, then re-run.")
quit(save = "no") }

# ==============================================================================

# STEP 8: Build Trinotate-like annotation table (final files in echogo_input/)

# ==============================================================================

message("STEP 8: Building Trinotate-like annotation table from
eggNOG...")

read_emapper_annotations \<- function(path) { emap_try \<-
readr::read_tsv(path, comment = "\#", show_col_types = FALSE) qcols \<-
intersect(c("#query", "query", "query_name"), names(emap_try)) if
(length(qcols)) { attr(emap_try, "query_col") \<- qcols[1]
return(emap_try) }

message(" No query column detected — re-reading as headerless emapper
output...") emap \<- readr::read_tsv(path, comment = "\#", col_names =
FALSE, show_col_types = FALSE)

canon_names \<- c( "query", "seed_ortholog", "evalue", "score",
"eggNOG_OGs", "best_tax_level", "COG_category", "Description",
"Preferred_name", "GOs", "EC", "KOs", "KEGG_Pathway", "KEGG_Module",
"KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy",
"BiGG_Reaction", "PFAMs" )

names(emap) \<- canon_names[seq_len(ncol(emap))] attr(emap, "query_col")
\<- "query" emap }

emap \<- read_emapper_annotations(emapper_path) query_col \<- attr(emap,
"query_col")

emap_small \<- emap \|\> dplyr::transmute( protein_id =
.data[[query_col]], EggNM.Preferred_name =
dplyr::coalesce(.data[["Preferred_name"]], NA_character\_), GOs = if
("GOs" %in% names(.)) .data[["GOs"]] else NA_character\_ ) \|\>
dplyr::filter(!is.na(protein_id), nzchar(protein_id)) \|\>
dplyr::distinct()

trinotate_for_echogo \<- gene_to_prot \|\> dplyr::inner_join(emap_small,
by = "protein_id") \|\> dplyr::group_by(gene_id) \|\> dplyr::summarise(
transcript_id = dplyr::first(gene_id), sprot_Top_BLASTX_hit =
NA_character\_, EggNM.Preferred_name = { vals \<-
EggNM.Preferred_name[!is.na(EggNM.Preferred_name)] if (length(vals))
vals[1] else NA_character\_ }, EggNM.max_annot_lvl = 33208L, .groups =
"drop" ) \|\> dplyr::select(transcript_id, sprot_Top_BLASTX_hit,
EggNM.Preferred_name, EggNM.max_annot_lvl)

# Write BOTH final annotation files (EchoGO bundle)

trinotate_path1 \<- file.path(echogo_input, "Trinotate_for_EchoGO.tsv")
trinotate_path2 \<- file.path(echogo_input, paste0(reference_label,
"\_eggNOG_for_EchoGO.tsv"))

readr::write_tsv(trinotate_for_echogo, trinotate_path1)
readr::write_tsv(trinotate_for_echogo, trinotate_path2)

message(" ✓ Wrote: ", trinotate_path1) message(" ✓ Wrote: ",
trinotate_path2) message(" ✓ Rows: ", nrow(trinotate_for_echogo), "\n")

# ==============================================================================

# STEP 9: Build gene → GO mapping from eggNOG (internal for GOseq)

# ==============================================================================

message("STEP 9: Building gene → GO mapping for GOseq...")

gene_go_tbl \<- gene_to_prot \|\> dplyr::inner_join(emap_small \|\>
dplyr::select(protein_id, GOs), by = "protein_id") \|\>
dplyr::mutate(GOs = dplyr::coalesce(GOs, NA_character\_)) \|\>
dplyr::filter(!is.na(GOs), nzchar(GOs)) \|\> dplyr::distinct(gene_id,
GOs)

GO_info_listed \<- gene_go_tbl \|\> dplyr::mutate(GO_list =
strsplit(GOs, ",")) \|\> dplyr::mutate(GO_list = purrr::map(GO_list, \~
trimws(.x))) \|\> dplyr::mutate(GO_list = purrr::map(GO_list, \~
.x[nzchar(.x)])) \|\> { stats::setNames(.$GO_list, .$gene_id) }

GO_to_gene_list \<- list() for (gid in intersect(names(GO_info_listed),
bg_universe)) { go_list \<- GO_info_listed[[gid]] for (go_id in go_list)
{ GO_to_gene_list[[go_id]] \<- c(GO_to_gene_list[[go_id]], gid) } }

message(" ✓ Genes with GO annotations: ", length(GO_info_listed), "\n")

# ==============================================================================

# STEP 10: Run GOseq per contrast and write final GOseq tables into echogo_input/

# ==============================================================================

message("STEP 10: Running GOseq per contrast...")

gene_lengths_tbl \<- readr::read_tsv(gene_len_file, show_col_types =
FALSE) \|\> dplyr::transmute( gene_id = strip_gene(.data[["gene_id"]]),
length = as.numeric(.data[["length"]]) ) \|\>
dplyr::filter(!is.na(gene_id), nzchar(gene_id), !is.na(length), length
\> 0) \|\> dplyr::distinct(gene_id, .keep_all = TRUE)

bg_with_len \<- intersect(bg_universe, gene_lengths_tbl\$gene_id) if
(!length(bg_with_len)) stop("No overlap between BG universe and gene
lengths table.")

gene_lengths_vec \<- gene_lengths_tbl$length
names(gene_lengths_vec) <- gene_lengths_tbl$gene_id

sample_set_gene_ids \<- bg_with_len gene_lengths_mat \<-
as.matrix(gene_lengths_vec[sample_set_gene_ids])
rownames(gene_lengths_mat) \<- sample_set_gene_ids

get_GO_term_descr \<- function(x) { go_info \<- GO.db::GOTERM[[x]] if
(is.null(go_info) \|\| length(go_info) == 0) return("none")
paste(AnnotationDbi::Ontology(go_info), AnnotationDbi::Term(go_info)) }

de_list_files \<- list.files(lists_dir, pattern = "\_DE_genes\\.txt\$",
full.names = TRUE)

goseq_summary \<- purrr::map_dfr(de_list_files, function(de_path) { base
\<- sub("\_DE_genes\\.txt\$", "", basename(de_path)) message(" - GOseq
for: ", base)

de_genes \<- readLines(de_path) de_genes \<- strip_gene(de_genes)
de_genes \<- unique(nzc(de_genes))

de_in_universe \<- intersect(de_genes, sample_set_gene_ids)

cat_genes_vec \<- as.integer(sample_set_gene_ids %in% de_in_universe)
names(cat_genes_vec) \<- sample_set_gene_ids

pwf \<- goseq::nullp(DEgenes = cat_genes_vec, bias.data =
gene_lengths_mat) rownames(pwf) \<- sample_set_gene_ids

res \<- goseq::goseq(pwf, gene2cat = GO_info_listed,
use_genes_without_cat = TRUE)

p_over \<- res$over_represented_pvalue
  p_over[p_over > 1 - 1e-10] <- 1 - 1e-10
  res$over_represented_FDR \<- qvalue::qvalue(p_over)\$qvalues

p_under \<- res$under_represented_pvalue
  p_under[p_under > 1 - 1e-10] <- 1 - 1e-10
  res$under_represented_FDR \<- qvalue::qvalue(p_under)\$qvalues

res_enriched \<- res \|\> dplyr::filter(over_represented_pvalue \<=
0.05) \|\> dplyr::mutate( go_term = vapply(category, get_GO_term_descr,
character(1)), gene_ids = vapply(category, function(go_id) { gene_list
\<- GO_to_gene_list[[go_id]] if (is.null(gene_list)) return("")
gene_list \<- unique(strip_gene(gene_list)) gene_list \<-
gene_list[gene_list %in% de_in_universe] paste(gene_list, collapse = ",
") }, character(1)) ) \|\> dplyr::arrange(over_represented_pvalue)

res_depleted \<- res \|\> dplyr::filter(under_represented_pvalue \<=
0.05) \|\> dplyr::mutate( go_term = vapply(category, get_GO_term_descr,
character(1)), gene_ids = vapply(category, function(go_id) { gene_list
\<- GO_to_gene_list[[go_id]] if (is.null(gene_list)) return("")
gene_list \<- unique(strip_gene(gene_list)) gene_list \<-
gene_list[gene_list %in% de_in_universe] paste(gene_list, collapse = ",
") }, character(1)) ) \|\> dplyr::arrange(under_represented_pvalue)

go_enrich_file \<- file.path(echogo_input, paste0(base,
".GOseq.enriched.tsv")) go_depleted_file \<- file.path(echogo_input,
paste0(base, ".GOseq.depleted.tsv"))

readr::write_tsv(res_enriched, go_enrich_file)

\# Only write depleted table if there are any rows (keeps folder clean)
if (nrow(res_depleted) \> 0) { readr::write_tsv(res_depleted,
go_depleted_file) }

tibble::tibble( contrast = base, n_bg_genes =
length(sample_set_gene_ids), n_de_genes = length(de_in_universe),
n_GO_enriched = nrow(res_enriched), n_GO_depleted = nrow(res_depleted) )
})

# Keep GOseq summary in echogo_work (not in final input bundle)

readr::write_csv(goseq_summary, file.path(goseq_dir,
"summary_goseq_per_contrast.csv")) print(goseq_summary) cat("\n")

# ==============================================================================

# STEP 11: Final sanity check — ensure final input folder matches expectation

# ==============================================================================

message("STEP 11: Final sanity check of EchoGO input folder...")

final_files \<- sort(list.files(echogo_input))

message("Final EchoGO input folder:\n ", echogo_input, "\n")
message("Files present:") print(final_files)

missing \<- character(0)

if (!file.exists(file.path(echogo_input, "allcounts_table.txt"))) {
missing \<- c(missing, "allcounts_table.txt") } if
(!any(grepl("\^dge\_.\*\\.csv$", final_files))) {
  missing <- c(missing, "dge_*.csv")
}
if (!any(grepl("\\.GOseq\\.enriched\\.tsv$", final_files))) { missing
\<- c(missing, "dge\_<CONTRAST>.GOseq.enriched.tsv") } if
(!file.exists(file.path(echogo_input, "Trinotate_for_EchoGO.tsv"))) {
missing \<- c(missing, "Trinotate_for_EchoGO.tsv") } if
(!file.exists(file.path(echogo_input, paste0(reference_label,
"\_eggNOG_for_EchoGO.tsv")))) { missing \<- c(missing,
paste0(reference_label, "\_eggNOG_for_EchoGO.tsv")) }

if (length(missing)) { message("\n⚠ Missing expected file types:\n  - ",
paste(missing, collapse = "\n  - ")) } else { message("\n✓ EchoGO input
bundle looks complete.") }

message("\n==============================================================================")
message("DONE — You can now run EchoGO using:") message(" input_dir =
"", echogo_input, """)
message("==============================================================================")

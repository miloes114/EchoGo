[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17658715.svg)](https://doi.org/10.5281/zenodo.17658715)

# 🧬 EchoGO

### Cross-Species Consensus Enrichment for Non-Model Organisms

**EchoGO** is a modular, end-to-end functional enrichment pipeline designed for **de novo transcriptomes**
(Trinity → Trinotate → GOseq) and **reference-based RNA-seq workflows** (e.g., HISAT2/STAR + featureCounts + DESeq2),
as long as you provide the expected input bundle.

It provides:

-   Independent **g:Profiler organism-context analyses** using one matched resolved gene set
-   Integration of **GOseq** and **g:Profiler** into **conservative** and **exploratory** results
-   Semantic similarity reduction using **RRvGO**
-   GO term overlap **network analysis**
-   Optional fully rendered **HTML reports**
------------------------------------------------------------------------

## 📥 Installation

### 1) Linux users — install required system libraries

(Windows and macOS users can skip this step.)

``` bash
sudo apt-get update
sudo apt-get install -y \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libfontconfig1-dev libfreetype6-dev libharfbuzz-dev \
  libfribidi-dev libpng-dev libtiff5-dev libjpeg-dev
```

### 2) Install base R helpers

``` r
install.packages(c("remotes", "BiocManager"))
```

### 3) Install EchoGO

``` r
remotes::install_github("miloes114/EchoGo", build_vignettes = TRUE)
library(EchoGO)
echogo_help()
```

For a faster installation without locally registered vignettes, set

`build_vignettes = FALSE`.
**Windows note**: If installation warns about curl.dll “Permission denied”, close all R/RStudio sessions and delete the 00LOCK folder in your user library, then reinstall curl and EchoGO.

If you want `browseVignettes("EchoGO")` and `vignette("EchoGO_workflow")` to work after installation, install EchoGO from a **built package tarball**:

- `remotes::install_github(..., build_vignettes = TRUE)` works because it builds the package first.
- A release asset such as `EchoGO_0.1.3.tar.gz` also works.
- A raw repository snapshot ZIP from GitHub/Zenodo is **not** the same thing as a built R package and does **not** include the installed vignette metadata that `browseVignettes()` uses.

In short: if you install directly from a source snapshot folder with `R CMD INSTALL <folder>`, the package can still work, but the vignettes will usually not be registered.

### 4) Install the OrgDb needed for RRvGO

Recommended automatic method:

``` r
EchoGO::echogo_install_orgdb("org.Mm.eg.db")
```

`GO.db` and `rrvgo` are installed with EchoGO. An organism-specific OrgDb,
such as `org.Mm.eg.db`, is still needed when RRvGO semantic reduction runs.

### 5) Optional: Dependencies for HTML report generation

``` r
install.packages(c(
  "rmarkdown", "knitr", "DT", "gt",
  "patchwork", "ggforce", "plotly",
  "visNetwork", "ggtext"
))
```

HTML reports and vignettes require Pandoc. RStudio includes it; when using R
from another terminal, confirm that `rmarkdown::pandoc_available()` returns
`TRUE` and install Quarto or Pandoc if it does not.

### 6) Test your installation

``` r
library(EchoGO)
EchoGO::echogo_quickstart(run_demo = TRUE)
```

This is a deterministic offline validation using cached demonstration responses
and the conservative/background-aware stream. It does not need an internet
connection. Use `echogo_quickstart(run_demo = FALSE)` to copy and inspect only
the inputs. Use `echogo_quickstart(run_demo = TRUE, live_gprofiler = TRUE)` for
an optional live integration test. To exercise every exploratory, RRvGO, and evaluation stage,
run `echogo_quickstart(run_demo = TRUE, full = TRUE)`; the full run can take
several minutes.

------------------------------------------------------------------------

## 🏁 Quickstart

``` r
# 1) Install one annotation database for RRvGO
EchoGO::echogo_install_orgdb("org.Mm.eg.db")

# 2) Validate species
species <- echogo_preflight_species(c("hsapiens", "mmusculus", "drerio"))

# 3) Create project scaffold
echogo_scaffold("my_project")

# 4) Place your input files into my_project/input/
#    De novo mode:
#      - gene.counts.matrix.tsv
#      - DE_*.tsv
#      - Trinotate.xls
#    Reference-based mode:
#      - allcounts_table.txt
#      - dge_<CONTRAST>.csv
#      - dge_<CONTRAST>.GOseq.enriched.tsv
#      - Trinotate_for_EchoGO.tsv
#      - <reference_label>_eggNOG_for_EchoGO.tsv
#    Generate these with echogo_prepare_reference_inputs().
#    See: doc/reference-based-inputs.html

# 5) Run the pipeline
echogo_run("my_project/input", "my_project/results")

# 6) Inspect consensus tables & plots in my_project/results/
```

------------------------------------------------------------------------

## 📚 Documentation

### Online (GitHub)
- Workflow vignette: [EchoGO_workflow](doc/EchoGO_workflow.html)
- Interpretation guide: [EchoGO_interpretation](doc/EchoGO_interpretation.html)
- Reference-based RNA-seq inputs: [reference-based-inputs](doc/reference-based-inputs.html)

### In R (after installation)
```r
browseVignettes("EchoGO")
vignette("EchoGO_workflow")
vignette("EchoGO_interpretation")
vignette("reference-based-inputs")
```

These commands require an installation that includes built vignette metadata. If they do not appear after install, reinstall from GitHub with:

```r
remotes::install_github("miloes114/EchoGo", build_vignettes = TRUE, force = TRUE)
```
------------------------------------------------------------------------

## ⚙️ Configuration

### Set default annotation database

``` r
options(EchoGO.default_orgdb = "org.Mm.eg.db")
```

For multi-species enrichment:

``` r
options(EchoGO.default_orgdb = c("org.Hs.eg.db", "org.Mm.eg.db", "org.Dr.eg.db"))
```

Check which OrgDb packages are installed:

``` r
EchoGO::echogo_list_orgdb()
```

------------------------------------------------------------------------

## 🔎 Species Selection

EchoGO provides tools to browse, search, validate, and programmatically select supported g:Profiler species.

### List supported organisms (interactive if DT is available)

``` r
echogo_list_species(view = TRUE)
```

### Search by name or NCBI ID

``` r
echogo_species_lookup(c("human", "mouse", "9606"))
```

### Fuzzy guessing

``` r
echogo_species_lookup(c("H. sapiens", "rat", "zebra fish"))
```

### Validate selected IDs

``` r
species <- echogo_preflight_species(c("hsapiens", "mmusculus", "drerio"))
```

### Filter by taxonomy & tags

``` r
# Interactive browsing
echogo_list_species(view = TRUE)

# Programmatic filter
ids <- echogo_resolve("tag:AnimalModels OR order:Perciformes")
```

------------------------------------------------------------------------

## 📂 Expected Inputs

If you create a scaffold:

``` r
echogo_scaffold("my_project")
```

Place these files under `my_project/input/`:

| File | Description |
|------------------------|------------------------------------------------|
| `gene.counts.matrix.tsv` | TSV format: first column = transcript/gene ID, remaining columns = sample counts |
| `DE_*.tsv` | One or more DE tables with columns: `id`, `log2FC`, `pvalue`, `padj` |
| `Trinotate.xls` | Standard Trinotate report with GO & KEGG annotation |

These correspond to the standard Trinity/Trinotate/GOseq workflow:\
<https://github.com/trinityrnaseq/trinityrnaseq/wiki>

### Reference-based RNA-seq (HISAT2/STAR + featureCounts + DESeq2)

EchoGO can build the reference-based bundle from a count matrix, DESeq2 CSV
files, the matching GFF3/GTF, and the matching protein FASTA. The preparation
helper is resumable around the external eggNOG-mapper step:

``` r
prep <- EchoGO::echogo_prepare_reference_inputs(
  root = "path/to/reference_project",
  gff_file = "genome.gff3.gz",
  protein_fasta = "proteins.faa.gz",
  reference_label = "MySpecies",
  orgdb = NULL
)

# Run the command stored in prep$emapper_command, then repeat the same call.
stopifnot(prep$status %in% c("awaiting_emapper", "complete"))
```

### Scientific gene-set contract (v0.1.3)

For each contrast, EchoGO derives the g:Profiler foreground from significant
genes in the DE table and the custom background from the genes tested in that
same experiment. The default threshold rule is `padj <= 0.05` and
`abs(log2FoldChange) >= 1`; an explicit logical `significant` column takes
precedence when present. Missing adjusted p-values do not pass the threshold.

The tested-universe precedence is: an explicit tested-gene list, all rows of a
full DE table, then the matching count matrix when the DE file is explicitly
declared significant-only. Foreground and background pass through the same
deterministic portable canonical-name resolver. Its fixed priority is a genuine
SwissProt-derived gene symbol, then `EggNM.Preferred_name`, then a portable
native symbol. Raw transcript, contig, seed-ortholog, species-prefixed Ensembl
protein, and arbitrary accession IDs are never submitted to g:Profiler; they
remain in `mapping_table.csv` with the reason `no portable canonical name`.
Blank values are ignored and duplicate portable names collapse by stable first
occurrence. Taxonomy labels are retained as provenance but do not numerically
filter this resolver. After mapping,
the foreground must be a proper non-degenerate subset of the background or the
run stops with a diagnostic.

The same resolved vectors are submitted independently under every selected
g:Profiler organism code. This is organism-context enrichment in a shared
identifier space, not construction of species-specific ortholog gene sets.
Exact vectors, the mapping table, per-run response metadata, effective sizes,
and hashes are saved under `gprofiler/`.

Use a species-native OrgDb, such as `org.Dr.eg.db`, in `orgdb` when one is
available. Otherwise, leave it `NULL` and GOseq uses eggNOG GO terms. The full
two-stage walkthrough and troubleshooting guide is here:

-   **Reference-based input preparation tutorial:** [reference-based-inputs](doc/reference-based-inputs.html)

Place these files under `my_project/input/`:

| File | Description |
|------------------------|------------------------------------------------|
| **allcounts_table.txt** | featureCounts gene count matrix (first column = gene ID, remaining columns = sample counts) |
| `dge_<CONTRAST>.csv` | One DESeq2 results table per contrast (e.g., `dge_Treatment_vs_Control.csv`) |
| `dge_<CONTRAST>.GOseq.enriched.tsv` | GOseq terms passing adjusted FDR 0.05 |
| `dge_<CONTRAST>.GOseq.depleted.tsv` | Optional: GOseq depleted table from a custom preparation workflow |
| **Trinotate_for_EchoGO.tsv** | Minimal Trinotate-like table used by EchoGO (created from eggNOG output) |
| `<reference_label>_eggNOG_for_EchoGO.tsv` | Species-labeled eggNOG table used by EchoGO (created from eggNOG output) |

At the end, your input folder should look like this:

```   text
my_project/input/
├── allcounts_table.txt
├── dge_<CONTRAST>.csv
├── dge_<CONTRAST>.GOseq.enriched.tsv
├── dge_<CONTRAST>.GOseq.depleted.tsv        (optional)
├── Trinotate_for_EchoGO.tsv
├── <reference_label>_eggNOG_for_EchoGO.tsv
```  

## ▶️ Running the Pipeline

### Basic scaffolded workflow
``` r
echogo_scaffold("my_project")

echogo_run(
  input_dir = "my_project/input",
  outdir = "my_project/results"
)
```

### Full explicit call

``` r
run_full_echogo(
  input_dir = "my_project/input",
  goseq_file = "my_project/input/dge_Treatment_vs_Control.GOseq.enriched.tsv",
  trinotate_file = "my_project/input/Trinotate_for_EchoGO.tsv",
  de_file = "my_project/input/dge_Treatment_vs_Control.csv",
  count_matrix_file = "my_project/input/allcounts_table.txt",
  species   = species,      # or species_expr = "..."
  orgdb     = "org.Mm.eg.db",
  outdir    = "my_project/results",
  make_report = TRUE,
  verbose     = TRUE
)
```

### Using species expression

``` r
run_full_echogo(
  input_dir = "my_project/input",
  species_expr = "tag:AnimalModels OR order:Perciformes",
  outdir = "my_project/results",
  make_report = TRUE
)
```

### Multiple OrgDb databases

``` r
run_full_echogo(
  input_dir = "my_project/input",
  species = species,
  orgdb   = c("org.Mm.eg.db", "org.Dr.eg.db"),
  outdir  = "my_project/results"
)
```

RRvGO and network modules will automatically use whichever annotation database is available.

------------------------------------------------------------------------

## 📤 Outputs

After running, the results directory includes:

```         
my_project/results/
├── consensus/
│   ├── consensus_enrichment_results_with_and_without_bg.xlsx
│   ├── plots_strict/
│   └── plots_exploratory/
├── diagnostics/
├── evaluation/
├── goseq/
├── gprofiler/
│   ├── run_manifest.json
│   ├── submitted_vectors/
│   ├── with_custom_background/
│   └── no_background_genome_wide/
├── rrvgo/
│   ├── rrvgo_true_consensus_with_bg/
│   └── rrvgo_exploratory_all_significant/
├── networks/
│   ├── with_bg/
│   └── with_bg_and_nobg/
└── report/
    └── (HTML report if make_report = TRUE)
```

`rrvgo_true_consensus_with_bg/` is a compatibility-only legacy folder name;
in v0.1.3 it contains the background-aware RRvGO stream.

------------------------------------------------------------------------

## 📘 Demo Datasets

EchoGO ships with frozen demo inputs and results.

``` r
# View demo paths
echogo_demo_path()            # demo inputs
echogo_demo_results_path()    # demo results

# Run small demo pipeline
echogo_quickstart(run_demo = TRUE)

# Open demo folders
echogo_open_demo()
```

------------------------------------------------------------------------

## 📑 Citation

If you use EchoGO in publications, please cite:

> Escobar-Sierra C., Langschied F., Inostroza P.A. (2026).\
> **EchoGO: A Cross-Species Consensus Framework for Functional Enrichment in Non-Model Organisms (v0.1.3).**\
> Zenodo. DOI: <https://doi.org/10.5281/zenodo.17658715>

Full citation entry is included in `inst/CITATION`.

------------------------------------------------------------------------

## 🐛 Issues & Support

Please report bugs, suggestions, or feature requests at:\
[**https://github.com/miloes114/EchoGo/issues**](https://github.com/miloes114/EchoGo/issues){.uri}

For general questions, contact the maintainers.

------------------------------------------------------------------------

## 💬 Acknowledgements

EchoGO was developed by **Camilo Escobar-Sierra**, **Felix Langschied**, and **Pedro A. Inostroza**, with additional input from collaborators and the community.

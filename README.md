[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17660708.svg)](https://doi.org/10.5281/zenodo.17660708)

# 🧬 EchoGO

### Cross-Species Consensus Enrichment for Non-Model Organisms

**EchoGO** is a modular, end-to-end functional enrichment pipeline designed for:

-   De novo transcriptomes and non-model species\
-   Multi-species orthology support via **g:Profiler**\
-   Integration of **GOseq** and **g:Profiler** into strict & exploratory consensus\
-   Semantic reduction using **RRvGO**\
-   GO term overlap **network analysis**\
-   Optional fully rendered **HTML reports**

EchoGO supports both **de novo transcriptomes (Trinity → Trinotate → GOseq)** and **reference-based RNA-seq workflows** (e.g., HISAT2/STAR + featureCounts + DESeq2), as long as you provide the expected input bundle.

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
remotes::install_github("miloes114/EchoGo", build_vignettes = FALSE)
library(EchoGO)
echogo_help()
```

To install with vignettes:

``` r
remotes::install_github("miloes114/EchoGo", build_vignettes = TRUE)
```

### 4) Install annotation packages (GO.db + OrgDb)

Recommended automatic method:

``` r
EchoGO::echogo_install_orgdb()
```

Or install explicitly:

``` r
EchoGO::echogo_install_orgdb(c("GO.db", "org.Mm.eg.db"))
```

EchoGO requires **GO.db** and at least one **OrgDb** package (e.g., `org.Mm.eg.db`).

### 5) Optional but recommended: RRvGO semantic reduction

``` r
BiocManager::install("rrvgo")
```

### 6) Optional: Dependencies for HTML report generation

``` r
install.packages(c(
  "rmarkdown", "knitr", "DT", "gt",
  "patchwork", "ggforce", "plotly",
  "visNetwork", "ggtext"
))
```

### 7) Test your installation

``` r
echogo_quickstart(run_demo = TRUE)
```

This confirms that consensus scoring, RRvGO, networks, and the HTML report all run successfully.

------------------------------------------------------------------------

## 🏁 Quickstart

``` r
# 1) Install annotation databases
EchoGO::echogo_install_orgdb()

# 2) Validate species
species <- echogo_validate_species(c("hsapiens", "mmusculus", "drerio"))

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
#    See: doc/reference-based-inputs.md

# 5) Run the pipeline
echogo_run("my_project/input", "my_project/results")

# 6) Inspect consensus tables & plots in my_project/results/
```

------------------------------------------------------------------------

## 📚 Documentation

### Online (GitHub)
- Workflow vignette: [EchoGO_workflow](doc/EchoGO_workflow.html)
- Interpretation guide: [EchoGO_interpretation](doc/EchoGO_interpretation.html)
- Reference-based RNA-seq inputs: [reference-based-inputs](doc/reference-based-inputs.md)

### In R (after installation)
```r
vignette("EchoGO_workflow")
vignette("EchoGO_interpretation")


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
echogo_guess_species(c("H. sapiens", "rat", "zebra fish"))
```

### Validate selected IDs

``` r
species <- echogo_validate_species(c("hsapiens", "mmusculus", "drerio"))
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
| **gene.counts.matrix.tsv** | TSV format: first column = transcript/gene ID, remaining columns = sample counts |
| \*\*DE\_\*.tsv\*\* | One or more DE tables with columns: `id`, `log2FC`, `pvalue`, `padj` |
| **Trinotate.xls** | Standard Trinotate report with GO & KEGG annotation |

These correspond to the standard Trinity/Trinotate/GOseq workflow:\
<https://github.com/trinityrnaseq/trinityrnaseq/wiki>

### Reference-based RNA-seq (HISAT2/STAR + featureCounts + DESeq2)

EchoGO can also run reference-based experiments, as long as you prepare an input folder that mirrors the expected bundle. A complete step-by-step guide (including the external eggNOG-mapper step) is here:

-   **Reference-based input preparation tutorial:** `doc/reference-based-inputs.md`

Place these files under `my_project/input/`:

| File | Description |
|------------------------|------------------------------------------------|
| **allcounts_table.txt** | featureCounts gene count matrix (first column = gene ID, remaining columns = sample counts) |
| **dge\_**<CONTRAST>.csv | One DESeq2 results table per contrast (e.g., `dge_Treatment_vs_Control.csv`) |
| **dge\_**<CONTRAST>.GOseq.enriched.tsv | GOseq enriched table produced during input preparation |
| **dge\_**<CONTRAST>.GOseq.depleted.tsv | Optional: GOseq depleted table (only if any depleted terms exist) |
| **Trinotate_for_EchoGO.tsv** | Minimal Trinotate-like table used by EchoGO (created from eggNOG output) |
| \*\*<reference_label>\_eggNOG_for_EchoGO.tsv\*\* | Species-labeled eggNOG table used by EchoGO (created from eggNOG output) |

At the end, your input folder should look like this:

``` text
my_project/input/
├── allcounts_table.txt
├── dge_<CONTRAST>.csv
├── dge_<CONTRAST>.GOseq.enriched.tsv
├── dge_<CONTRAST>.GOseq.depleted.tsv        (optional)
├── Trinotate_for_EchoGO.tsv
├── <reference_label>_eggNOG_for_EchoGO.tsv

---

## ▶️ Running the Pipeline

### Basic scaffolded workflow
```r
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

> Escobar-Sierra C., Langschied F., Inostroza P.A. (2025).\
> **EchoGO: Cross-Species Consensus Functional Enrichment Analysis.**\
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

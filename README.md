[![DOI](https://zenodo.org/badge/803971735.svg)](https://doi.org/10.5281/zenodo.17660708)

# 🧬 EchoGO  
### Cross-Species Consensus Enrichment for Non-Model Organisms

**EchoGO** is a modular, end-to-end functional enrichment pipeline designed for:

- De novo transcriptomes and non-model species  
- Multi-species orthology support via **g:Profiler**  
- Integration of **GOseq** and **g:Profiler** into strict & exploratory consensus  
- Semantic reduction using **RRvGO**  
- GO term overlap **network analysis**  
- Optional fully rendered **HTML reports**

EchoGO expects that users have already run a standard **Trinity → Trinotate → GOseq** workflow or an equivalent RNA-seq pipeline.

---

## 📥 Installation

### 1) Linux users — install required system libraries
(Windows and macOS users can skip this step.)
```bash
sudo apt-get update
sudo apt-get install -y \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libfontconfig1-dev libfreetype6-dev libharfbuzz-dev \
  libfribidi-dev libpng-dev libtiff5-dev libjpeg-dev
```

### 2) Install base R helpers
```r
install.packages(c("remotes", "BiocManager"))
```

### 3) Install EchoGO
```r
remotes::install_github("miloes114/EchoGo", build_vignettes = FALSE)
library(EchoGO)
echogo_help()
```

To install with vignettes:
```r
remotes::install_github("miloes114/EchoGo", build_vignettes = TRUE)
```

### 4) Install annotation packages (GO.db + OrgDb)

Recommended automatic method:
```r
EchoGO::echogo_install_orgdb()
```

Or install explicitly:
```r
EchoGO::echogo_install_orgdb(c("GO.db", "org.Mm.eg.db"))
```

EchoGO requires **GO.db** and at least one **OrgDb** package (e.g., `org.Mm.eg.db`).

### 5) Optional but recommended: RRvGO semantic reduction
```r
BiocManager::install("rrvgo")
```

### 6) Optional: Dependencies for HTML report generation
```r
install.packages(c(
  "rmarkdown", "knitr", "DT", "gt",
  "patchwork", "ggforce", "plotly",
  "visNetwork", "ggtext"
))
```

### 7) Test your installation
```r
echogo_quickstart(run_demo = TRUE)
```

This confirms that consensus scoring, RRvGO, networks, and the HTML report all run successfully.

---

## 🏁 Quickstart
```r
# 1) Install annotation databases
EchoGO::echogo_install_orgdb()

# 2) Validate species
species <- echogo_validate_species(c("hsapiens", "mmusculus", "drerio"))

# 3) Create project scaffold
echogo_scaffold("my_project")

# 4) Place your input files into my_project/input/
#    - gene.counts.matrix.tsv
#    - DE_*.tsv (differential expression results)
#    - Trinotate.xls

# 5) Run the pipeline
echogo_run("my_project/input", "my_project/results")

# 6) Inspect consensus tables & plots in my_project/results/
```

---

## ⚙️ Configuration

### Set default annotation database
```r
options(EchoGO.default_orgdb = "org.Mm.eg.db")
```

For multi-species enrichment:
```r
options(EchoGO.default_orgdb = c("org.Hs.eg.db", "org.Mm.eg.db", "org.Dr.eg.db"))
```

Check which OrgDb packages are installed:
```r
EchoGO::echogo_list_orgdb()
```

---

## 🔎 Species Selection

EchoGO provides tools to browse, search, validate, and programmatically select supported g:Profiler species.

### List supported organisms (interactive if DT is available)
```r
echogo_list_species(view = TRUE)
```

### Search by name or NCBI ID
```r
echogo_species_lookup(c("human", "mouse", "9606"))
```

### Fuzzy guessing
```r
echogo_guess_species(c("H. sapiens", "rat", "zebra fish"))
```

### Validate selected IDs
```r
species <- echogo_validate_species(c("hsapiens", "mmusculus", "drerio"))
```

### Filter by taxonomy & tags
```r
# Interactive browsing
echogo_list_species(view = TRUE)

# Programmatic filter
ids <- echogo_resolve("tag:AnimalModels OR order:Perciformes")
```

---

## 📂 Expected Inputs

If you create a scaffold:
```r
echogo_scaffold("my_project")
```

Place these files under `my_project/input/`:

| File | Description |
|------|-------------|
| **gene.counts.matrix.tsv** | TSV format: first column = transcript/gene ID, remaining columns = sample counts |
| **DE_*.tsv** | One or more DE tables with columns: `id`, `log2FC`, `pvalue`, `padj` |
| **Trinotate.xls** | Standard Trinotate report with GO & KEGG annotation |

These correspond to the standard Trinity/Trinotate/GOseq workflow:  
https://github.com/trinityrnaseq/trinityrnaseq/wiki

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
```r
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
```r
run_full_echogo(
  input_dir = "my_project/input",
  species_expr = "tag:AnimalModels OR order:Perciformes",
  outdir = "my_project/results",
  make_report = TRUE
)
```

### Multiple OrgDb databases
```r
run_full_echogo(
  input_dir = "my_project/input",
  species = species,
  orgdb   = c("org.Mm.eg.db", "org.Dr.eg.db"),
  outdir  = "my_project/results"
)
```

RRvGO and network modules will automatically use whichever annotation database is available.

---

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

---

## 📘 Demo Datasets

EchoGO ships with frozen demo inputs and results.
```r
# View demo paths
echogo_demo_path()            # demo inputs
echogo_demo_results_path()    # demo results

# Run small demo pipeline
echogo_quickstart(run_demo = TRUE)

# Open demo folders
echogo_open_demo()
```

---

## 📑 Citation

If you use EchoGO in publications, please cite:

> Escobar-Sierra C., Langschied F., Inostroza P.A. (2025).  
> **EchoGO: Cross-Species Consensus Functional Enrichment Analysis.**  
> Zenodo. DOI: https://doi.org/10.5281/zenodo.17658715

Full citation entry is included in `inst/CITATION`.

---

## 🐛 Issues & Support

Please report bugs, suggestions, or feature requests at:  
**https://github.com/miloes114/EchoGo/issues**

For general questions, contact the maintainers.

---

## 💬 Acknowledgements

EchoGO was developed by **Camilo Escobar-Sierra**, **Felix Langschied**, and **Pedro A. Inostroza**, with additional input from collaborators and the community.

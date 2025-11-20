[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17658715.svg)](https://doi.org/10.5281/zenodo.17658715)

# 🧬 EchoGO  
### Cross-Species Consensus Enrichment for Non-Model Organisms

**EchoGO** is a modular, end-to-end functional enrichment pipeline designed for:

- de novo transcriptomes and non-model species  
- multi-species orthology support via **g:Profiler**  
- integration of **GOseq** and **g:Profiler** into strict & exploratory consensus  
- semantic reduction using **RRvGO**  
- GO term overlap **network analysis**  
- optional fully rendered **HTML reports**

EchoGO expects that users have already run a standard  
**Trinity → Trinotate → GOseq** workflow or an equivalent RNA-seq pipeline.

From those inputs, EchoGO produces:

- bias-aware GOseq enrichment summaries  
- cross-species (background-corrected + exploratory) enrichment  
- strict and exploratory consensus scoring tables  
- RRvGO semantic similarity clusters  
- GO term overlap networks (static and interactive)  
- optional HTML reports

---

# 📥 Installation

```r
# Install annotation packages (GO.db + OrgDb)
EchoGO::echogo_install_orgdb()

# Or specify explicitly:
EchoGO::echogo_install_orgdb(c("GO.db", "org.Mm.eg.db"))
```

EchoGO requires R and at least one **OrgDb** package (e.g., `org.Mm.eg.db`).

---

# ⚙️ Setting default annotation database

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

# 🔎 Species exploration & selection

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

### Use validated species in the pipeline
```r
run_full_echogo(
  input_dir = "my_project/input",
  species = species,
  outdir = "my_project/results",
  make_report = TRUE
)
```

---

# 🧬 Species selection by taxonomy & tags

EchoGO includes taxonomy and high-level “tags” for programmatic species selection.

### Interactive browsing
```r
echogo_list_species(view = TRUE)
```

### Programmatic filter
```r
ids <- echogo_resolve("tag:AnimalModels OR order:Perciformes")
```

Use directly in the pipeline:

```r
run_full_echogo(
  input_dir = "my_project/input",
  species_expr = "tag:AnimalModels OR order:Perciformes",
  make_report = FALSE
)
```

---

# 📂 Expected inputs (`input/`)

If you create a scaffold:

```r
echogo_scaffold("my_project")
```

Then place these files under `my_project/input/`:

| File | Description |
|------|-------------|
| **gene.counts.matrix.tsv** | TSV, first column = transcript/gene ID, remaining columns = sample counts |
| **DE_*.tsv** | One or more DE tables with columns `id`, `log2FC`, `pvalue`, `padj` |
| **Trinotate.xls** | Standard Trinotate report with GO & KEGG annotation |

These correspond to the standard Trinity/Trinotate/GOseq workflow:  
https://github.com/trinityrnaseq/trinityrnaseq/wiki

---

# ▶️ Running the pipeline

### Scaffolded project workflow
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

# 📤 Outputs

After running, the results directory includes:

```
consensus/
    consensus_enrichment_results_with_and_without_bg.xlsx
    plots_strict/
    plots_exploratory/

diagnostics/
evaluation/

goseq/
gprofiler/
    with_custom_background/
    no_background_genome_wide/

rrvgo/
    rrvgo_true_consensus_with_bg/
    rrvgo_exploratory_all_significant/

networks/
    with_bg/
    with_bg_and_nobg/

report/
    (HTML report if make_report = TRUE)
```

---

# 📘 Demo datasets

EchoGO ships with frozen demo inputs and results.

```r
echogo_demo_path()            # demo inputs
echogo_demo_results_path()    # demo results

echogo_quickstart(run_demo = TRUE)  # run small demo pipeline
echogo_open_demo()                 # open folders
```

---

# 🏁 Minimal quickstart

```r
# 1) Install annotation DBs
EchoGO::echogo_install_orgdb()

# 2) Choose species
species <- echogo_validate_species(c("hsapiens","mmusculus","drerio"))

# 3) Create scaffold
echogo_scaffold("my_project")

# 4) Place inputs into my_project/input/

# 5) Run
echogo_run("my_project/input", "my_project/results")

# 6) Inspect consensus tables & plots
```

---

# 📑 Citation

If you use EchoGO in publications, please cite:

> Escobar-Sierra C., Langschied F., Inostroza P.A. (2025).  
> **EchoGO: Cross-Species Consensus Functional Enrichment Analysis.**  
> Zenodo. DOI: "https://doi.org/10.5281/zenodo.17658715"

Full citation entry is included in `inst/CITATION`.

---

# 🐛 Issues & Contact

Please report bugs, suggestions, or feature requests at:  
**https://github.com/miloes114/EchoGo/issues**

For general questions, contact the maintainers.

---

# 💬 Acknowledgements

EchoGO was developed by **Camilo Escobar-Sierra**, **Felix Langschied**, and **Pedro A. Inostroza**, with additional input from collaborators and the community.

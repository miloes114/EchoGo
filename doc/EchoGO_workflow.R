## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse   = TRUE,
  comment    = "#>",
  fig.width  = 7,
  fig.height = 4.5,
  error      = TRUE
)

RUN_HEAVY <- identical(Sys.getenv("NOT_CRAN"), "true")

suppressPackageStartupMessages({
  library(EchoGO)
  library(readxl)
  library(dplyr)
})

## ----quick_start--------------------------------------------------------------
# See the welcome message with defaults and helper functions
echogo_help()

## ----eval=FALSE---------------------------------------------------------------
# suppressPackageStartupMessages(library(EchoGO))

## ----demo_paths---------------------------------------------------------------
demo_in  <- echogo_demo_path()          # Location of demo inputs
demo_out <- echogo_demo_results_path()  # Location of frozen demo results

cat("Demo input path:\n", demo_in, "\n\n")
cat("Demo output path:\n", demo_out, "\n")

## ----eval=FALSE---------------------------------------------------------------
# echogo_open_demo()  # Opens demo input/results and lists key files

## ----show_defaults------------------------------------------------------------
# Check current defaults
default_orgdb <- getOption("EchoGO.default_orgdb", "org.Mm.eg.db")
default_species <- getOption(
  "EchoGO.default_species",
  c("hsapiens", "mmusculus", "drerio",
    "btaurus", "rnorvegicus",
    "dmelanogaster", "ggallus")
)

cat("Default OrgDb:", default_orgdb, "\n")
cat("Default species:", paste(default_species, collapse = ", "), "\n")

## ----eval=FALSE---------------------------------------------------------------
# # Use human genome + fewer species for a quick run
# options(EchoGO.default_orgdb   = "org.Hs.eg.db")
# options(EchoGO.default_species = c("hsapiens", "mmusculus", "drerio"))

## ----list_orgdb---------------------------------------------------------------
# Check which OrgDb packages are installed
echogo_list_orgdb()

## ----eval=FALSE---------------------------------------------------------------
# # Get installation instructions (BiocManager::install commands)
# echogo_install_orgdb_instructions()
# 
# # Or let EchoGO install it automatically
# echogo_install_orgdb(orgdb = "org.Hs.eg.db")
# 
# # Ensure OrgDb is available before running
# echogo_require_orgdb(orgdb = "org.Mm.eg.db")

## ----eval=FALSE---------------------------------------------------------------
# echogo_pick_species()

## ----eval=FALSE---------------------------------------------------------------
# my_species <- c("hsapiens", "mmusculus", "drerio")
# 
# # Check if all species are recognized by g:Profiler
# echogo_preflight_species(my_species)

## ----quickstart_demo, eval=RUN_HEAVY------------------------------------------
# # Run the built-in demo (small subset, 3 species + default OrgDb)
# echogo_quickstart(run_demo = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# # Create a new project directory with input/output structure
# echogo_scaffold("my_project")

## ----eval=FALSE---------------------------------------------------------------
# echogo_run(
#   input_dir = "my_project/input",
#   outdir    = "my_project/results"
# )

## ----eval=FALSE---------------------------------------------------------------
# run_full_echogo(
#   input_dir   = "my_project/input",
#   species     = c("hsapiens", "mmusculus", "drerio"),
#   orgdb       = "org.Hs.eg.db",
#   outdir      = "my_project/results",
#   make_report = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# input_dir <- "E:/.../dge_RRf_BBf/input"
# outdir    <- "E:/.../dge_RRf_BBf/results"
# 
# fish_species <- c(
#   "drerio", "strutta", "gaculeatus", "olatipes", "trubripes", "amexicanus",
#   "oniloticus", "ssalar", "omykiss", "okisutch", "otshawytscha"
# )
# 
# res <- run_full_echogo(
#   input_dir              = input_dir,
#   species                = fish_species,
#   orgdb                  = "org.Dr.eg.db",
#   outdir                 = outdir,
#   strict_only            = FALSE,
#   run_evaluation         = TRUE,
#   use_trinotate_universe = TRUE,
#   make_report            = TRUE,
#   verbose                = TRUE
# )

## ----explore_outputs----------------------------------------------------------
# See what's in the frozen demo results
list.files(demo_out, recursive = TRUE, max.depth = 2)

## ----peek_consensus-----------------------------------------------------------
consensus_file <- file.path(
  demo_out, "consensus",
  "consensus_enrichment_results_with_and_without_bg.xlsx"
)

stopifnot(file.exists(consensus_file))

cons <- read_xlsx(consensus_file)

# View key columns
head(cons[, c(
  "term_id", "term_name", "ontology",
  "consensus_score", "consensus_score_all",
  "source_origin", "significant_in_any"
)])

## ----peek_goseq---------------------------------------------------------------
goseq_file <- file.path(demo_out, "goseq", "GOseq_enrichment_full_annotated.csv")
if (file.exists(goseq_file)) {
  goseq <- read.csv(goseq_file)
  head(goseq[, intersect(c("category","term","ontology","numDEInCat","numInCat","foldEnrichment","over_represented_FDR"), names(goseq))])
}

## ----peek_gprofiler-----------------------------------------------------------
# Results are stored per species and background mode
gprofiler_dir <- file.path(demo_out, "gprofiler")
list.files(gprofiler_dir, recursive = TRUE)

## ----example_load-------------------------------------------------------------
consensus_file <- file.path(
  demo_out, "consensus",
  "consensus_enrichment_results_with_and_without_bg.xlsx"
)
cons <- read_xlsx(consensus_file)

top_terms <- cons %>%
  dplyr::filter(significant_in_any) %>%
  dplyr::arrange(dplyr::desc(consensus_score_all)) %>%
  dplyr::slice(1:20)

# Show key columns
top_terms %>%
  dplyr::select(
    term_name, ontology,
    consensus_score, consensus_score_all,
    source_origin
  ) %>%
  head(10)

## ----example_plot, fig.width=8, fig.height=5----------------------------------
library(ggplot2)
library(stringr)

plot_df <- cons %>%
  dplyr::filter(significant_in_any, !is.na(consensus_score_all)) %>%
  dplyr::slice_max(consensus_score_all, n = 15, with_ties = FALSE) %>%
  dplyr::arrange(consensus_score_all) %>%
  dplyr::mutate(
    term_label = stringr::str_trunc(
      paste0(term_name, " (", term_id, ")"),
      45
    ),
    term_label = factor(term_label, levels = term_label)
  )

ggplot(plot_df, aes(
  y     = term_label,
  x     = consensus_score_all,
  color = ontology
)) +
  geom_segment(aes(
    yend = term_label,
    xend = 0
  ),
  linewidth = 0.6,
  color     = "grey70"
  ) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Dark2", name = "Ontology") +
  labs(
    x     = "Exploratory consensus score",
    y     = NULL,
    title = "Top 15 EchoGO terms by exploratory consensus"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y     = element_text(size = 9),
    legend.position = "bottom",
    plot.title      = element_text(hjust = 0.5, size = 12)
  )

## ----eval=FALSE---------------------------------------------------------------
# run_full_echogo(
#   input_dir = NULL,
#   goseq_file = NULL,
#   trinotate_file = NULL,
#   de_file = NULL,
#   count_matrix_file = NULL,
#   species = getOption("EchoGO.default_species", c("hsapiens","mmusculus","drerio")),
#   species_expr = NULL,
#   orgdb   = getOption("EchoGO.default_orgdb", "org.Dr.eg.db"),
#   outdir  = "echogo_out",
#   make_report = TRUE,
#   report_sections = c("overview","goseq","gprofiler","consensus","rrvgo","networks"),
#   report_top_n = 25,
#   report_theme = "flatly",
#   report_template = NULL,
#   strict_only = FALSE,
#   run_evaluation = TRUE,
#   use_trinotate_universe = FALSE,
#   verbose = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# echogo_preflight_species(
#   species = NULL,               # Character vector of organism codes
#   suggest = TRUE                # Suggest corrections for invalid codes?
# )

## ----eval=FALSE---------------------------------------------------------------
# echogo_install_orgdb("org.Mm.eg.db")
# echogo_require_orgdb("org.Mm.eg.db")

## ----eval=FALSE---------------------------------------------------------------
# echogo_preflight_species("hsapiens")  # Use proper code, not common name

## ----eval=FALSE---------------------------------------------------------------
# options(EchoGO.gprofiler_timeout = 120)

## ----session-info, echo=FALSE-------------------------------------------------
sessionInfo()


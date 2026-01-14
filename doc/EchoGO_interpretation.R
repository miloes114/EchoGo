## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  error    = TRUE
)
options(rmarkdown.html_vignette.check_title = FALSE)

suppressPackageStartupMessages({
  library(EchoGO)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})

# Frozen demo results shipped with the package
out <- echogo_demo_results_path()

## ----load_consensus-----------------------------------------------------------
cons_xlsx <- file.path(
  out, "consensus", 
  "consensus_enrichment_results_with_and_without_bg.xlsx"
)
stopifnot(file.exists(cons_xlsx))

cons <- read_xlsx(cons_xlsx)
dplyr::glimpse(cons)

## ----prev_pseudocode, eval=FALSE----------------------------------------------
# bg_prev <- num_species_gprof_bg / n_bg_cols
# nobg_prev <- num_species_gprof_nobg / n_nobg_cols

## ----logp_pseudocode, eval=FALSE----------------------------------------------
# cap_logp <- function(p, P_CAP) {
#   p <- dplyr::coalesce(p, 1)
#   p <- pmax(p, .Machine$double.xmin)
#   pmin(-log10(p), P_CAP) / P_CAP  # scaled to [0, 1]
# }

## ----pval_components, eval=FALSE----------------------------------------------
# # Strict p-value component: only g:Profiler with background
# comp_p_strict <- cap_logp(min_pval_gprof_bg, P_CAP)
# 
# # Exploratory p-value component: best signal across all methods
# comp_p_all <- pmax(
#   cap_logp(min_pval_gprof_bg, P_CAP),
#   cap_logp(min_pval_gprof_nobg, P_CAP),
#   cap_logp(min_pval_goseq, P_CAP)
# )

## ----fold_pseudocode, eval=FALSE----------------------------------------------
# log2p1_cap <- function(x, cap) {
#   log2(1 + pmin(dplyr::coalesce(x, 0), cap))
# }
# 
# depth_w <- function(d, maxd = 12L) {
#   d <- suppressWarnings(as.integer(d))
#   ifelse(is.na(d), 1, pmin(d / maxd, 1))
# }

## ----fold_components, eval=FALSE----------------------------------------------
# comp_fold_bg <- depth_w(depth) * log2p1_cap(avg_fold_gprof_bg, FOLD_CAP)
# comp_fold_nb <- depth_w(depth) * log2p1_cap(avg_fold_gprof_nobg, FOLD_CAP)
# comp_fold_gs <- depth_w(depth) * log2p1_cap(fold_enrichment_goseq, FOLD_CAP)

## ----goseq_pseudocode, eval=FALSE---------------------------------------------
# comp_goseq <- as.integer(isTRUE(in_goseq))

## ----weights_pseudocode, eval=FALSE-------------------------------------------
# W_STRICT <- list(
#   w_goseq = 1.0,
#   w_bgprev = 0.8,
#   w_foldbg = 0.40,
#   w_foldgs = 0.40,
#   w_p = 0.40
# )
# 
# W_ALL <- list(
#   w_goseq = 1.0,
#   w_bgprev = 0.7,
#   w_nbprev = 0.3,
#   w_foldbg = 0.35,
#   w_foldnb = 0.25,
#   w_foldgs = 0.35,
#   w_p = 0.35
# )

## ----strict_score, eval=FALSE-------------------------------------------------
# consensus_score <- W_STRICT$w_goseq * comp_goseq +
#   W_STRICT$w_bgprev * bg_prev +
#   W_STRICT$w_foldbg * comp_fold_bg +
#   W_STRICT$w_foldgs * comp_fold_gs +
#   W_STRICT$w_p * comp_p_strict

## ----all_score, eval=FALSE----------------------------------------------------
# consensus_score_all <- W_ALL$w_goseq * comp_goseq +
#   W_ALL$w_bgprev * bg_prev +
#   W_ALL$w_nbprev * nobg_prev +
#   W_ALL$w_foldbg * comp_fold_bg +
#   W_ALL$w_foldnb * comp_fold_nb +
#   W_ALL$w_foldgs * comp_fold_gs +
#   W_ALL$w_p * comp_p_all

## ----read_consensus-----------------------------------------------------------
top_consensus <- cons %>%
  dplyr::filter(significant_in_any) %>%
  dplyr::arrange(dplyr::desc(consensus_score)) %>%
  dplyr::select(
    term_id, term_name, ontology,
    num_species_gprof_bg, num_species_gprof_nobg,
    consensus_score, consensus_score_all,
    origin, source_origin
  ) %>%
  head(25)

top_consensus

## ----rrvgo_check--------------------------------------------------------------
rr_true <- file.path(out, "rrvgo", "rrvgo_true_consensus_with_bg")
rr_all <- file.path(out, "rrvgo", "rrvgo_exploratory_all_significant")

data.frame(
  mode = c("True_Consensus_with_BG", "Exploratory_all"),
  exists = c(dir.exists(rr_true), dir.exists(rr_all))
)

## ----rrvgo_explore------------------------------------------------------------
if (dir.exists(rr_true)) {
  head(list.files(rr_true, recursive = TRUE), 10)
}

## ----network_explore----------------------------------------------------------
net_dir <- file.path(out, "networks")
head(list.files(net_dir, recursive = TRUE), 10)

## ----minimal_plot, fig.width=6, fig.height=5, dpi=120-------------------------

# Improved minimal plotting example (exploratory)

plot_df <- cons %>%
  dplyr::filter(significant_in_any, !is.na(consensus_score_all)) %>%
  dplyr::slice_max(consensus_score_all, n = 15, with_ties = FALSE) %>%   # smaller = cleaner
  dplyr::mutate(
    term = stringr::str_trunc(
      paste0(term_name, " (", term_id, ")"), 
      50
    )
  ) %>%
  dplyr::arrange(consensus_score_all) %>%
  dplyr::mutate(term = factor(term, levels = term))

ggplot(plot_df, aes(
  y = term,
  x = consensus_score_all
)) +
  geom_segment(aes(
    yend = term,
    xend = 0
  ), linewidth = 0.6, color = "grey70") +
  geom_point(
    aes(color = ontology),
    size = 2.5
  ) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = "Consensus score (exploratory)",
    y = NULL,
    title = "Top 15 exploratory GO/KEGG terms by consensus score"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )


## ----session-info, echo=FALSE-------------------------------------------------
sessionInfo()


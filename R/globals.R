# Imports and non-standard-evaluation bindings used across EchoGO modules.
#'
#' @importFrom grDevices dev.off pdf
#' @importFrom stats alias ave runif setNames
#' @importFrom utils capture.output tail
NULL

utils::globalVariables(c(
  ".", ".bg_sum", ".depth_calc", ".nobg_sum", ".p_kegg", ".src_calc",
  ":=", "EQI", "alias", "all_genes", "avg_degree_bg",
  "avg_degree_exploratory", "avg_fold_gprof_bg", "avg_fold_gprof_nobg",
  "alternative_context_support_fraction", "alternative_context_support_n",
  "alternative_queried_context_n",
  "bg_prev", "category", "clean_go_term", "cluster", "community",
  "community_label", "comp_fold_bg", "comp_fold_gs", "comp_fold_nb",
  "comp_goseq", "comp_p_all", "comp_p_strict", "consensus_score",
  "delta_avg_degree", "delta_edges", "delta_terms", "depth", "enrichment", "evidence_profile",
  "foldEnrichment", "fold_enrichment_goseq", "gene_count", "gene_names",
  "genes_flat", "go", "hiddenLabel", "id", "in_goseq", "in_gsq",
  "method", "min_pval_goseq", "min_pval_gprof_bg", "min_pval_gprof_nobg",
  "gprof_custom_context_support_n", "name", "ncbi", "nobg_prev", "numDEInCat", "numInCat",
  "num_species_gprof_bg", "num_species_gprof_nobg", "ontology",
  "ontology_gprof", "organism", "origin", "over_represented_FDR",
  "primary_name", "score", "sig_any", "significant_in_any", "size",
  "queried_context_n", "size_attr", "sources_count", "target_goseq_supported", "term", "term_display", "term_id",
  "term_name", "top_term", "total_edges_bg", "total_edges_exploratory",
  "total_terms_bg", "total_terms_exploratory", "transcript_id", "x", "y"
))

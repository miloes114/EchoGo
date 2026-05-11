.echogo_truthy <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x)) return(!is.na(x) & x != 0)
  if (is.character(x)) return(tolower(trimws(x)) %in% c("true", "1", "yes", "y"))
  rep(FALSE, length(x))
}

.echogo_norm_go_ontology <- function(x) {
  dplyr::case_when(
    x %in% c("GO:BP", "BP") ~ "BP",
    x %in% c("GO:MF", "MF") ~ "MF",
    x %in% c("GO:CC", "CC") ~ "CC",
    TRUE ~ as.character(x)
  )
}

.echogo_rrvgo_mode_tables <- function(consensus_df) {
  if (!is.data.frame(consensus_df) || !nrow(consensus_df)) {
    empty <- consensus_df[0, , drop = FALSE]
    return(list(
      true_consensus_with_bg = empty,
      conservative_bg_supported = empty,
      exploratory_all_significant = empty
    ))
  }

  df <- consensus_df %>%
    dplyr::mutate(
      ontology = .echogo_norm_go_ontology(.data$ontology),
      significant_in_any = .echogo_truthy(.data$significant_in_any),
      origin = as.character(.data$origin)
    ) %>%
    dplyr::filter(.data$ontology %in% c("BP", "MF", "CC"), .data$significant_in_any)

  list(
    true_consensus_with_bg = df %>%
      dplyr::filter(.data$origin == "GO terms - Consensus (with BG)"),
    conservative_bg_supported = df %>%
      dplyr::filter(.data$origin %in% c(
        "GO terms - GOseq only",
        "GO terms - g:Profiler only (with BG)",
        "GO terms - Consensus (with BG)"
      )),
    exploratory_all_significant = df
  )
}

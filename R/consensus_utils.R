#' @keywords internal
.echogo_add_consensus_scores <- function(df) {
  if (is.null(df) || !nrow(df)) return(df)
  # Delegate to the canonical scorer so formulas & column names match the Rmd/Guide
  score_consensus_terms(df)
}

#' @keywords internal
.echogo_add_consensus_scores <- function(df) {
  if (is.null(df) || !nrow(df)) return(df)
  # Historical wide-table compatibility only; canonical v0.1.4 execution does
  # not call this adapter and it deliberately does not recreate a score.
  score_consensus_terms(df)
}

# Report-safe aliases for hypothesis drill-down -------------------------------
# Preserve the existing returned components and add short aliases consumed by
# the biological report. This is presentation compatibility only.

.echogo_hypothesis_drilldown_base <- echogo_hypothesis_drilldown

echogo_hypothesis_drilldown <- function(evidence_table, source_provenance, term_id,
                                        annotation_provenance = NULL) {
  out <- .echogo_hypothesis_drilldown_base(
    evidence_table = evidence_table,
    source_provenance = source_provenance,
    term_id = term_id,
    annotation_provenance = annotation_provenance
  )
  out$term <- out$term_evidence
  out$sources <- out$source_provenance
  out$annotation <- out$annotation_provenance
  out
}

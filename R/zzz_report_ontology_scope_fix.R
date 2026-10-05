# Report ontology-scope compatibility fix ------------------------------------
#
# The v0.1.4 report layer uses function arguments named `ontology` alongside
# evidence tables that also contain a column named `ontology`. Inside dplyr's
# data mask, an unqualified bare `ontology` can therefore resolve to the data
# column rather than the function argument. In the affected helpers this made
# expressions such as
#
#   .data$ontology %in% c(ontology, paste0("GO:", ontology))
#
# effectively self-referential, allowing BP/MF/CC rows to leak into one
# another. The canonical evidence table itself is unchanged; this file fixes
# only report/network subset selection.
#
# This file is intentionally collated late (`zzz_`) so it can wrap the three
# already-established report helpers without duplicating their scientific
# logic. A future refactor may move the same explicit `.env`/base-subset logic
# into the source helpers directly, but the behaviour below is canonical for
# the v0.1.4 release.

.echogo_normalize_ontology_code <- function(x) {
  if (is.null(x) || !length(x)) return(NA_character_)
  raw <- toupper(trimws(as.character(x)))
  raw <- sub("^GO[[:space:]_:-]*", "", raw)
  key <- gsub("[^A-Z]", "", raw)
  out <- rep(NA_character_, length(key))
  out[key %in% c("BP", "BIOLOGICALPROCESS")] <- "BP"
  out[key %in% c("MF", "MOLECULARFUNCTION")] <- "MF"
  out[key %in% c("CC", "CELLULARCOMPONENT")] <- "CC"
  out
}

.echogo_subset_evidence_ontology <- function(evidence, ontology) {
  if (is.null(evidence) || !is.data.frame(evidence) || !nrow(evidence)) {
    return(evidence)
  }
  if (!"ontology" %in% names(evidence)) {
    stop("Ontology-scoped report products require an `ontology` column.", call. = FALSE)
  }

  ont <- .echogo_normalize_ontology_code(ontology[[1]])
  if (is.na(ont)) {
    stop(
      "Unsupported GO ontology code: ", as.character(ontology[[1]]),
      ". Expected BP, MF, or CC.",
      call. = FALSE
    )
  }

  evidence_codes <- .echogo_normalize_ontology_code(evidence$ontology)
  evidence[!is.na(evidence_codes) & evidence_codes == ont, , drop = FALSE]
}

# Preserve the established implementations once, then pre-scope their input
# explicitly. This avoids relying on tidy-evaluation name resolution.
if (!exists(".echogo_report_network_data_unscoped", inherits = FALSE)) {
  .echogo_report_network_data_unscoped <- .echogo_report_network_data
}

.echogo_report_network_data <- function(evidence, profiles, ontology,
                                        min_shared_genes = 2L,
                                        min_gene_count = 2L,
                                        max_terms = 120L) {
  ont <- .echogo_normalize_ontology_code(ontology[[1]])
  scoped <- .echogo_subset_evidence_ontology(evidence, ont)
  .echogo_report_network_data_unscoped(
    evidence = scoped,
    profiles = profiles,
    ontology = ont,
    min_shared_genes = min_shared_genes,
    min_gene_count = min_gene_count,
    max_terms = max_terms
  )
}

if (!exists(".echogo_key_findings_data_unscoped", inherits = FALSE)) {
  .echogo_key_findings_data_unscoped <- .echogo_key_findings_data
}

.echogo_key_findings_data <- function(evidence,
                                      ontology = NULL,
                                      max_target_recovered = 8L,
                                      max_target_only = 4L,
                                      max_hypotheses = 8L) {
  scoped <- evidence
  if (!is.null(ontology)) {
    scoped <- .echogo_subset_evidence_ontology(evidence, ontology)
  }

  # The underlying helper no longer needs to evaluate its own ontology filter;
  # the input has already been scoped with ordinary vector comparison.
  .echogo_key_findings_data_unscoped(
    evidence = scoped,
    ontology = NULL,
    max_target_recovered = max_target_recovered,
    max_target_only = max_target_only,
    max_hypotheses = max_hypotheses
  )
}

# This helper is short enough to define directly with explicit ontology
# scoping, avoiding the ambiguous bare-name filter completely.
.echogo_landscape_terms <- function(evidence, ontology, max_terms_per_profile = 12L) {
  profile_order <- c("TARGET_ONLY", "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT")
  scoped <- .echogo_subset_evidence_ontology(evidence, ontology)
  if (is.null(scoped) || !is.data.frame(scoped) || !nrow(scoped)) {
    return(tibble::tibble())
  }

  scoped |>
    dplyr::filter(.data$evidence_profile %in% profile_order) |>
    dplyr::mutate(profile_order = match(.data$evidence_profile, profile_order)) |>
    dplyr::arrange(.data$profile_order, .data$display_order, .data$term_id) |>
    dplyr::group_by(.data$evidence_profile) |>
    dplyr::slice_head(n = max_terms_per_profile) |>
    dplyr::ungroup()
}

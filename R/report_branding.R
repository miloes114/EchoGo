# EchoGO report branding and biological summary helpers ----------------------
#
# These helpers affect report presentation only. They must never change the
# exact-term evidence table, source-local statistics, recurrence, or semantic
# products.

.echogo_interpretation_counts <- function(evidence, contexts = NULL) {
  if (is.null(evidence) || !is.data.frame(evidence) || !nrow(evidence)) {
    return(tibble::tibble(
      target_supported_n = 0L,
      target_plus_context_n = 0L,
      target_only_n = 0L,
      hypothesis_n = 0L,
      alternative_reference_n = if (is.data.frame(contexts)) {
        sum(contexts$context_role %in% "ALTERNATIVE", na.rm = TRUE)
      } else 0L
    ))
  }

  profile <- as.character(evidence$evidence_profile)
  tibble::tibble(
    target_supported_n = sum(profile %in% c("TARGET_ONLY", "TARGET_PLUS_CONTEXT"), na.rm = TRUE),
    target_plus_context_n = sum(profile %in% "TARGET_PLUS_CONTEXT", na.rm = TRUE),
    target_only_n = sum(profile %in% "TARGET_ONLY", na.rm = TRUE),
    hypothesis_n = sum(profile %in% "ALTERNATIVE_CONTEXT", na.rm = TRUE),
    alternative_reference_n = if (is.data.frame(contexts) && "context_role" %in% names(contexts)) {
      sum(contexts$context_role %in% "ALTERNATIVE", na.rm = TRUE)
    } else if ("alternative_queried_context_n" %in% names(evidence)) {
      x <- suppressWarnings(as.integer(evidence$alternative_queried_context_n))
      x <- x[is.finite(x)]
      if (length(x)) max(x) else 0L
    } else {
      0L
    }
  )
}

.echogo_primary_metric_cards <- function(evidence, contexts = NULL) {
  x <- .echogo_interpretation_counts(evidence, contexts)
  tibble::tibble(
    value = c(
      x$target_supported_n,
      x$target_plus_context_n,
      x$target_only_n,
      x$hypothesis_n,
      x$alternative_reference_n
    ),
    label = c(
      "GO terms supported by Target GOseq",
      "Target + context",
      "Target only",
      "Context-derived hypotheses",
      "Selected annotation contexts"
    )
  )
}

.echogo_representability_statement <- function(mapping) {
  s <- .echogo_mapping_summary(mapping)
  if (!nrow(s)) return(NA_character_)
  tested <- as.integer(s$tested_entities[[1]])
  submitted <- as.integer(s$submitted_background_names[[1]])
  if (!is.finite(tested) || tested <= 0L) return(NA_character_)
  pct <- 100 * submitted / tested
  paste0(
    format(submitted, big.mark = ","), " of ", format(tested, big.mark = ","),
    " tested target entities (", sprintf("%.1f", pct),
    "%) were represented by unique names in the matched-background comparative analysis."
  )
}

.echogo_logo_source_path <- function(explicit = NULL) {
  if (!is.null(explicit) && length(explicit) && !is.na(explicit[[1]]) && nzchar(explicit[[1]])) {
    explicit_path <- as.character(explicit[[1]])
    if (!file.exists(explicit_path)) return(NA_character_)
    return(normalizePath(explicit_path, winslash = "/", mustWork = TRUE))
  }

  candidates <- character()

  installed <- tryCatch(
    system.file("reports", "assets", "echogo-logo.png", package = "EchoGO"),
    error = function(e) ""
  )
  if (nzchar(installed)) candidates <- c(candidates, installed)

  candidates <- c(
    candidates,
    file.path("inst", "reports", "assets", "echogo-logo.png"),
    file.path(getwd(), "inst", "reports", "assets", "echogo-logo.png")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) normalizePath(hit[[1]], winslash = "/", mustWork = TRUE) else NA_character_
}

.echogo_stage_report_logo <- function(report_dir, logo_path = NULL,
                                      staged_name = "echogo-logo.png") {
  source <- .echogo_logo_source_path(logo_path)
  if (is.na(source) || !file.exists(source)) return(NA_character_)
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(report_dir, staged_name)
  ok <- file.copy(source, dest, overwrite = TRUE, copy.date = TRUE)
  if (!isTRUE(ok) || !file.exists(dest)) return(NA_character_)
  normalizePath(dest, winslash = "/", mustWork = TRUE)
}

.echogo_html_escape_attr <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

.echogo_report_logo_html <- function(src = "echogo-logo.png",
                                     alt = "EchoGO logo") {
  if (is.null(src) || length(src) == 0L || is.na(src[[1]]) || !nzchar(src[[1]])) return("")
  paste0(
    "<div class='hero-logo-shell'>",
    "<img class='hero-echogo-logo' src='", .echogo_html_escape_attr(src[[1]]),
    "' alt='", .echogo_html_escape_attr(alt), "'/>",
    "</div>"
  )
}

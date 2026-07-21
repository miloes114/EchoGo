.echogo_detect_delim <- function(file,
                                 candidates = c("\t", ";", ","),
                                 n = 50L,
                                 expected = NULL) {
  header <- tryCatch(
    readLines(file, n = 1L, warn = FALSE),
    error = function(e) ""
  )
  header <- enc2utf8(header)

  score_header <- function(sep) {
    parts <- strsplit(header, sep, fixed = TRUE)[[1]]
    parts <- gsub("^\ufeff", "", parts)
    parts <- gsub("\u00a0", " ", parts, fixed = TRUE)
    parts <- tolower(trimws(parts))

    if (length(parts) < 2L) return(-Inf)
    if (is.null(expected)) return(0)
    sum(parts %in% tolower(expected))
  }

  if (!is.null(expected)) {
    scores <- vapply(candidates, score_header, numeric(1))
    best <- max(scores, na.rm = TRUE)
    if (is.finite(best) && best > 0) {
      return(candidates[which.max(scores)])
    }
  }

  best <- candidates[1]
  best_score <- -Inf

  for (sep in candidates) {
    fields <- tryCatch(
      utils::count.fields(
        file,
        sep = sep,
        quote = "\"",
        comment.char = "",
        skipNul = TRUE
      ),
      error = function(e) integer(0)
    )
    if (!length(fields)) next

    fields <- fields[seq_len(min(length(fields), n))]
    median_fields <- stats::median(fields)
    if (is.na(median_fields) || median_fields < 2) next

    penalty <- if (median_fields > 50) 50 else 0
    score <- (median_fields - stats::sd(fields)) - penalty
    if (is.finite(score) && score > best_score) {
      best_score <- score
      best <- sep
    }
  }

  best
}

.echogo_read_delim_robust <- function(file,
                                      candidates = c("\t", ";", ","),
                                      dec = ".",
                                      expected = NULL,
                                      ...) {
  sep <- .echogo_detect_delim(
    file,
    candidates = candidates,
    expected = expected
  )

  if (identical(sep, ";")) {
    first <- readLines(file, n = 2L, warn = FALSE)
    if (any(grepl("\\d,\\d", first)) && !any(grepl("\\d\\.\\d", first))) {
      dec <- ","
    }
  }

  out <- tryCatch(
    utils::read.delim(
      file,
      sep = sep,
      dec = dec,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "\"",
      fill = TRUE,
      comment.char = "",
      ...
    ),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)

  if (requireNamespace("data.table", quietly = TRUE)) {
    out <- tryCatch(
      data.table::fread(
        file,
        sep = sep,
        dec = dec,
        data.table = FALSE,
        fill = TRUE,
        quote = "\"",
        showProgress = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(out)) return(out)
  }

  stop("Failed to read file robustly: ", file, call. = FALSE)
}

.echogo_check_goseq_parse <- function(df) {
  spill_columns <- grep("^V[0-9]+$", names(df), value = TRUE)

  if ("gene_ids" %in% names(df) && length(spill_columns) >= 3L) {
    stop(
      "The GOseq table appears to have been parsed incorrectly. ",
      "The gene_ids field may have been split into extra columns: ",
      paste(utils::head(spill_columns, 10), collapse = ", "),
      ". Use a tab-delimited file or a correctly quoted CSV.",
      call. = FALSE
    )
  }

  invisible(df)
}

.echogo_find_one <- function(patterns,
                             must = TRUE,
                             label = "input",
                             verbose = FALSE) {
  for (pattern in patterns) {
    hits <- sort(Sys.glob(pattern))

    if (length(hits) == 1L) {
      if (isTRUE(verbose)) {
        message("Resolved ", label, ": ", basename(hits))
      }
      return(hits)
    }

    if (length(hits) > 1L) {
      stop(
        "Multiple candidate ", label, " files matched:\n  ",
        paste(basename(hits), collapse = "\n  "),
        "\nProvide the file explicitly.",
        call. = FALSE
      )
    }
  }

  if (isTRUE(must)) {
    stop("Could not resolve ", label, ".", call. = FALSE)
  }

  NULL
}

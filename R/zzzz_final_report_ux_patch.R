# Report UX and asset-integrity helpers -------------------------------------
#
# This late-loaded compatibility layer implements final report and asset
# refinements without changing the frozen EchoGO evidence model.

# PNG integrity -------------------------------------------------------------

.echogo_png_is_renderable <- function(path) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(FALSE)
  }
  info <- file.info(path)
  if (is.na(info$size) || info$size < 16L) return(FALSE)
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  sig <- readBin(con, what = "raw", n = 8L)
  identical(as.integer(sig), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L))
}

.echogo_rrvgo_treemap_png_from_clusters <- function(cluster_file, png_file,
                                                     width = 1800L,
                                                     height = 1500L,
                                                     res = 160L) {
  if (!file.exists(cluster_file)) return(FALSE)
  if (!requireNamespace("treemap", quietly = TRUE)) return(FALSE)
  if (!requireNamespace("scales", quietly = TRUE)) return(FALSE)

  reduced <- tryCatch(utils::read.csv(cluster_file, stringsAsFactors = FALSE), error = function(e) NULL)
  required <- c("parentTerm", "term", "size", "score")
  if (is.null(reduced) || !nrow(reduced) || !all(required %in% names(reduced))) return(FALSE)

  top <- reduced |>
    dplyr::arrange(dplyr::desc(.data$score)) |>
    utils::head(300L) |>
    dplyr::filter(!is.na(.data$size), is.finite(.data$size), .data$size > 0) |>
    dplyr::mutate(
      parentTerm = dplyr::if_else(is.na(.data$parentTerm) | !nzchar(.data$parentTerm), "Other", as.character(.data$parentTerm)),
      term = as.character(.data$term),
      scaled_size = log1p(.data$size)
    )
  if (!nrow(top)) return(FALSE)

  dir.create(dirname(png_file), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(png_file, ".tmp.png")
  if (file.exists(tmp)) unlink(tmp, force = TRUE)

  ok <- tryCatch({
    grDevices::png(tmp, width = width, height = height, res = res, bg = "white")
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    treemap::treemap(
      top,
      index = c("parentTerm", "term"),
      vSize = "scaled_size",
      type = "index",
      palette = scales::hue_pal()(max(1L, length(unique(top$parentTerm)))),
      fontcolor.labels = c("#FFFFFFDD", "#00000080"),
      bg.labels = 0,
      border.col = "#00000080"
    )
    grDevices::dev.off()
    TRUE
  }, error = function(e) FALSE)

  if (!isTRUE(ok) || !.echogo_png_is_renderable(tmp)) {
    if (file.exists(tmp)) unlink(tmp, force = TRUE)
    return(FALSE)
  }
  if (file.exists(png_file)) unlink(png_file, force = TRUE)
  file.rename(tmp, png_file)
  .echogo_png_is_renderable(png_file)
}

.echogo_ensure_rrvgo_report_pngs <- function(rrvgo_dir) {
  if (!dir.exists(rrvgo_dir)) return(invisible(tibble::tibble()))
  clusters <- list.files(
    rrvgo_dir,
    pattern = "rrvgo_(BP|MF|CC)_clusters\\.csv$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (!length(clusters)) return(invisible(tibble::tibble()))

  rows <- lapply(clusters, function(cluster_file) {
    stem <- sub("_clusters\\.csv$", "", cluster_file, ignore.case = TRUE)
    treemap_png <- paste0(stem, "_treemap.png")
    before <- .echogo_png_is_renderable(treemap_png)
    repaired <- FALSE
    if (!before) {
      repaired <- isTRUE(.echogo_rrvgo_treemap_png_from_clusters(cluster_file, treemap_png))
    }
    tibble::tibble(
      cluster_file = normalizePath(cluster_file, winslash = "/", mustWork = FALSE),
      png_file = normalizePath(treemap_png, winslash = "/", mustWork = FALSE),
      valid_before = before,
      repaired = repaired,
      valid_after = .echogo_png_is_renderable(treemap_png)
    )
  })
  dplyr::bind_rows(rows)
}

# Key Findings ordering ------------------------------------------------------
# Target-supported findings are selected from the target experiment first:
# GOseq adjusted p-value -> reference recurrence -> deterministic order.
# New hypotheses remain recurrence-first.  No cross-reference p-values enter
# either stream.

.echogo_key_findings_data <- function(evidence,
                                      ontology = NULL,
                                      max_target_recovered = 8L,
                                      max_target_only = 4L,
                                      max_hypotheses = 8L) {
  if (is.null(evidence) || !is.data.frame(evidence) || !nrow(evidence)) {
    return(tibble::tibble())
  }

  required <- c(
    "term_id", "term_name", "ontology", "evidence_profile",
    "target_goseq_adjusted_p", "alternative_context_support_n",
    "alternative_queried_context_n", "alternative_context_support_fraction",
    "contributing_genes"
  )
  missing <- setdiff(required, names(evidence))
  if (length(missing)) {
    stop("Key-findings display requires columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!"display_order" %in% names(evidence)) evidence <- .echogo_add_display_order(evidence)

  allowed_profiles <- c("TARGET_PLUS_CONTEXT", "TARGET_ONLY", "ALTERNATIVE_CONTEXT")
  profile_order <- c("TARGET_PLUS_CONTEXT", "TARGET_ONLY", "ALTERNATIVE_CONTEXT")
  profile_levels <- unname(.echogo_key_findings_group_labels[profile_order])

  d <- tibble::as_tibble(evidence) |>
    dplyr::filter(.data$evidence_profile %in% allowed_profiles)

  if (!is.null(ontology)) {
    ontology_code <- toupper(as.character(ontology[[1]]))
    if (!ontology_code %in% c("BP", "MF", "CC")) stop("ontology must be BP, MF, or CC", call. = FALSE)
    keep <- sub("^GO:", "", toupper(as.character(d$ontology))) == ontology_code
    d <- d[keep, , drop = FALSE]
  }
  if (!nrow(d)) return(tibble::tibble())

  d <- d |>
    dplyr::mutate(
      ontology = sub("^GO:", "", toupper(as.character(.data$ontology))),
      alternative_context_support_n = dplyr::coalesce(suppressWarnings(as.integer(.data$alternative_context_support_n)), 0L),
      alternative_queried_context_n = dplyr::coalesce(suppressWarnings(as.integer(.data$alternative_queried_context_n)), 0L),
      alternative_context_support_fraction = dplyr::if_else(
        .data$alternative_queried_context_n > 0L,
        .data$alternative_context_support_n / .data$alternative_queried_context_n,
        0
      ),
      target_goseq_adjusted_p = suppressWarnings(as.numeric(.data$target_goseq_adjusted_p)),
      contributing_gene_n = .echogo_key_findings_gene_count(.data$contributing_genes),
      finding_group = unname(.echogo_key_findings_group_labels[.data$evidence_profile]),
      profile_order = match(.data$evidence_profile, profile_order),
      support_label = dplyr::if_else(
        .data$alternative_queried_context_n > 0L,
        paste0(.data$alternative_context_support_n, "/", .data$alternative_queried_context_n),
        "no other references"
      )
    )

  plus <- d |>
    dplyr::filter(.data$evidence_profile == "TARGET_PLUS_CONTEXT") |>
    dplyr::arrange(
      .data$target_goseq_adjusted_p,
      dplyr::desc(.data$alternative_context_support_n),
      .data$display_order,
      .data$term_id
    ) |>
    dplyr::slice_head(n = as.integer(max_target_recovered)) |>
    dplyr::mutate(within_stream_order = dplyr::row_number())

  target_only <- d |>
    dplyr::filter(.data$evidence_profile == "TARGET_ONLY") |>
    dplyr::arrange(.data$target_goseq_adjusted_p, .data$display_order, .data$term_id) |>
    dplyr::slice_head(n = as.integer(max_target_only)) |>
    dplyr::mutate(within_stream_order = dplyr::row_number())

  hypothesis <- d |>
    dplyr::filter(.data$evidence_profile == "ALTERNATIVE_CONTEXT") |>
    dplyr::arrange(
      dplyr::desc(.data$alternative_context_support_n),
      .data$display_order,
      .data$term_id
    ) |>
    dplyr::slice_head(n = as.integer(max_hypotheses)) |>
    dplyr::mutate(within_stream_order = dplyr::row_number())

  dplyr::bind_rows(plus, target_only, hypothesis) |>
    dplyr::mutate(
      finding_group = factor(.data$finding_group, levels = profile_levels)
    ) |>
    dplyr::arrange(.data$profile_order, .data$within_stream_order) |>
    dplyr::mutate(key_findings_display_order = dplyr::row_number())
}

# Evidence Landscape density ------------------------------------------------
# Cap the detailed audit view at 10 rows per evidence profile/ontology even
# when an older report template passes 12.  The complete table remains linked.

.echogo_landscape_terms <- function(evidence, ontology, max_terms_per_profile = 10L) {
  ontology_code <- toupper(as.character(ontology[[1]]))
  if (!ontology_code %in% c("BP", "MF", "CC")) stop("ontology must be BP, MF, or CC", call. = FALSE)
  max_terms_per_profile <- min(10L, as.integer(max_terms_per_profile))
  profile_order <- c("TARGET_ONLY", "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT")
  keep <- sub("^GO:", "", toupper(as.character(evidence$ontology))) == ontology_code &
    evidence$evidence_profile %in% profile_order
  evidence[keep, , drop = FALSE] |>
    dplyr::mutate(profile_order = match(.data$evidence_profile, profile_order)) |>
    dplyr::arrange(.data$profile_order, .data$display_order, .data$term_id) |>
    dplyr::group_by(.data$evidence_profile) |>
    dplyr::slice_head(n = max_terms_per_profile) |>
    dplyr::ungroup()
}

# HTML / print polish --------------------------------------------------------

.echogo_final_report_html_polish <- function(html_path) {
  if (!file.exists(html_path)) return(invisible(FALSE))
  txt <- paste(readLines(html_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

  # Normalize a few known presentation-only fragments in report copy.
  txt <- gsub("whatthe", "what the", txt, fixed = TRUE)
  txt <- gsub("Key FindingsMap", "Key Findings Map", txt, fixed = TRUE)
  txt <- gsub("experimentto", "experiment to", txt, fixed = TRUE)
  txt <- gsub("role:target reference", "role: target reference", txt, fixed = TRUE)
  txt <- gsub("biologicalthemes", "biological themes", txt, fixed = TRUE)
  txt <- gsub("(GO:[0-9]{7}) \\((GO:[0-9]{7})\\)", "\\1", txt, perl = TRUE)

  # In HTML, reference-species galleries stay collapsed.  During print/PDF we
  # temporarily open every details element and restore the prior state after
  # printing so source-level g:Profiler plots are actually included.
  print_script <- paste0(
    "<script id='echogo-print-details'>",
    "(function(){var prior=[];",
    "window.addEventListener('beforeprint',function(){prior=[];document.querySelectorAll('details').forEach(function(d){prior.push([d,d.open]);d.open=true;});});",
    "window.addEventListener('afterprint',function(){prior.forEach(function(x){x[0].open=x[1];});});",
    "})();</script>"
  )
  if (!grepl("echogo-print-details", txt, fixed = TRUE)) {
    txt <- sub("</body>", paste0(print_script, "\n</body>"), txt, fixed = TRUE)
  }

  writeLines(txt, html_path, useBytes = TRUE)
  invisible(TRUE)
}

.echogo_validate_report_image_assets <- function(html_path, report_dir) {
  if (!file.exists(html_path)) stop("Report HTML does not exist: ", html_path, call. = FALSE)
  txt <- paste(readLines(html_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  refs <- unlist(regmatches(txt, gregexpr("<img[^>]+src=['\"][^'\"]+['\"]", txt, perl = TRUE)))
  if (!length(refs)) return(invisible(TRUE))
  srcs <- sub("^.*src=['\"]([^'\"]+)['\"].*$", "\\1", refs, perl = TRUE)
  srcs <- srcs[!grepl("^(data:|https?:|file:)", srcs, ignore.case = TRUE)]
  if (!length(srcs)) return(invisible(TRUE))

  failures <- character()
  for (src in unique(srcs)) {
    path <- file.path(report_dir, gsub("/", .Platform$file.sep, src, fixed = TRUE))
    if (!file.exists(path)) {
      failures <- c(failures, paste0(src, " [missing]"))
    } else if (tolower(tools::file_ext(path)) == "png" && !.echogo_png_is_renderable(path)) {
      failures <- c(failures, paste0(src, " [invalid PNG]"))
    }
  }
  if (length(failures)) {
    stop(
      "Report references missing or invalid image assets:\n- ",
      paste(failures, collapse = "\n- "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Rarefaction benchmark (exploratory/consensus) with priority ordering & plateau
#'
#' Builds the same cumulative and permutation rarefaction curves as
#' \code{run_rarefaction_curves()}, but:
#' \itemize{
#'   \item pins a \strong{priority species list} first (e.g., core top models),
#'   \item adds remaining species via \strong{greedy contribution} (stepwise, ALL ontology),
#'   \item computes a \strong{statistical plateau} per panel (window gain + zero-streak; optional segmented + knee),
#'   \item writes \code{species_order_final_*.csv} and \code{plateau_detection_summary_*.csv},
#'         and draws a dashed vertical line at the suggested stop in each panel.
#' }
#'
#' It \emph{does not} re-run g:Profiler; it consumes your already-computed consensus table
#' (and gracefully handles species with zero results by treating missing columns as 0/FALSE).
#'
#' @param consensus_file Path to consensus enrichment Excel/CSV for this run
#'   (e.g., \code{.../consensus_enrichment_exploration_clean.xlsx}).
#' @param output_dir Directory where outputs are written.
#' @param species_priority Character vector of species labels to pin first
#'   (e.g., \code{c("hsapiens","mmusculus","rnorvegicus","drerio","dmelanogaster","celegans","ggallus","xtropicalis")}).
#'   Labels must match the column suffixes used in the consensus (e.g., \code{in_<label>_nobg}).
#' @param use_modes Character vector of modes to build. Any of \code{c("true_consensus","exploratory")}.
#'   Default \code{"exploratory"}.
#' @param min_steps Minimum steps before testing plateau (default 5).
#' @param window_k Moving window (steps) for gain test (default 3).
#' @param min_abs_gain Absolute new-term gain over last \code{window_k} steps below which we stop (default 5).
#' @param min_rel_gain Relative gain over last \code{window_k} steps (sum/new total) below which we stop (default 0.01).
#' @param zero_streak_k Stop if the last \code{k} steps each added 0 terms (default 3).
#'
#' @return Invisibly returns a list with:
#'   \item{species_order}{final order used (named by mode)}
#'   \item{plateau}{per-ontology stop summaries (named by mode)}
#'   \item{paths}{written artifact paths}
#'
#' @export
#' @keywords internal
#' @importFrom dplyr %>% filter mutate select across arrange desc pull bind_rows n_distinct
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_detect str_remove_all
#' @importFrom purrr map map_dfr
#' @importFrom openxlsx read.xlsx write.xlsx
#' @importFrom readr read_csv write_csv
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon labs theme_minimal element_text ggsave facet_wrap
run_rarefaction_benchmark <- function(
    consensus_file,
    output_dir,
    species_priority = c("hsapiens","mmusculus","rnorvegicus","drerio","dmelanogaster","celegans","ggallus","xtropicalis"),
    use_modes = c("exploratory"),
    min_steps = 5L,
    window_k = 3L,
    min_abs_gain = 5L,
    min_rel_gain = 0.01,
    zero_streak_k = 3L
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  comparison_dir <- output_dir
  suppressWarnings(dir.create(file.path(comparison_dir, "exploratory_no_bg"), recursive = TRUE, showWarnings = FALSE))

  # ---- Load consensus (same guards as your working version) ----
  if (!file.exists(consensus_file)) stop("consensus_file not found: ", consensus_file)
  consensus <- if (grepl("\\.xlsx$", consensus_file, ignore.case = TRUE)) {
    openxlsx::read.xlsx(consensus_file)
  } else {
    readr::read_csv(consensus_file, show_col_types = FALSE)
  } %>%
    dplyr::filter(.data$significant_in_any == TRUE, .data$ontology %in% c("BP","MF","CC"))

  # presence/absence transform (unchanged)
  species_cols <- grep("^in_", names(consensus), value = TRUE)
  species_cols <- setdiff(species_cols, "in_goseq")
  consensus_pa <- consensus %>%
    dplyr::select(.data$term_id, .data$ontology, dplyr::all_of(species_cols), .data$in_goseq) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(species_cols), ~ as.integer(. == TRUE), .names = "pa_{.col}")) %>%
    dplyr::mutate(pa_in_goseq = as.integer(.data$in_goseq == TRUE))
  pa_species_cols <- grep("^pa_in_", names(consensus_pa), value = TRUE)
  consensus_pa[pa_species_cols] <- lapply(consensus_pa[pa_species_cols], function(x) ifelse(is.na(x), 0, x))

  # species_info table (unchanged)
  species_info <- tibble::tibble(col = pa_species_cols) %>%
    dplyr::mutate(
      species      = stringr::str_remove_all(.data$col, "^pa_in_|_nobg$"),
      nobg         = stringr::str_detect(.data$col, "_nobg$"),
      species_mode = .data$col
    ) %>%
    dplyr::filter(.data$species != "goseq")

  # helper: build per-species term sets beyond GOseq baseline (optionally per-ontology)
  build_termsets <- function(entries, ont = "ALL") {
    if (ont == "ALL") {
      df <- consensus_pa
    } else {
      df <- consensus_pa %>% dplyr::filter(.data$ontology == ont)
    }
    base_terms <- df %>% dplyr::filter(.data$pa_in_goseq == 1) %>% dplyr::pull(.data$term_id) %>% unique()
    sets <- purrr::map(entries$species_mode, ~ {
      colname <- .
      df %>%
        dplyr::filter(.data[[colname]] == 1, .data$pa_in_goseq == 0) %>%
        dplyr::pull(.data$term_id) %>% unique()
    })
    names(sets) <- entries$species
    list(base = base_terms, sets = sets)
  }

  # ---- ORDERING per mode: priority fixed, then GREEDY remainder by marginal contribution (ALL) ----
  greedy_order <- function(entries) {
    tsALL <- build_termsets(entries, "ALL")
    # lock only priority species present in this mode
    pri <- species_priority[species_priority %in% names(tsALL$sets)]
    remaining <- setdiff(names(tsALL$sets), pri)
    order_final <- pri
    cumulative <- tsALL$base
    while (length(remaining)) {
      gains <- sapply(remaining, function(sp) length(setdiff(tsALL$sets[[sp]], cumulative)))
      sp_pick <- remaining[which.max(gains)]
      cumulative <- union(cumulative, tsALL$sets[[sp_pick]])
      order_final <- c(order_final, sp_pick)
      remaining <- setdiff(remaining, sp_pick)
    }
    order_final
  }

  # ---- Plateau detection helpers (per panel) ----
  plateau_window <- function(curve,
                             window_sz    = window_k,
                             min_abs      = min_abs_gain,
                             min_rel      = min_rel_gain,
                             zero_streak  = zero_streak_k,
                             min_steps_sz = min_steps) {
    if (nrow(curve) <= min_steps_sz)
      return(list(stop_at = NA_integer_, reason = "insufficient_steps"))

    for (i in seq(min_steps_sz, nrow(curve) - 1L)) {
      last_idxs <- (i - window_sz + 1L):i
      last_idxs <- last_idxs[last_idxs > 0]
      new_terms <- curve$n_new_terms[last_idxs]
      sum_gain  <- sum(new_terms, na.rm = TRUE)
      rel_gain  <- ifelse(curve$n_terms[i + 1L] > 0, sum_gain / curve$n_terms[i + 1L], 0)
      zs <- length(new_terms) >= zero_streak && all(utils::tail(new_terms, zero_streak) == 0)
      if (zs || sum_gain < min_abs || rel_gain < min_rel) {
        return(list(
          stop_at = i,
          reason  = sprintf("window[zs=%s,sum=%d,rel=%.3f]", zs, sum_gain, rel_gain)
        ))
      }
    }
    list(stop_at = NA_integer_, reason = "no_plateau_by_window")
  }
  plateau_segmented <- function(curve) {
    if (!requireNamespace("segmented", quietly = TRUE)) return(list(stop_at = NA_integer_, reason = "segmented_not_installed"))
    df <- curve[curve$step >= 0, , drop = FALSE]
    m0 <- stats::lm(n_terms ~ step, data = df)
    seg <- try(segmented::segmented(m0, seg.Z = ~ step, npsi = 1), silent = TRUE)
    if (inherits(seg, "try-error")) return(list(stop_at = NA_integer_, reason = "segmented_failed"))
    bp <- try(as.numeric(segmented::broken.line(seg)$psi[, "Est."]), silent = TRUE)
    if (inherits(bp, "try-error") || !length(bp)) return(list(stop_at = NA_integer_, reason = "segmented_no_bp"))
    list(stop_at = as.integer(round(bp[1])), reason = "segmented_breakpoint")
  }
  plateau_knee <- function(curve, span = 0.3) {
    if (nrow(curve) < 5L) return(list(stop_at = NA_integer_, reason = "too_short_for_knee"))
    lo <- stats::loess(n_terms ~ step, data = curve, span = span)
    xs <- curve$step; ys <- stats::predict(lo, xs)
    d1 <- diff(ys); d2 <- diff(d1)
    k  <- which.min(d2) + 1L
    list(stop_at = as.integer(k), reason = "knee_second_difference")
  }
  pick_stop <- function(curve) {
    a <- plateau_window(curve); b <- plateau_segmented(curve); c <- plateau_knee(curve)
    # choose median of available methods to be robust
    stops <- c(a$stop_at, b$stop_at, c$stop_at); stops <- stops[!is.na(stops) & stops > 0]
    sel <- if (length(stops)) as.integer(stats::median(stops)) else NA_integer_
    list(stop = sel, methods = list(window = a, segmented = b, knee = c))
  }

  # ---- Core builder (your original plotting style, with greedy order & per-panel vlines) ----
  build_mode <- function(bg_mode = TRUE, species_order = NULL) {
    tag <- if (bg_mode) "true_consensus" else "exploratory"
    out_dir <- if (bg_mode) comparison_dir else file.path(comparison_dir, "exploratory_no_bg")
    entries <- species_info %>% dplyr::filter(.data$nobg == !bg_mode)
    if (!nrow(entries)) return(invisible(NULL))

    # Determine order (priority fixed + greedy remainder) if not provided
    if (is.null(species_order)) species_order <- greedy_order(entries)

    # Save order with source labels
    readr::write_csv(
      tibble::tibble(
        order = seq_along(species_order),
        species = species_order,
        source = ifelse(species_order %in% species_priority, "priority", "greedy_by_contribution")
      ),
      file.path(out_dir, paste0("species_order_final_", tag, ".csv"))
    )

    # Build per-ontology term sets once
    tsets <- list(
      BP  = build_termsets(entries, "BP"),
      MF  = build_termsets(entries, "MF"),
      CC  = build_termsets(entries, "CC"),
      ALL = build_termsets(entries, "ALL")
    )

    mk_curve <- function(ont) {
      base <- tsets[[ont]]$base
      sets <- tsets[[ont]]$sets
      cumulative <- base
      curve <- tibble::tibble(species = "GOseq", ontology = ont,
                              n_terms = length(cumulative), n_new_terms = length(base),
                              step = 0L)
      for (i in seq_along(species_order)) {
        sp <- species_order[i]
        add_now <- setdiff(sets[[sp]], cumulative)
        cumulative <- union(cumulative, add_now)
        curve <- dplyr::bind_rows(curve, tibble::tibble(
          species = sp, ontology = ont,
          n_terms = length(cumulative), n_new_terms = length(add_now),
          step = i
        ))
      }
      curve
    }

    curve_BP  <- mk_curve("BP")
    curve_MF  <- mk_curve("MF")
    curve_CC  <- mk_curve("CC")
    curve_ALL <- mk_curve("ALL")
    cumulative_df <- dplyr::bind_rows(curve_BP, curve_MF, curve_CC, curve_ALL)

    # ---- Plateau per panel (BP/MF/CC/ALL) ----
    stops_list <- list(
      BP  = pick_stop(curve_BP),
      MF  = pick_stop(curve_MF),
      CC  = pick_stop(curve_CC),
      ALL = pick_stop(curve_ALL)
    )
    stop_tbl <- tibble::tibble(
      ontology = names(stops_list),
      stop_step = sapply(stops_list, function(x) x$stop),
      method_window = sapply(stops_list, function(x) x$methods$window$reason),
      method_segmented = sapply(stops_list, function(x) x$methods$segmented$reason),
      method_knee = sapply(stops_list, function(x) x$methods$knee$reason)
    )
    readr::write_csv(stop_tbl, file.path(out_dir, paste0("plateau_detection_summary_", tag, ".csv")))

    # ---- Exports & plots (same look; dashed vline per panel using that panel's stop) ----
    if (nrow(cumulative_df) > 0) {
      openxlsx::write.xlsx(
        cumulative_df,
        file.path(out_dir, paste0("cumulative_by_origin_", tag, ".xlsx")),
        asTable = TRUE
      )

      ontology_cols <- c(BP = "#1b9e77", MF = "#d95f02", CC = "#7570b3", ALL = "#e7298a")
      mk_panel <- function(df, ont) {
        df_ont <- df %>% dplyr::filter(.data$ontology == ont)
        if (all(df_ont$species == "GOseq")) {
          df_ont <- df_ont %>% dplyr::mutate(species = factor(.data$species, levels = "GOseq"))
        } else {
          species_levels <- df_ont %>%
            dplyr::filter(.data$species != "GOseq") %>%
            dplyr::arrange(.data$step, .data$species) %>%
            dplyr::pull(.data$species) %>% unique()
          species_levels <- c("GOseq", species_levels)
          df_ont <- df_ont %>%
            dplyr::mutate(species = factor(.data$species, levels = species_levels)) %>%
            dplyr::arrange(.data$step)
        }
        p <- ggplot2::ggplot(df_ont, ggplot2::aes(x = .data$species, y = .data$n_terms)) +
          ggplot2::geom_line(group = 1, linewidth = 1.2, color = ontology_cols[[ont]]) +
          ggplot2::geom_point(size = 3, color = ontology_cols[[ont]]) +
          ggplot2::labs(
            title = paste("Cumulative GO Term Curve -", ont),
            x = "Species Added (priority fixed, greedy remainder)", y = "Cumulative Unique GO Terms"
          ) +
          ggplot2::theme_minimal(base_size = 13) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
        s <- stop_tbl$stop_step[stop_tbl$ontology == ont]
        if (length(s) == 1L && !is.na(s) && s > 0) {
          p <- p + ggplot2::geom_vline(xintercept = s + 1, linetype = 2) # +1: after GOseq tick
        }
        p
      }

      plots <- list(
        BP  = mk_panel(cumulative_df, "BP"),
        MF  = mk_panel(cumulative_df, "MF"),
        CC  = mk_panel(cumulative_df, "CC"),
        ALL = mk_panel(cumulative_df, "ALL")
      )
      combo <- plots$BP + plots$MF + plots$CC + plots$ALL +
        patchwork::plot_layout(ncol = 2) +
        patchwork::plot_annotation(
          title = if (bg_mode) "Cumulative GO Term Curves (True Consensus)"
          else                 "Cumulative GO Term Curves (Exploratory)"
        )

      ggplot2::ggsave(file.path(out_dir, paste0("cumulative_curve_combined_", tag, ".pdf")),
                      combo, width = 14, height = 10)
    }

    # ---- Permutation rarefaction (unchanged visuals) ----
    set.seed(1L)
    n_perm <- 1000L
    rarefied_all <- list()
    for (ont in c("BP","MF","CC","ALL")) {
      base_terms <- if (identical(ont, "ALL")) tsets$ALL$base else {
        consensus_pa %>% dplyr::filter(.data$ontology == ont, .data$pa_in_goseq == 1) %>%
          dplyr::pull(.data$term_id) %>% unique()
      }
      term_sets <- if (identical(ont, "ALL")) tsets$ALL$sets else {
        entries$species_mode %>%
          stats::setNames(stringr::str_remove_all(entries$species_mode, "^pa_in_|_nobg$")) %>%
          purrr::map(function(colname) {
            consensus_pa %>%
              dplyr::filter(if (ont != "ALL") .data$ontology == ont else TRUE,
                            .data[[colname]] == 1, .data$pa_in_goseq == 0) %>%
              dplyr::pull(.data$term_id) %>% unique()
          })
      }
      species_list <- names(term_sets); if (!length(species_list)) next

      rarem <- replicate(n_perm, {
        sp_order <- sample(species_list)
        cumulative <- base_terms
        cumsum <- numeric(length(sp_order))
        for (i in seq_along(sp_order)) {
          new_terms <- setdiff(term_sets[[sp_order[i]]], cumulative)
          cumulative <- union(cumulative, new_terms)
          cumsum[i] <- length(cumulative)
        }
        cumsum
      }, simplify = "matrix")

      stats_df <- tibble::tibble(
        species_index = seq_along(species_list) + 1L,
        mean_terms = rowMeans(rarem),
        ci_lower   = apply(rarem, 1, stats::quantile, probs = 0.025),
        ci_upper   = apply(rarem, 1, stats::quantile, probs = 0.975),
        ontology   = ont
      )
      rarefied_all[[ont]] <- stats_df
    }

    if (length(rarefied_all)) {
      rarefied_df <- dplyr::bind_rows(rarefied_all)
      openxlsx::write.xlsx(rarefied_df,
                           file.path(out_dir, paste0("rarefaction_by_origin_", tag, ".xlsx")),
                           asTable = TRUE)

      if (nrow(rarefied_df) > 0) {
        p2 <- ggplot2::ggplot(rarefied_df, ggplot2::aes(x = .data$species_index, y = .data$mean_terms, color = .data$ontology, group = .data$ontology)) +
          ggplot2::geom_line(linewidth = 1.2) +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper, fill = .data$ontology), alpha = 0.2, color = NA) +
          ggplot2::facet_wrap(~ontology, scales = "free_y") +
          ggplot2::labs(
            title = if (bg_mode) "Rarefaction Curve (1000 Perms): GOseq + g:Profiler (True Consensus)"
            else                 "Rarefaction Curve (1000 Perms): GOseq + All Significant Terms (Exploratory)",
            x = "Species Randomly Added (Index)", y = "Mean Cumulative Unique GO Terms"
          ) +
          ggplot2::theme_minimal(base_size = 13)
        ggplot2::ggsave(file.path(out_dir, paste0("rarefaction_curve_by_origin_", tag, ".pdf")),
                        p2, width = 10, height = 6)
      }
    }

    # presence/absence table
    openxlsx::write.xlsx(consensus_pa,
                         file = file.path(out_dir, "presence_absence_table.xlsx"),
                         asTable = TRUE)

    invisible(list(
      species_order = species_order,
      stops = stop_tbl
    ))
  } # end build_mode

  out <- list(paths = list())
  if ("true_consensus" %in% use_modes) {
    out_bg <- build_mode(TRUE, species_order = NULL)  # greedy order for BG mode
    out$species_order$true_consensus <- out_bg$species_order
    out$plateau$true_consensus <- out_bg$stops
  }
  if ("exploratory" %in% use_modes) {
    out_nb <- build_mode(FALSE, species_order = NULL) # greedy order for noBG mode
    out$species_order$exploratory <- out_nb$species_order
    out$plateau$exploratory <- out_nb$stops
  }

  invisible(out)
}


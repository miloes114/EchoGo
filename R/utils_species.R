# ---- file: R/utils_species.R -----------------------------------------------
# File-level imports (put these at the very top of the file; one block only)
#' @importFrom dplyr mutate select distinct arrange bind_rows left_join if_else
#' @importFrom dplyr .data
#' @importFrom tibble tibble
#' @importFrom utils data head
#' @importFrom readr read_csv write_csv
#' @importFrom xml2 read_html
#' @importFrom rvest html_elements html_table html_attr html_text2
#' @importFrom httr GET RETRY content timeout user_agent accept_json status_code
#' @importFrom stringr str_match_all
#' @importFrom rappdirs user_cache_dir
#' @importFrom jsonlite fromJSON write_json read_json
NULL

# small infix
`%||%` <- function(a,b) if (!is.null(a)) a else b

# ---- Cache helpers & metadata ----------------------------------------------

.echogo_cache_dir    <- function() rappdirs::user_cache_dir("EchoGO")
.echogo_species_csv  <- function() file.path(.echogo_cache_dir(), "gprofiler_species.csv")
.echogo_species_meta <- function() file.path(.echogo_cache_dir(), "gprofiler_species_meta.json")

.echogo_write_meta <- function(source, ok = TRUE) {
  dir.create(.echogo_cache_dir(), recursive = TRUE, showWarnings = FALSE)
  meta <- list(
    source    = source,
    timestamp = as.integer(Sys.time()),
    ok        = isTRUE(ok),
    version   = as.character(utils::packageVersion("EchoGO"))
  )
  jsonlite::write_json(meta, .echogo_species_meta(), auto_unbox = TRUE)
}

.echogo_read_meta <- function() {
  f <- .echogo_species_meta()
  if (!file.exists(f)) return(NULL)
  tryCatch(jsonlite::read_json(f, simplifyVector = TRUE), error = function(e) NULL)
}

.echogo_is_stale <- function(max_age_days = getOption("EchoGO.species_max_age_days", 30L)) {
  m <- .echogo_read_meta(); if (is.null(m)) return(TRUE)
  age_days <- (as.integer(Sys.time()) - as.integer(m$timestamp)) / 86400
  isTRUE(age_days > max_age_days)
}

# provenance tagger
.set_prov <- function(df, src) {
  attr(df, "echogo_source")    <- src
  attr(df, "echogo_timestamp") <- Sys.time()
  df
}

# gatekeeper: seed cache from packaged data, try online refresh if allowed/stale
#' Ensure a fresh-enough species cache (offline-first)
#' @return invisibly TRUE
#' @export
echogo_ensure_species_cache <- function(force = FALSE) {
  dir.create(.echogo_cache_dir(), recursive = TRUE, showWarnings = FALSE)
  csv   <- .echogo_species_csv()
  auto  <- isTRUE(getOption("EchoGO.species_autoupdate", TRUE))
  stale <- .echogo_is_stale()

  # If cache exists and neither forced nor stale (when auto-update), keep it.
  if (file.exists(csv) && !(force || (auto && stale))) return(invisible(TRUE))

  # Try online refresh (non-fatal if offline)
  if (auto || force) {
    df <- try(.echogo_fetch_species_online(), silent = TRUE)
    if (!inherits(df, "try-error") && is.data.frame(df) && nrow(df)) {
      readr::write_csv(df, csv)
      .echogo_write_meta("api")
      return(invisible(TRUE))
    }
  }

  # Seed from packaged fallback if cache missing
  if (!file.exists(csv)) {
    fb <- .echogo_species_fallback()
    if (nrow(fb)) {
      readr::write_csv(fb, csv)
      .echogo_write_meta("packaged")
    }
  }
  invisible(TRUE)
}

# ---- Core species table fetchers --------------------------------------------

# packaged fallback loader
.echogo_species_fallback <- function() {
  fb <- tryCatch({
    utils::data("echogo_species", package = "EchoGO", envir = environment())
    get("echogo_species", envir = environment())
  }, error = function(e) tibble::tibble())
  .normalize_species_tbl(fb)
}

# online fetcher (API preferred, HTML fallback)
.echogo_fetch_species_online <- function() {
  base_timeout <- getOption("EchoGO.net_timeout_sec", 25L)
  bases <- c(
    getOption("EchoGO.gprofiler_base", "https://biit.cs.ut.ee/gprofiler"),
    "https://biit.cs.ut.ee/gprofiler_beta"
  )
  # A) Try API first
  for (b in unique(bases)) {
    api_df <- try(.fetch_api_species(b, timeout_sec = base_timeout), silent = TRUE)
    if (!inherits(api_df, "try-error") && is.data.frame(api_df) && nrow(api_df)) {
      return(.set_prov(api_df, "api"))
    }
  }
  # B) HTML fallback last
  live <- .fetch_html_species()
  if (nrow(live)) return(.set_prov(live, "html"))
  tibble::tibble(organism=character(), name=character(), ncbi=integer(), alias=character())
}

# standardizer
.normalize_species_tbl <- function(x) {
  if (!is.data.frame(x) || !nrow(x)) {
    return(tibble::tibble(
      organism = character(), name = character(), ncbi = integer(), alias = character()
    ))
  }
  x <- x |>
    dplyr::mutate(
      organism = as.character(.data$organism %||% NA_character_),
      name     = as.character(.data$name     %||% NA_character_),
      ncbi     = suppressWarnings(as.integer(.data$ncbi %||% NA_integer_)),
      alias    = as.character(.data$alias    %||% NA_character_)
    ) |>
    dplyr::select(organism, name, ncbi, alias) |>
    dplyr::distinct(organism, .keep_all = TRUE) |>
    dplyr::arrange(.data$organism)
  x[!is.na(x$organism) & nzchar(x$organism), , drop = FALSE]
}

# API branch
.fetch_api_species <- function(base_url, timeout_sec = 25L) {
  dbg  <- isTRUE(getOption("EchoGO.debug"))
  base <- sub("/+$", "", base_url %||% "https://biit.cs.ut.ee/gprofiler")
  paths <- c("/api/util/organisms_list", "/api/util/organisms_list/")
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    if (dbg) message("[EchoGO] jsonlite not available; skipping API branch.")
    return(tibble::tibble())
  }
  for (p in paths) {
    api <- paste0(base, p)
    if (dbg) message("[EchoGO] fetch_api(): GET ", api)
    resp <- try(httr::RETRY("GET", api,
                            httr::user_agent("EchoGO/1.0 (+https://github.com/echogo/echogo)"),
                            httr::accept_json(),
                            httr::timeout(as.numeric(timeout_sec)),
                            times = 3),
                silent = TRUE)
    if (inherits(resp, "try-error")) next
    if (httr::status_code(resp) != 200) next
    txt <- try(httr::content(resp, as = "text", encoding = "UTF-8"), silent = TRUE)
    if (inherits(txt, "try-error") || !nzchar(txt)) next

    j <- try(jsonlite::fromJSON(txt, simplifyVector = TRUE), silent = TRUE)
    if (inherits(j, "try-error") || !is.data.frame(j) || !nrow(j)) next

    cn <- tolower(names(j))
    id_col    <- if ("organism" %in% cn) "organism" else if ("id" %in% cn) "id" else NA_character_
    name_col  <- if ("name" %in% cn) "name" else if ("display_name" %in% cn) "display_name"
    else if ("scientific_name" %in% cn) "scientific_name" else NA_character_
    ncbi_col  <- if ("ncbi" %in% cn) "ncbi" else if ("taxonomy_id" %in% cn) "taxonomy_id"
    else if ("taxon_id" %in% cn) "taxon_id" else NA_character_
    alias_col <- if ("alias" %in% cn) "alias" else if ("aliases" %in% cn) "aliases" else NA_character_
    if (is.na(id_col) || is.na(name_col)) next

    org   <- j[[id_col]]
    name  <- j[[name_col]]
    ncbi  <- if (!is.na(ncbi_col)) j[[ncbi_col]] else NA_integer_
    alias <- if (!is.na(alias_col)) j[[alias_col]] else NA_character_
    if (is.list(alias)) alias <- vapply(alias, function(x) paste0(x, collapse=";"), "")

    out <- tibble::tibble(
      organism = trimws(as.character(org)),
      name     = trimws(as.character(name)),
      ncbi     = suppressWarnings(as.integer(ncbi)),
      alias    = as.character(alias)
    )
    out <- out[!is.na(out$organism) & nzchar(out$organism), , drop = FALSE]
    if (nrow(out)) return(.normalize_species_tbl(out))
  }
  if (dbg) message("[EchoGO] API branch yielded 0 rows; falling back to HTML.")
  tibble::tibble()
}

# HTML branch (very best-effort)
.fetch_html_species <- function() {
  url <- "https://biit.cs.ut.ee/gprofiler/page/organism-list"
  if (requireNamespace("xml2", quietly = TRUE) && requireNamespace("rvest", quietly = TRUE)) {
    doc <- tryCatch(xml2::read_html(url), error = function(e) NULL)
    if (!is.null(doc)) {
      nodes <- tryCatch(rvest::html_elements(doc, "[data-organism]"), error = function(e) NULL)
      if (!is.null(nodes) && length(nodes)) {
        org  <- rvest::html_attr(nodes, "data-organism")
        name <- tryCatch(rvest::html_text2(nodes), error = function(e) rvest::html_text(nodes))
        df <- tibble::tibble(
          organism = trimws(as.character(org)),
          name     = trimws(as.character(name)),
          ncbi     = NA_integer_,
          alias    = NA_character_
        )
        df <- df[!is.na(df$organism) & nzchar(df$organism), , drop = FALSE]
        if (nrow(df)) return(.normalize_species_tbl(df))
      }
      nodes_tbl <- tryCatch(rvest::html_elements(doc, "table"), error = function(e) NULL)
      if (!is.null(nodes_tbl) && length(nodes_tbl)) {
        tbs <- lapply(nodes_tbl, function(n)
          tryCatch(rvest::html_table(n, fill = TRUE), error = function(e) NULL))
        tbs <- unlist(tbs, recursive = FALSE)
        tbs <- Filter(function(z) is.data.frame(z) && nrow(z) > 0, tbs)
        if (length(tbs)) {
          pieces <- lapply(tbs, function(tb) {
            cn <- tolower(gsub("\\s+", "_", names(tb))); names(tb) <- cn
            nm  <- if ("display_name" %in% cn) tb[["display_name"]] else
              if ("scientific_name" %in% cn) tb[["scientific_name"]] else
                if ("name" %in% cn) tb[["name"]] else NA_character_
            org <- if ("id" %in% cn) tb[["id"]] else
              if ("organism" %in% cn) tb[["organism"]] else NA_character_
            tibble::tibble(
              organism = trimws(as.character(org)),
              name     = trimws(as.character(nm)),
              ncbi     = NA_integer_,
              alias    = NA_character_
            )
          })
          df <- dplyr::bind_rows(pieces)
          df <- df[!duplicated(df$organism), , drop = FALSE]
          if (nrow(df)) return(.normalize_species_tbl(df))
        }
      }
    }
  }
  # last-resort regex
  txt <- tryCatch(httr::content(httr::GET(url), as = "text", encoding = "UTF-8"),
                  error = function(e) "")
  if (!nzchar(txt)) return(tibble::tibble())
  m <- stringr::str_match_all(
    txt, "data-organism=\"([A-Za-z0-9_\\.\\-]+)\"[^>]*>\\s*([^<\\n]+)\\s*<"
  )[[1]]
  if (!nrow(m)) return(tibble::tibble())
  .normalize_species_tbl(tibble::tibble(
    organism = trimws(m[,2]), name = trimws(m[,3]),
    ncbi = NA_integer_, alias = NA_character_
  ))
}

# ---- Public: Get species table (cached; offline-first with fallback) --------

#' Get g:Profiler species table (cached, with packaged fallback + API/HTML scrape)
#'
#' Tries the JSON API first, then HTML, else packaged fallback. Seeds/refreshes
#' the user cache via \code{echogo_ensure_species_cache()} (non-blocking).
#'
#' @param refresh logical; if TRUE, ignores the user cache and re-fetches online.
#' @return tibble with columns: organism, name, ncbi, alias
#' @export
echogo_gprofiler_species <- function(refresh = FALSE) {
  cache_dir  <- .echogo_cache_dir()
  cache_file <- .echogo_species_csv()
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  # Ensure cache exists (offline-first; online if allowed and stale)
  if (!isTRUE(refresh)) echogo_ensure_species_cache(force = FALSE)

  read_cache <- function() {
    if (file.exists(cache_file)) {
      tryCatch(readr::read_csv(cache_file, show_col_types = FALSE),
               error = function(e) tibble::tibble())
    } else tibble::tibble()
  }

  # 1) Cached unless refresh
  if (!isTRUE(refresh)) {
    cached <- .normalize_species_tbl(read_cache())
    if (nrow(cached)) return(.set_prov(cached, "cache"))
  }

  # 2) Live online (API → HTML)
  online <- .echogo_fetch_species_online()
  if (nrow(online)) {
    try(readr::write_csv(online, cache_file), silent = TRUE)
    .echogo_write_meta(attr(online, "echogo_source") %||% "api")
    return(online)
  }

  # 3) Packaged fallback
  fb <- .echogo_species_fallback()
  .set_prov(fb, "packaged")
}

#' Refresh and persist the g:Profiler species list
#' @inherit echogo_gprofiler_species return
#' @export
echogo_update_species_cache <- function(update_package_data = FALSE, verbose = TRUE) {
  sp <- echogo_gprofiler_species(refresh = TRUE)
  if (verbose) message(sprintf("Fetched %s species (unique organism IDs).",
                               format(nrow(sp), big.mark=",")))
  if (!nrow(sp)) {
    warning("No species fetched; aborting.")
    return(invisible(sp))
  }

  if (isTRUE(update_package_data)) {
    pkg_root <- NULL
    if (requireNamespace("usethis", quietly = TRUE)) {
      pkg_root <- tryCatch(usethis::proj_get(), error = function(e) NULL)
    }
    if (is.null(pkg_root) && requireNamespace("rprojroot", quietly = TRUE)) {
      pkg_root <- tryCatch(rprojroot::find_package_root_file("."), error = function(e) NULL)
    }
    if (is.null(pkg_root)) stop("Could not locate package root. Run inside the EchoGO package project.")

    echogo_species <- sp
    data_dir <- file.path(pkg_root, "data")
    dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
    save_path <- file.path(data_dir, "echogo_species.rda")

    if (requireNamespace("usethis", quietly = TRUE)) {
      usethis::use_data(echogo_species, overwrite = TRUE, compress = "xz")
    } else {
      save(echogo_species, file = save_path, compress = "xz")
    }

    if (verbose) {
      message("Wrote packaged fallback: ",
              normalizePath(save_path, winslash = "/", mustWork = file.exists(save_path)))
      message("Re-run devtools::document(); devtools::install() to ship the new fallback.")
    }
  }

  invisible(sp)
}

# ---- Taxonomy enrichment (fallback + optional online) -----------------------

#' Enrich species table with taxonomic ranks (superkingdom → genus)
#' Uses packaged fallback (data/echogo_taxonomy_fallback.rda) and fills gaps
#' from a user cache (taxonomy_lineage.rds). If auto-update is enabled and
#' 'taxize' is installed, missing ranks are fetched from NCBI classification().
#' @export
echogo_enrich_taxonomy <- function(df, refresh = FALSE, rate_limit = 0.3, verbose = interactive()) {
  ranks <- c("superkingdom","kingdom","phylum","class","order","family","genus")

  # 1) packaged fallback join by ncbi (suppress warning if dataset not shipped yet)
  tax_fb <- suppressWarnings(tryCatch({
    utils::data("echogo_taxonomy_fallback", package = "EchoGO", envir = environment())
    get("echogo_taxonomy_fallback", envir = environment())
  }, error = function(e) NULL))

  if (!is.null(tax_fb) && is.data.frame(tax_fb) && nrow(tax_fb)) {
    df <- dplyr::left_join(df, tax_fb, by = "ncbi")
  } else {
    for (r in ranks) if (!r %in% names(df)) df[[r]] <- NA_character_
  }

  # 2) fill missing ranks from user cache + optional online (taxize)
  cache_file <- file.path(.echogo_cache_dir(), "taxonomy_lineage.rds")
  lin_cache  <- if (file.exists(cache_file)) readRDS(cache_file) else list()

  missing_ix <- which(!is.na(df$ncbi) &
                        !Reduce(`&`, lapply(ranks, function(r) !is.na(df[[r]]) & nzchar(df[[r]]))))
  todo_ids <- unique(df$ncbi[missing_ix]); todo_ids <- todo_ids[!is.na(todo_ids)]
  if (!length(todo_ids)) return(df)

  allow_net  <- isTRUE(getOption("EchoGO.species_autoupdate", TRUE)) &&
    isTRUE(getOption("EchoGO.taxonomy_online", TRUE))
  has_taxize <- requireNamespace("taxize", quietly = TRUE)
  has_key    <- nzchar(Sys.getenv("ENTREZ_KEY", ""))

  if (allow_net && has_taxize && has_key) {
    for (tx in setdiff(todo_ids, as.integer(names(lin_cache)))) {
      cl <- try(taxize::classification(tx, db = "ncbi"), silent = TRUE)
      if (!inherits(cl, "try-error")) {
        # unwrap list result if needed
        if (is.list(cl) && length(cl) > 0) cl <- cl[[1]]
        # keep only what we need
        if (is.data.frame(cl) && all(c("name","rank") %in% names(cl))) {
          lin_cache[[as.character(tx)]] <- cl[, c("name","rank")]
        }
      }
      if (rate_limit > 0) Sys.sleep(rate_limit)
    }
    saveRDS(lin_cache, cache_file)
  } else if (allow_net && has_taxize && !has_key && isTRUE(verbose)) {
    message("[EchoGO] Skipping online taxonomy fill: no ENTREZ_KEY set. ",
            "Set one with usethis::edit_r_environ() (ENTREZ_KEY=...) or Sys.setenv(ENTREZ_KEY='...'). ",
            "Or disable with options(EchoGO.taxonomy_online = FALSE).")
  }


  fill_one <- function(tx, r) {
    cl <- lin_cache[[as.character(tx)]]
    if (is.null(cl)) return(NA_character_)
    v <- cl$name[match(r, cl$rank)]
    ifelse(is.na(v), NA_character_, v)
  }
  for (r in ranks) {
    miss <- which(is.na(df[[r]]) | !nzchar(df[[r]]))
    if (length(miss)) df[[r]][miss] <- vapply(df$ncbi[miss], fill_one, character(1), r = r)
  }
  df
}


# ---- Tags: curated panels from YAML -----------------------------------------

#' Add curated tag panels from inst/extdata/species_tags.yaml
#' @export
echogo_add_tags <- function(df, yaml_path = system.file("extdata","species_tags.yaml", package="EchoGO")) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    df$tags <- I(vector("list", nrow(df))); return(df)
  }
  tags <- tryCatch(yaml::read_yaml(yaml_path), error = function(e) NULL)
  if (is.null(tags)) { df$tags <- I(vector("list", nrow(df))); return(df) }

  tag_map <- lapply(tags, function(v) unique(as.character(v)))
  revmap <- list()
  for (tg in names(tag_map)) for (id in tag_map[[tg]]) revmap[[id]] <- c(revmap[[id]], tg)
  df$tags <- I(lapply(df$organism, function(id) unique(revmap[[id]] %||% character())))
  df
}

# ---- Unified public builder + interactive view ------------------------------

#' Build the augmented species table (cache → taxonomy → tags)
#' @export
echogo_species_table <- function(refresh = FALSE, include_taxonomy = TRUE, include_tags = TRUE) {
  echogo_ensure_species_cache(force = isTRUE(refresh))
  sp <- echogo_gprofiler_species(refresh = isTRUE(refresh))
  if (include_taxonomy) sp <- echogo_enrich_taxonomy(sp, refresh = refresh)
  if (include_tags)     sp <- echogo_add_tags(sp)
  sp
}

#' List supported g:Profiler species (compact/interactive)
#' @param view logical; if TRUE, shows an interactive table when DT is available.
#' @param refresh logical; passed to [echogo_species_table()].
#' @param n integer; number of rows to show in console preview.
#' @return tibble of species (invisibly when `view = TRUE`)
#' @export
echogo_list_species <- function(view = interactive(), refresh = FALSE, n = 30) {
  df <- echogo_species_table(refresh = refresh)

  if (!nrow(df)) {
    message("No species available. Try echogo_gprofiler_species(refresh = TRUE) ",
            "or ensure a packaged fallback data/echogo_species.rda exists.")
    return(invisible(df))
  }

  if (isTRUE(view) && requireNamespace("DT", quietly = TRUE)) {
    show <- df
    show$tags_display <- vapply(show$tags, function(x) paste(x, collapse = ", "), "")
    keep <- c("organism","name","ncbi","superkingdom","kingdom","phylum",
              "class","order","family","genus","tags_display")
    show <- show[, intersect(keep, names(show)), drop = FALSE]

    DT::datatable(
      show,
      rownames   = FALSE,
      extensions = c("Buttons","FixedHeader","ColReorder"),
      options = list(
        pageLength = 15, dom = "Bfrtip",
        buttons = c("copy","csv","excel","colvis"),
        fixedHeader = TRUE, colReorder = TRUE,
        searchHighlight = TRUE
      ),
      caption = htmltools::tags$caption(
        style = 'caption-side: top; text-align: left;',
        "Tip: filter ranks via column search; copy organism IDs from the first column."
      )
    ) |>
      DT::formatStyle("organism", `font-weight` = "bold")
  } else {
    print(utils::head(df, n))
  }
  invisible(df)
}

# ---- Search / validate / pick (kept; augmented below) -----------------------

#' Search EchoGO's species list by name/ID/alias/NCBI
#' @export
echogo_find_species <- function(query,
                                fields = c("organism","name","alias","ncbi"),
                                fuzzy = TRUE,
                                refresh = FALSE) {
  stopifnot(is.character(query), length(query) >= 1)
  fields <- intersect(fields, c("organism","name","alias","ncbi"))
  if (!length(fields)) fields <- c("organism","name","alias","ncbi")

  df <- echogo_gprofiler_species(refresh = refresh)
  if (!nrow(df)) return(df)

  col_to_char <- function(x) if (is.numeric(x)) format(x, trim = TRUE, scientific = FALSE) else as.character(x %||% "")
  L <- lapply(fields, function(f) tolower(col_to_char(df[[f]])))
  names(L) <- fields

  hits <- logical(nrow(df))
  for (q in query) {
    q <- tolower(q)
    rx <- vapply(L, function(col) grepl(q, col, fixed = TRUE) | grepl(q, col), logical(nrow(df)))
    hits <- hits | apply(rx, 1, any)

    if (isTRUE(fuzzy) && !grepl("^\\s*\\d+\\s*$", q)) {
      approx_idx <- unique(unlist(lapply(L, function(col) {
        which(vapply(col, function(x) length(base::agrep(q, x, ignore.case = TRUE, max.distance = 0.1)) > 0, logical(1)))
      })))
      if (length(approx_idx)) hits[approx_idx] <- TRUE
    }
  }
  out <- df[hits, , drop = FALSE]
  dplyr::arrange(out, .data$organism)
}

#' Validate and/or suggest species IDs for `run_full_echogo()`
#' @export
echogo_validate_species <- function(species, refresh = FALSE) {
  stopifnot(is.character(species))
  df <- echogo_gprofiler_species(refresh = refresh)
  if (!nrow(df)) {
    return(tibble::tibble(
      input = species, status = "unknown",
      organism = NA_character_, name = NA_character_, ncbi = NA_integer_,
      alias = NA_character_, suggestion = NA_character_
    ))
  }

  exact_idx <- match(species, df$organism)
  status <- ifelse(!is.na(exact_idx), "ok", "unknown")

  lower_fields <- list(
    organism = tolower(df$organism),
    name     = tolower(df$name    %||% ""),
    alias    = tolower(df$alias   %||% "")
  )

  suggest_one <- function(x) {
    x0 <- tolower(x)
    sw <- startsWith(lower_fields$organism, x0)
    if (any(sw)) return(df$organism[which(sw)[1]])
    sub_name  <- grepl(x0, lower_fields$name,  fixed = TRUE)
    sub_alias <- grepl(x0, lower_fields$alias, fixed = TRUE)
    if (any(sub_name | sub_alias)) return(df$organism[which(sub_name | sub_alias)[1]])
    score <- function(candidate) as.numeric(utils::adist(x0, candidate, ignore.case = TRUE))
    cand_vec <- unique(c(lower_fields$organism, lower_fields$name, lower_fields$alias))
    cand_vec <- cand_vec[nzchar(cand_vec)]
    if (!length(cand_vec)) return(NA_character_)
    best <- cand_vec[which.min(vapply(cand_vec, score, numeric(1)))]
    row <- which(lower_fields$organism == best |
                   lower_fields$name == best |
                   lower_fields$alias == best)[1]
    if (length(row) && !is.na(row)) df$organism[row] else NA_character_
  }

  suggestion <- vapply(seq_along(species), function(i) {
    if (status[i] == "ok") return(df$organism[exact_idx[i]])
    suggest_one(species[i])
  }, character(1))

  row_pick <- ifelse(!is.na(exact_idx), exact_idx, match(suggestion, df$organism))
  out <- tibble::tibble(
    input  = species,
    status = ifelse(!is.na(exact_idx), "ok",
                    ifelse(!is.na(suggestion), "suggested", "unknown")),
    organism = dplyr::if_else(!is.na(row_pick), df$organism[row_pick], NA_character_),
    name     = dplyr::if_else(!is.na(row_pick), df$name[row_pick], NA_character_),
    ncbi     = dplyr::if_else(!is.na(row_pick), df$ncbi[row_pick], NA_integer_),
    alias    = dplyr::if_else(!is.na(row_pick), df$alias[row_pick], NA_character_),
    suggestion = dplyr::if_else(status == "ok", NA_character_, suggestion)
  )
  dplyr::arrange(out, factor(.data$status, levels = c("ok","suggested","unknown")), .data$input)
}

#' Pick species interactively (always shows a view)
#' @export
echogo_pick_species <- function(refresh = FALSE, page_len = 25, force_browser = FALSE) {
  df <- echogo_species_table(refresh = refresh)
  if (!nrow(df)) {
    message("No species available. Try echogo_gprofiler_species(refresh = TRUE).")
    return(invisible(df))
  }

  if (interactive() && requireNamespace("DT", quietly = TRUE)) {
    show <- df
    show$tags_display <- vapply(show$tags, function(x) paste(x, collapse = ", "), "")
    keep <- c("organism","name","ncbi","superkingdom","kingdom","phylum",
              "class","order","family","genus","tags_display")
    show <- show[, intersect(keep, names(show)), drop = FALSE]

    dt <- DT::datatable(
      show,
      rownames   = FALSE,
      extensions = c("Buttons","FixedHeader","ColReorder"),
      options = list(
        pageLength = page_len,
        dom = "Bfrtip",
        buttons = c("copy","csv","excel","colvis"),
        fixedHeader = TRUE,
        colReorder  = TRUE,
        searchHighlight = TRUE
      ),
      caption = htmltools::tags$caption(
        style = 'caption-side: top; text-align: left;',
        "Tip: search by name, ID, NCBI, or filter by taxonomy; copy IDs from 'organism'."
      )
    )
    try(print(dt), silent = TRUE)

    if (isTRUE(force_browser) && requireNamespace("htmlwidgets", quietly = TRUE)) {
      tmp <- tempfile(fileext = ".html")
      htmlwidgets::saveWidget(dt, tmp, selfcontained = TRUE)
      utils::browseURL(tmp)
    }

  } else if (interactive()) {
    message("DT not available; opening a spreadsheet view.")
    try(utils::View(df), silent = TRUE)

  } else {
    message("Non-interactive session. Showing a preview; install DT for an interactive table.")
    print(utils::head(df, 30))
  }

  invisible(df)
}

#' Check species argument and provide actionable feedback (preflight)
#' @export
echogo_preflight_species <- function(species, refresh = FALSE, error_if_unknown = TRUE) {
  v <- echogo_validate_species(species, refresh = refresh)
  adoptable <- which(v$status == "suggested" & !is.na(v$suggestion))
  if (length(adoptable)) {
    message("[EchoGO] Interpreting: ",
            paste0(v$input[adoptable], "→", v$suggestion[adoptable], collapse = ", "))
    v$organism[adoptable] <- v$suggestion[adoptable]
    v$status[adoptable] <- "ok"
  }
  unknown <- v$input[v$status == "unknown"]
  if (length(unknown)) {
    msg <- paste0(
      "Unknown species input: ", paste(unknown, collapse = ", "),
      "\nTry: echogo_find_species('<name or ID>') or call echogo_pick_species()."
    )
    if (isTRUE(error_if_unknown)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  unique(v$organism[v$status == "ok"])
}

# ---- Tag/taxonomy aware selection helpers -----------------------------------

#' Build a species set by IDs, tags, and/or taxonomy filters
#' @param ids character vector of g:Profiler organism IDs (optional)
#' @param tags character vector of tag names, e.g. c("AnimalModels") (optional)
#' @param taxon named list of rank=value pairs, e.g. list(order="Perciformes")
#' @param include_related if TRUE and tags are given, union tag set with taxon/ids
#' @export
echogo_select_species <- function(ids = NULL, tags = NULL, taxon = NULL,
                                  include_related = TRUE, refresh = FALSE) {
  sp <- echogo_species_table(refresh = refresh)
  pick <- rep(TRUE, nrow(sp))

  if (!is.null(ids) && length(ids)) {
    pick <- pick & (sp$organism %in% ids)
  }
  if (!is.null(tags) && length(tags)) {
    has_tag <- vapply(sp$tags, function(x) any(x %in% tags), logical(1))
    pick <- if (isTRUE(include_related)) (pick | has_tag) else (pick & has_tag)
  }
  if (!is.null(taxon) && length(taxon)) {
    for (rk in names(taxon)) {
      rk <- tolower(rk)
      if (!rk %in% c("superkingdom","kingdom","phylum","class","order","family","genus")) next
      target <- tolower(as.character(taxon[[rk]]))
      pick <- pick & tolower(sp[[rk]]) %in% target
    }
  }
  unique(sp$organism[pick])
}

#' Resolve species expressions like "tag:AnimalModels OR order:Perciformes"
#' Supports tokens:
#'   - tag:NAME                    (matches list-column 'tags')
#'   - organism:ID or id:ID        (exact g:Profiler organism code)
#'   - superkingdom:/kingdom:/phylum:/class:/order:/family:/genus:VALUE (exact, case-insensitive)
#'
#' Operators:
#'   - OR  (union)
#'   - AND (intersection)
#' If you provide multiple tokens with NO operator, the default is AND (intersection).
#'
#' Returns a character vector of valid organism IDs present in the species table.
#' Never returns NA.
#' @export
echogo_resolve <- function(expr, refresh = FALSE) {
  sp <- echogo_species_table(refresh = refresh)

  # guard: ensure needed cols exist
  rank_cols <- c("superkingdom","kingdom","phylum","class","order","family","genus")
  for (rc in rank_cols) if (!rc %in% names(sp)) sp[[rc]] <- NA_character_
  if (!"tags" %in% names(sp)) sp$tags <- replicate(nrow(sp), character(0), simplify = FALSE)

  # helper: return a set of organism codes for one token
  to_set <- function(tok) {
    tok <- trimws(tok)
    if (!nzchar(tok)) return(character(0))

    if (grepl("^tag:", tok, ignore.case = TRUE)) {
      tg <- sub("^tag:", "", tok, ignore.case = TRUE)
      idx <- vapply(sp$tags, function(x) isTRUE(tg %in% x), logical(1))
      idx[is.na(idx)] <- FALSE
      return(sp$organism[!is.na(sp$organism) & idx])

    } else if (grepl("^[a-z]+:", tok)) {
      k <- tolower(sub(":.*", "", tok))
      v <- sub("^[^:]+:", "", tok)
      if (k %in% c("id","organism")) k <- "organism"

      valid_keys <- c("organism", rank_cols)
      if (!k %in% valid_keys) return(character(0))

      col <- sp[[k]]
      idx <- !is.na(col) & tolower(as.character(col)) == tolower(v)
      return(sp$organism[!is.na(sp$organism) & idx])

    } else {
      # bare organism code
      return(sp$organism[!is.na(sp$organism) & sp$organism == tok])
    }
  }

  # tokenize
  parts <- strsplit(gsub("\\s+", " ", trimws(expr)), "\\s+")[[1]]
  if (!length(parts)) return(character(0))

  sets <- list()
  ops  <- character(0)
  for (w in parts) {
    up <- toupper(w)
    if (up %in% c("OR","AND")) {
      ops <- c(ops, up)
    } else {
      sets <- c(sets, list(unique(to_set(w))))
    }
  }

  if (!length(sets)) return(character(0))

  # combine sets
  if (!length(ops)) {
    # default: AND across all tokens
    out <- sets[[1]]
    if (length(sets) > 1) {
      for (i in 2:length(sets)) out <- intersect(out, sets[[i]])
    }
  } else {
    # left-associative evaluation using provided ops between consecutive sets
    out <- sets[[1]]
    for (i in seq_along(ops)) {
      if (ops[i] == "OR")  out <- union(out,     sets[[i+1]])
      if (ops[i] == "AND") out <- intersect(out, sets[[i+1]])
    }
  }

  # sanitize: drop NA/empty and keep only supported organism IDs
  out <- unique(out)
  out <- out[!is.na(out) & nzchar(out)]
  out <- out[out %in% sp$organism]
  sort(out)
}

# ---- Smart lookup (kept; minor polish) --------------------------------------

#' Suggest best-matching g:Profiler IDs for each input (smart ranking)
#' @export
echogo_species_lookup <- function(query,
                                  fields  = c("organism","name","alias","ncbi"),
                                  top_n   = 5,
                                  fuzzy   = TRUE,
                                  refresh = FALSE) {
  stopifnot(is.character(query), length(query) >= 1)
  sp <- echogo_gprofiler_species(refresh = refresh)
  if (!nrow(sp)) return(sp[0, c("organism","name","ncbi","alias")])

  model_ids <- c("hsapiens","mmusculus","rnorvegicus","drerio","dmelanogaster","celegans","scerevisiae","ggallus")
  synonyms <- c(
    "human"="hsapiens","man"="hsapiens","mouse"="mmusculus","mice"="mmusculus",
    "rat"="rnorvegicus","zebrafish"="drerio","zebra fish"="drerio",
    "fruit fly"="dmelanogaster","fly"="dmelanogaster","worm"="celegans",
    "yeast"="scerevisiae","chicken"="ggallus"
  )

  to_char <- function(x) if (is.numeric(x)) format(x, trim = TRUE, scientific = FALSE) else as.character(x %||% "")
  L <- list(
    organism = tolower(sp$organism),
    name     = tolower(sp$name),
    alias    = tolower(to_char(sp$alias)),
    ncbi     = tolower(to_char(sp$ncbi))
  )
  fields <- intersect(fields, names(L)); if (!length(fields)) fields <- names(L)
  ad <- function(q, x) utils::adist(q, x, partial = TRUE, ignore.case = TRUE)
  is_digits <- function(s) grepl("^\\d+$", s)

  rank_one <- function(q) {
    ql <- tolower(trimws(q)); if (!nzchar(ql)) return(NULL)
    syn_id <- unname(synonyms[ql]); syn_ix <- if (!is.na(syn_id)) which(sp$organism == syn_id) else integer()
    exact_id_ix   <- which(L$organism == ql)
    exact_ncbi_ix <- if (is_digits(ql)) which(L$ncbi == ql) else integer()
    exact_name_ix <- which(L$name == ql | L$alias == ql)
    word_pat <- paste0("\\b", gsub("([.^$|()*+?{}\\[\\]\\\\])", "\\\\\\1", ql), "\\b")
    word_ix  <- unique(which(Reduce(`|`, lapply(L[fields], function(col) grepl(word_pat, col)))))
    starts_ix <- unique(which(Reduce(`|`, lapply(L[fields], function(col) startsWith(col, ql)))))
    sub_ix    <- unique(which(Reduce(`|`, lapply(L[fields], function(col) grepl(ql, col, fixed = TRUE)))))
    fuzzy_ix <- integer()
    if (isTRUE(fuzzy) && !is_digits(ql)) {
      fuzzy_ix <- unique(unlist(lapply(L[fields], function(col) {
        which(vapply(col, function(x) length(base::agrep(ql, x, ignore.case = TRUE, max.distance = 0.1)) > 0, logical(1)))
      })))
    }
    tiers <- list(
      `-1_synonym` = syn_ix,
      `0_exact_id_or_ncbi` = unique(c(exact_id_ix, exact_ncbi_ix)),
      `1_exact_name_or_alias` = setdiff(exact_name_ix, c(exact_id_ix, exact_ncbi_ix)),
      `2_whole_word` = setdiff(word_ix, c(exact_id_ix, exact_ncbi_ix, exact_name_ix)),
      `3_starts_with` = setdiff(starts_ix, c(exact_id_ix, exact_ncbi_ix, exact_name_ix, word_ix)),
      `4_substring` = setdiff(sub_ix, c(exact_id_ix, exact_ncbi_ix, exact_name_ix, word_ix, starts_ix)),
      `5_fuzzy` = setdiff(fuzzy_ix, c(exact_id_ix, exact_ncbi_ix, exact_name_ix, word_ix, starts_ix, sub_ix))
    )
    any_ix <- unique(unlist(tiers)); if (!length(any_ix)) return(NULL)
    cand <- sp[any_ix, , drop = FALSE]
    tier_map <- setNames(seq_along(tiers)-2L, names(tiers))
    tier_vec <- integer(length(any_ix))
    for (nm in names(tiers)) if (length(tiers[[nm]])) tier_vec[match(tiers[[nm]], any_ix)] <- tier_map[[nm]]
    d_name  <- ad(ql, tolower(to_char(cand$name)))
    d_org   <- ad(ql, tolower(to_char(cand$organism)))
    d_alias <- ad(ql, tolower(to_char(cand$alias)))
    d_ncbi  <- ad(ql, tolower(to_char(cand$ncbi)))
    dist_sc <- pmin(d_name, d_org, d_alias, d_ncbi, na.rm = TRUE)
    model_boost <- ifelse(cand$organism %in% model_ids, 0L, 1L)

    out <- tibble::tibble(
      query = q,
      organism = cand$organism,
      name = cand$name,
      ncbi = cand$ncbi,
      alias = cand$alias,
      tier = tier_vec,
      score = as.numeric(dist_sc),
      model = model_boost
    )
    out <- out[order(out$tier, out$score, out$model, out$organism), , drop = FALSE]
    out$rank <- ave(out$score, out$query, FUN = function(z) rank(z, ties.method = "first"))
    head(out, top_n)
  }

  pieces <- lapply(unique(query), rank_one)
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) return(sp[0, c("organism","name","ncbi","alias")])
  out <- dplyr::bind_rows(pieces)
  rownames(out) <- NULL
  out[, c("query","organism","name","ncbi","alias","score","rank")]
}

# ---- Front-door UX wrapper --------------------------------------------------

#' Compact species / OrgDb guide for EchoGO
#'
#' This is a convenience wrapper around [echogo_list_species()] that gives users
#' a quick overview of supported g:Profiler organism IDs and their taxonomy.
#' In interactive sessions with the DT package installed, it will show an
#' interactive table; otherwise it prints a console preview.
#'
#' @param view logical; passed to [echogo_list_species()], default = interactive().
#' @param refresh logical; passed to [echogo_list_species()], default = FALSE.
#' @param n integer; number of rows to show in console preview.
#'
#' @return Invisibly returns the species tibble.
#' @export
echogo_species_guide <- function(view = interactive(), refresh = FALSE, n = 30) {
  echogo_list_species(view = view, refresh = refresh, n = n)
}


test_that("target-supported Key Findings are target-first while hypotheses stay recurrence-first", {
  ev <- tibble::tibble(
    term_id = c("GO:0000001", "GO:0000002", "GO:0000003", "GO:0000004"),
    term_name = c("target strong", "target recurrent", "hyp recurrent", "hyp narrow"),
    ontology = "BP",
    evidence_profile = c("TARGET_PLUS_CONTEXT", "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT", "ALTERNATIVE_CONTEXT"),
    target_goseq_adjusted_p = c(0.001, 0.02, NA, NA),
    alternative_context_support_n = c(2L, 8L, 7L, 2L),
    alternative_queried_context_n = 8L,
    alternative_context_support_fraction = c(2/8, 1, 7/8, 2/8),
    contributing_genes = c("a;b", "c;d", "e;f", "g;h"),
    display_order = 1:4
  )

  out <- EchoGO:::.echogo_key_findings_data(
    ev,
    ontology = "BP",
    max_target_recovered = 2L,
    max_target_only = 0L,
    max_hypotheses = 2L
  )

  plus <- out[out$evidence_profile == "TARGET_PLUS_CONTEXT", , drop = FALSE]
  hyp <- out[out$evidence_profile == "ALTERNATIVE_CONTEXT", , drop = FALSE]
  expect_identical(plus$term_id, c("GO:0000001", "GO:0000002"))
  expect_identical(hyp$term_id, c("GO:0000003", "GO:0000004"))
})

test_that("Evidence Landscape is capped at ten terms per profile", {
  ev <- tibble::tibble(
    term_id = sprintf("GO:%07d", 1:36),
    term_name = paste("term", 1:36),
    ontology = "BP",
    evidence_profile = rep(c("TARGET_ONLY", "TARGET_PLUS_CONTEXT", "ALTERNATIVE_CONTEXT"), each = 12),
    display_order = 1:36
  )
  out <- EchoGO:::.echogo_landscape_terms(ev, "BP", max_terms_per_profile = 12L)
  counts <- table(out$evidence_profile)
  expect_true(all(counts <= 10L))
  expect_equal(nrow(out), 30L)
})

test_that("PNG integrity helper rejects missing and corrupt assets", {
  td <- tempfile("echogo-png-")
  dir.create(td)
  bad <- file.path(td, "bad.png")
  writeLines("not a png", bad)
  expect_false(EchoGO:::.echogo_png_is_renderable(file.path(td, "missing.png")))
  expect_false(EchoGO:::.echogo_png_is_renderable(bad))

  good <- file.path(td, "good.png")
  grDevices::png(good, width = 200, height = 200)
  graphics::plot.new()
  grDevices::dev.off()
  expect_true(EchoGO:::.echogo_png_is_renderable(good))
})

test_that("final report polish repairs copy and adds print expansion", {
  td <- tempfile("echogo-html-")
  dir.create(td)
  html <- file.path(td, "report.html")
  writeLines(c(
    "<html><body>",
    "<p>whatthe Key FindingsMap experimentto role:target reference biologicalthemes</p>",
    "<h3>GO:0004656 (GO:0004656): example</h3>",
    "<details><summary>reference</summary><div>plot</div></details>",
    "</body></html>"
  ), html)
  expect_true(EchoGO:::.echogo_final_report_html_polish(html))
  txt <- paste(readLines(html), collapse = "\n")
  expect_match(txt, "what the")
  expect_match(txt, "Key Findings Map")
  expect_match(txt, "experiment to")
  expect_match(txt, "role: target reference")
  expect_match(txt, "biological themes")
  expect_false(grepl("GO:0004656 \\(GO:0004656\\)", txt))
  expect_match(txt, "beforeprint")
  expect_match(txt, "d.open=true", fixed = TRUE)
})

test_that("report asset validator fails on missing or corrupt relative image", {
  td <- tempfile("echogo-assets-")
  dir.create(td)
  dir.create(file.path(td, "assets"))
  html <- file.path(td, "report.html")
  writeLines("<html><body><img src='assets/missing.png'/></body></html>", html)
  expect_error(
    EchoGO:::.echogo_validate_report_image_assets(html, td),
    "missing or invalid image assets"
  )

  good <- file.path(td, "assets", "good.png")
  grDevices::png(good, width = 200, height = 200)
  graphics::plot.new()
  grDevices::dev.off()
  writeLines("<html><body><img src='assets/good.png'/></body></html>", html)
  expect_silent(EchoGO:::.echogo_validate_report_image_assets(html, td))
})

test_that("RRvGO cluster tables can regenerate a valid treemap PNG", {
  skip_if_not_installed("treemap")
  skip_if_not_installed("scales")
  td <- tempfile("echogo-rrvgo-")
  dir.create(td)
  clusters <- file.path(td, "rrvgo_BP_clusters.csv")
  utils::write.csv(
    data.frame(
      go = c("GO:0000001", "GO:0000002", "GO:0000003"),
      term = c("alpha process", "beta process", "gamma process"),
      parentTerm = c("theme A", "theme A", "theme B"),
      size = c(5, 3, 4),
      score = c(3, 2, 1)
    ),
    clusters,
    row.names = FALSE
  )
  png <- file.path(td, "rrvgo_BP_treemap.png")
  expect_true(EchoGO:::.echogo_rrvgo_treemap_png_from_clusters(clusters, png))
  expect_true(EchoGO:::.echogo_png_is_renderable(png))
})

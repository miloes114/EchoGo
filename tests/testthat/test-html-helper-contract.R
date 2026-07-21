test_that("missing report assets return well-formed diagnostic HTML", {
  html <- embed_html_toggle_external_plus(
    title = "Interactive network",
    path_abs = NA_character_,
    report_dir = tempdir(),
    base_dir = tempdir()
  )
  expect_match(html, "Missing HTML", fixed = TRUE)
  expect_equal(lengths(regmatches(html, gregexpr("<em>", html, fixed = TRUE))), 1L)
  expect_equal(lengths(regmatches(html, gregexpr("</em>", html, fixed = TRUE))), 1L)
})

test_that("HTML widgets retain their companion dependency directory", {
  root <- file.path(tempdir(), paste0("echogo-widget-", Sys.getpid()))
  unlink(root, recursive = TRUE, force = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  widget_dir <- file.path(root, "networks", "with_bg")
  dependency_dir <- file.path(widget_dir, "network_BP_files", "vis")
  report_dir <- file.path(root, "report")
  dir.create(dependency_dir, recursive = TRUE)

  widget <- file.path(widget_dir, "network_BP.html")
  script <- file.path(dependency_dir, "vis-network.min.js")
  stylesheet <- file.path(dependency_dir, "vis-network.min.css")
  pdf <- file.path(widget_dir, "network_BP.pdf")
  writeLines(
    c(
      '<link href="network_BP_files/vis/vis-network.min.css" rel="stylesheet">',
      '<script src="network_BP_files/vis/vis-network.min.js"></script>'
    ),
    widget
  )
  writeLines("window.visNetworkLoaded = true;", script)
  writeLines("#network { width: 100%; }", stylesheet)
  writeBin(charToRaw("%PDF-1.4\n%%EOF\n"), pdf)

  html <- embed_html_toggle_external_plus(
    title = "Interactive network",
    path_abs = widget,
    report_dir = report_dir,
    base_dir = root
  )

  staged <- file.path(report_dir, "assets", "networks", "with_bg")
  expect_true(file.exists(file.path(staged, "network_BP.html")))
  expect_true(file.exists(file.path(staged, "network_BP_files", "vis", "vis-network.min.js")))
  expect_true(file.exists(file.path(staged, "network_BP_files", "vis", "vis-network.min.css")))
  expect_identical(readLines(file.path(staged, "network_BP_files", "vis", "vis-network.min.js")), readLines(script))
  expect_true(file.exists(file.path(staged, "network_BP.pdf")))
  expect_match(html, "assets/networks/with_bg/network_BP.html", fixed = TRUE)
  expect_match(html, "Open the static PDF version", fixed = TRUE)
})

test_that("external HTML widget dependencies are staged together", {
  root <- file.path(tempdir(), paste0("echogo-report-", Sys.getpid()))
  external <- file.path(tempdir(), paste0("echogo-external-", Sys.getpid()))
  unlink(c(root, external), recursive = TRUE, force = TRUE)
  on.exit(unlink(c(root, external), recursive = TRUE, force = TRUE), add = TRUE)

  dep_dir <- file.path(external, "widget_files", "htmlwidgets")
  dir.create(dep_dir, recursive = TRUE)
  widget <- file.path(external, "widget.html")
  dependency <- file.path(dep_dir, "htmlwidgets.js")
  writeLines('<script src="widget_files/htmlwidgets/htmlwidgets.js"></script>', widget)
  writeLines("window.HTMLWidgets = {};", dependency)

  html <- embed_html_toggle_external_plus(
    title = "External widget",
    path_abs = widget,
    report_dir = file.path(root, "report"),
    base_dir = root
  )

  staged <- file.path(root, "report", "assets", "external")
  expect_true(file.exists(file.path(staged, "widget.html")))
  expect_true(file.exists(file.path(staged, "widget_files", "htmlwidgets", "htmlwidgets.js")))
  expect_match(html, "assets/external/widget.html", fixed = TRUE)
})

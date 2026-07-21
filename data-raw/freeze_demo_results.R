# Freeze corrected demo results into a new versioned directory.

if (!requireNamespace("rappdirs", quietly = TRUE)) {
  stop("Please install 'rappdirs' to locate the demo directory.")
}

demo_root <- getOption("EchoGO.demo_root", rappdirs::user_data_dir("EchoGO"))
src <- getOption(
  "EchoGO.demo_results_source",
  file.path(demo_root, "echogo_demo", "results")
)
dst <- file.path("inst", "extdata", "echogo_demo_results_v0.1.3")

if (!dir.exists(src)) {
  stop("Run echogo_quickstart(run_demo = TRUE) once to produce demo results at:\n", src)
}

src_files <- list.files(src, all.files = TRUE, full.names = TRUE, recursive = TRUE, no.. = TRUE)
src_files <- src_files[basename(src_files) != "__report_render.log"]
if (length(src_files) == 0) {
  stop("Demo results folder exists but is empty:\n", src, "\nRe-run echogo_quickstart(run_demo = TRUE).")
}

# Clean only the new v0.1.3 destination. The pre-correction snapshot remains
# available under inst/extdata/echogo_demo_results.
dir.create(dst, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(dst, all.files = TRUE, full.names = TRUE, no.. = TRUE),
       recursive = TRUE, force = TRUE)

# Compute relative paths in a platform-safe way
rel <- substring(src_files, nchar(normalizePath(src, winslash = "/", mustWork = TRUE)) + 2L)
rel <- gsub("/", .Platform$file.sep, rel, fixed = TRUE)   # normalize separators for Windows
dest_files <- file.path(dst, rel)

# Create destination directories (one by one; dir.create can't take a vector)
dest_dirs <- unique(dirname(dest_files))
for (d in dest_dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Copy files
ok <- file.copy(from = src_files, to = dest_files, overwrite = TRUE)
if (!all(ok)) {
  warning("Some demo results could not be copied. First few failures:\n",
          paste(head(src_files[!ok], 10), collapse = "\n"))
}

# Remove machine-specific paths from indexes copied into the package.
index_files <- list.files(
  dst,
  pattern = "^__file_index\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
for (index_file in index_files) {
  index <- utils::read.csv(index_file, check.names = FALSE)
  if (all(c("rel_path", "full_path") %in% names(index))) {
    index$full_path <- paste0("<output_dir>/", index$rel_path)
    utils::write.csv(index, index_file, row.names = FALSE, quote = FALSE)
  }
}

# Replace local roots in text artifacts while retaining useful relative paths.
text_files <- list.files(
  dst,
  pattern = "\\.(csv|html|log|md|Rmd|txt)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)
local_roots <- unique(c(
  normalizePath(src, winslash = "/", mustWork = TRUE),
  normalizePath(path.expand("~"), winslash = "/", mustWork = TRUE)
))
replacements <- c("<output_dir>", "<user_home>")
for (text_file in text_files) {
  text <- readLines(text_file, warn = FALSE, encoding = "UTF-8")
  for (i in seq_along(local_roots)) {
    text <- gsub(local_roots[[i]], replacements[[i]], text, fixed = TRUE)
    text <- gsub(
      chartr("/", "\\", local_roots[[i]]),
      replacements[[i]],
      text,
      fixed = TRUE
    )
  }
  writeLines(text, text_file, useBytes = TRUE)
}

message("[EchoGO] Frozen demo results at: ", normalizePath(dst, winslash = "/"))


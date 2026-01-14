# data-raw/freeze_demo_results.R
# Freeze a local demo run into inst/extdata/echogo_demo_results for offline vignettes.

if (!requireNamespace("rappdirs", quietly = TRUE)) {
  stop("Please install 'rappdirs' to locate the demo directory.")
}

demo_root <- getOption("EchoGO.demo_root", rappdirs::user_data_dir("EchoGO"))
src <- file.path(demo_root, "echogo_demo", "results")
dst <- file.path("inst", "extdata", "echogo_demo_results")

if (!dir.exists(src)) {
  stop("Run echogo_quickstart(run_demo = TRUE) once to produce demo results at:\n", src)
}

src_files <- list.files(src, all.files = TRUE, full.names = TRUE, recursive = TRUE, no.. = TRUE)
if (length(src_files) == 0) {
  stop("Demo results folder exists but is empty:\n", src, "\nRe-run echogo_quickstart(run_demo = TRUE).")
}

# Clean destination
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

message("[EchoGO] Frozen demo results at: ", normalizePath(dst, winslash = "/"))


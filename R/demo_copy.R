# Internal filesystem helper for packaged demo material -----------------------
#
# `file.copy()` behaves differently across platforms when a single call mixes
# ordinary files and directories in its `from` vector.  In particular, the
# installed Windows package quickstart can fail when the packaged demo contains
# both top-level input files and the nested demo cache directory.
#
# Copy directory trees as matched file-to-file paths instead.  This keeps the
# packaged cache byte-identical, preserves nested paths, and avoids depending on
# recursive directory-copy semantics.

.echogo_copy_tree <- function(from, to, overwrite = TRUE) {
  if (length(from) != 1L || is.na(from) || !nzchar(from) || !dir.exists(from)) {
    stop("EchoGO copy source must be one existing directory: ", from, call. = FALSE)
  }
  if (length(to) != 1L || is.na(to) || !nzchar(to)) {
    stop("EchoGO copy destination must be one directory path.", call. = FALSE)
  }

  dir.create(to, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(to)) {
    stop("Could not create EchoGO copy destination: ", to, call. = FALSE)
  }

  # Recreate directories first, including empty directories.  `list.dirs()`
  # may represent the source root as either "." or an empty string depending on
  # platform/R version, so both are removed explicitly.
  rel_dirs <- list.dirs(from, full.names = FALSE, recursive = TRUE)
  rel_dirs <- rel_dirs[!is.na(rel_dirs) & nzchar(rel_dirs) & rel_dirs != "."]
  if (length(rel_dirs)) {
    for (rel_dir in rel_dirs) {
      dir.create(file.path(to, rel_dir), recursive = TRUE, showWarnings = FALSE)
    }
  }

  rel_files <- list.files(
    from,
    recursive = TRUE,
    full.names = FALSE,
    all.files = TRUE,
    no.. = TRUE,
    include.dirs = FALSE
  )
  if (!length(rel_files)) return(invisible(character()))

  source_files <- file.path(from, rel_files)
  target_files <- file.path(to, rel_files)

  # Defensive creation of parent directories also covers unusual directory
  # enumeration behaviour and makes the helper safe for deeply nested caches.
  parent_dirs <- unique(dirname(target_files))
  for (parent_dir in parent_dirs) {
    dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)
  }

  copied <- file.copy(
    from = source_files,
    to = target_files,
    overwrite = isTRUE(overwrite),
    copy.mode = TRUE,
    copy.date = TRUE
  )

  if (length(copied) != length(source_files) || !all(copied)) {
    failed <- rel_files[seq_along(copied) <= length(copied) & !copied]
    if (!length(failed) && length(copied) != length(source_files)) {
      failed <- rel_files
    }
    stop(
      "Could not copy the packaged EchoGO demo tree. Failed paths: ",
      paste(failed, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(target_files)
}

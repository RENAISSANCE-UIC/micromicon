#' File Utility Functions
#'
#' @description
#' Low-level file operation utilities used across the framework.
#'
#' @keywords internal


#' Check if file exists
#'
#' @description
#' Checks if a file exists and optionally throws an error if not found.
#'
#' @param path Character string file path
#' @param error_on_missing Logical; throw error if file missing (default TRUE)
#' @return Logical indicating if file exists
#' @keywords internal
check_file_exists <- function(path, error_on_missing = TRUE) {
  exists <- file.exists(path)

  if (!exists && error_on_missing) {
    cli::cli_abort("File not found: {path}")
  }

  exists
}


#' Normalize file path
#'
#' @description
#' Normalizes a file path to absolute path with forward slashes.
#'
#' @param path Character string file path
#' @return Normalized absolute path
#' @keywords internal
normalize_path <- function(path) {
  if (length(path) == 0 || !nzchar(path)) {
    return(path)
  }

  # Expand ~ and normalize
  path <- path.expand(path)
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  path
}


#' Create temporary file with extension
#'
#' @description
#' Creates a temporary file with specified extension.
#'
#' @param fileext Character string file extension (e.g., ".fasta")
#' @param pattern Character string pattern for tempfile name
#' @return Path to temporary file
#' @keywords internal
create_temp_file <- function(fileext = "", pattern = "micromicon_") {
  tempfile(pattern = pattern, fileext = fileext)
}


#' Read lines from file safely
#'
#' @description
#' Reads lines from a file with error handling.
#'
#' @param path Character string file path
#' @param warn Logical; show warnings (default FALSE)
#' @return Character vector of lines
#' @keywords internal
read_lines_safe <- function(path, warn = FALSE) {
  check_file_exists(path, error_on_missing = TRUE)

  tryCatch(
    readLines(path, warn = warn),
    error = function(e) {
      cli::cli_abort(c(
        "Failed to read file: {path}",
        "x" = e$message
      ))
    }
  )
}


#' Write lines to file safely
#'
#' @description
#' Writes lines to a file with error handling.
#'
#' @param text Character vector of lines
#' @param path Character string file path
#' @return Invisibly returns the path
#' @keywords internal
write_lines_safe <- function(text, path) {
  tryCatch(
    {
      writeLines(text, path)
      invisible(path)
    },
    error = function(e) {
      cli::cli_abort(c(
        "Failed to write file: {path}",
        "x" = e$message
      ))
    }
  )
}


#' Get file extension
#'
#' @description
#' Extracts file extension from path (including the dot).
#'
#' @param path Character string file path
#' @return Character string extension (e.g., ".fasta")
#' @keywords internal
get_file_ext <- function(path) {
  parts <- strsplit(basename(path), ".", fixed = TRUE)[[1]]

  if (length(parts) <= 1) {
    return("")
  }

  paste0(".", parts[length(parts)])
}


#' Detect file format from extension
#'
#' @description
#' Detects file format (genbank, fasta, gff3) from file extension.
#'
#' @param path Character string file path
#' @return Character string format name or NA if unknown
#' @keywords internal
detect_format_from_extension <- function(path) {
  ext <- tolower(get_file_ext(path))

  if (ext %in% c(".gbk", ".gb", ".genbank")) {
    return("genbank")
  } else if (ext %in% c(".fasta", ".fa", ".fna", ".faa")) {
    return("fasta")
  } else if (ext %in% c(".gff", ".gff3")) {
    return("gff3")
  } else {
    return(NA_character_)
  }
}

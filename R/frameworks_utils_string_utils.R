#' String Utility Functions
#'
#' @description
#' Low-level string manipulation utilities used across the framework.
#'
#' @keywords internal


#' Infix operator for default values
#'
#' @description
#' A tiny infix helper to mirror rlang's %||%, avoiding hard dependency.
#' Returns y if x is NULL, otherwise returns x.
#'
#' @param x First value to check
#' @param y Default value if x is NULL
#' @return x if x is not NULL, otherwise y
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


#' Scrub common sequence name prefixes
#'
#' @description
#' Removes common prefixes from sequence names like "lcl|" and "chr".
#' Extend as your corpus dictates.
#'
#' @param x Character vector of sequence names
#' @return Character vector with prefixes removed
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' .scrub_prefixes(c("lcl|seq1", "chr1", "CHR2"))
#' # Returns: c("seq1", "1", "2")
#' }
.scrub_prefixes <- function(x) {
  x <- as.character(x)
  x <- sub("^lcl\\|", "", x, perl = TRUE)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x
}


#' Trim whitespace from strings
#'
#' @description
#' Removes leading and trailing whitespace from character vectors.
#'
#' @param x Character vector
#' @return Character vector with whitespace trimmed
#' @keywords internal
trim <- function(x) {
  sub("\\s+$", "", sub("^\\s+", "", x))
}


#' Wrap sequence string at specified width
#'
#' @description
#' Splits a long sequence string into lines of specified width.
#'
#' @param seq Character string (DNA/protein sequence)
#' @param width Integer; line width (default 80)
#' @return Character vector of wrapped lines
#' @keywords internal
wrap_sequence <- function(seq, width = 80) {
  if (nchar(seq) <= width) {
    return(seq)
  }

  positions <- seq(1, nchar(seq), by = width)
  vapply(positions, function(start) {
    end <- min(start + width - 1, nchar(seq))
    substr(seq, start, end)
  }, character(1))
}

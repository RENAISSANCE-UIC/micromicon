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

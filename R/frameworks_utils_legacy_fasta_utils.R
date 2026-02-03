#' Legacy FASTA Utility Functions
#'
#' @description
#' Legacy functions for FASTA file operations. These are maintained for
#' backward compatibility with the old API but are not used in the new
#' Clean Architecture layers.
#'
#' @keywords internal


#' Get FASTA index information (LEGACY)
#'
#' @description
#' Return FASTA names and lengths from either DNAStringSet or FaFile.
#' This is a legacy function; new code should use bioconductor_gateway methods.
#'
#' @param fa_obj DNAStringSet or FaFile object
#' @return List with 'names' and 'lengths' components
#' @keywords internal
.get_fasta_index <- function(fa_obj) {
  # Require string utils for %||%
  `%||%` <- function(x, y) if (is.null(x)) y else x

  if (inherits(fa_obj, "DNAStringSet")) {
    nm <- names(fa_obj) %||% character(length(fa_obj))
    len <- Biostrings::width(fa_obj)
    return(list(names = nm, lengths = setNames(as.integer(len), nm)))
  }

  if (inherits(fa_obj, "FaFile")) {
    if (!requireNamespace("Rsamtools", quietly = TRUE)) {
      stop("Rsamtools is required to inspect FaFile indices.")
    }
    idx <- Rsamtools::scanFaIndex(fa_obj)
    nm <- names(idx)
    len <- vapply(idx, function(rec) rec$seqlength, integer(1))
    return(list(names = nm, lengths = setNames(as.integer(len), nm)))
  }

  # Fallback if you only have headers
  list(names = character(), lengths = integer())
}

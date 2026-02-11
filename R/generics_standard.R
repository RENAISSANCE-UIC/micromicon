#' Summary of genome_entity
#'
#' @param object A genome_entity object
#' @param ... Additional arguments (ignored)
#' @export
summary.genome_entity <- function(object, ...) {
  cat("genome_entity Summary\n")
  cat("=====================\n\n")

  # Sequences
  n_seqs <- length(object$sequences$dna_raw)
  total_bp <- if (n_seqs > 0) sum(nchar(object$sequences$dna_raw)) else 0

  cat("Sequences:", n_seqs, "sequence")
  if (n_seqs != 1) cat("s")
  cat(", ", format(total_bp, big.mark = ","), " bp total\n", sep = "")

  # Features
  n_feat <- nrow(object$features)
  cat("Features:", n_feat, "\n")

  # Records
  n_records <- nrow(object$metadata)
  cat("Records:", n_records, "\n\n")

  # Feature type breakdown
  if (n_feat > 0 && "type" %in% names(object$features)) {
    cat("Feature Types:\n\n")
    print(table(object$features$type))
    cat("\n")
  }

  # Sequence details
  if (n_seqs > 0) {
    cat("Sequences:\n")
    seq_df <- data.frame(
      seqname = names(object$sequences$dna_raw),
      length_bp = nchar(object$sequences$dna_raw)
    )
    print(seq_df, row.names = FALSE)
  }

  invisible(object)
}

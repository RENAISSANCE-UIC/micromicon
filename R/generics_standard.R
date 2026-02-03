#' Print genome_entity
#'
#' @param x A genome_entity object
#' @param ... Additional arguments (ignored)
#' @export
print.genome_entity <- function(x, ...) {
  cat("genome_entity\n")
  cat("=============\n\n")

  # Sequences
  n_seqs <- length(x$sequences$dna_raw)
  if (n_seqs > 0) {
    total_bp <- sum(nchar(x$sequences$dna_raw))
    cat("Sequences:", n_seqs, "sequence")
    if (n_seqs != 1) cat("s")
    cat(", ", format(total_bp, big.mark = ","), " bp total\n", sep = "")

    if (n_seqs <= 5) {
      seq_names <- names(x$sequences$dna_raw)
      seq_lens <- nchar(x$sequences$dna_raw)
      for (i in seq_along(seq_names)) {
        cat("  -", seq_names[i], ":", format(seq_lens[i], big.mark = ","), "bp\n")
      }
    } else {
      seq_names <- names(x$sequences$dna_raw)
      seq_lens <- nchar(x$sequences$dna_raw)
      for (i in 1:3) {
        cat("  -", seq_names[i], ":", format(seq_lens[i], big.mark = ","), "bp\n")
      }
      cat("  ... and", n_seqs - 4, "more\n")
      i <- n_seqs
      cat("  -", seq_names[i], ":", format(seq_lens[i], big.mark = ","), "bp\n")
    }
  } else {
    cat("Sequences: none\n")
  }

  cat("\n")

  # Features
  n_feat <- nrow(x$features)
  if (n_feat > 0) {
    cat("Features:", n_feat, "feature")
    if (n_feat != 1) cat("s")
    cat("\n")

    if ("type" %in% names(x$features)) {
      type_counts <- table(x$features$type)
      top_types <- head(sort(type_counts, decreasing = TRUE), 5)
      for (i in seq_along(top_types)) {
        cat("  -", names(top_types)[i], ":", top_types[i], "\n")
      }
      if (length(type_counts) > 5) {
        cat("  ... and", length(type_counts) - 5, "more type")
        if (length(type_counts) - 5 != 1) cat("s")
        cat("\n")
      }
    }
  } else {
    cat("Features: none\n")
  }

  cat("\n")

  # Records
  n_records <- nrow(x$metadata)
  if (n_records > 0) {
    cat("Records:", n_records, "record")
    if (n_records != 1) cat("s")
    cat("\n")
  }

  invisible(x)
}


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

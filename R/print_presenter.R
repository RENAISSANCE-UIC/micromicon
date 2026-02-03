#' Print Method for genome_entity
#'
#' @description
#' Presenter for formatting genome_entity objects for console display.
#' This is part of the Interface Adapters layer, responsible for output formatting.
#'
#' @param x A genome_entity object
#' @param ... Additional arguments (ignored)
#' @return The object (invisibly)
#' @export
print.genome_entity <- function(x, ...) {
  cat("genome_entity\n")
  cat("=============\n\n")

  # Source information (if available in metadata)
  if ("source" %in% names(x$metadata)) {
    if (length(x$metadata$source) > 0 && nchar(x$metadata$source[1]) > 0) {
      cat("Source:", x$metadata$source[1], "\n")
    }
  }

  if ("import_date" %in% names(x$metadata)) {
    if (length(x$metadata$import_date) > 0) {
      cat("Import date:", format(x$metadata$import_date[1], '%Y-%m-%d %H:%M:%S'), "\n")
    }
  }

  cat("\n")

  # Sequences
  n_seqs <- length(x$sequences$dna_raw)
  if (n_seqs > 0) {
    total_bp <- sum(nchar(x$sequences$dna_raw))
    cat("Sequences:", n_seqs, "sequence")
    if (n_seqs != 1) cat("s")
    cat(", ", format(total_bp, big.mark = ","), " bp total\n", sep = "")

    # Show individual sequences if <= 5
    if (n_seqs <= 5) {
      seq_names <- names(x$sequences$dna_raw)
      seq_lens <- nchar(x$sequences$dna_raw)
      for (i in seq_along(seq_names)) {
        cat("  -", seq_names[i], ":", format(seq_lens[i], big.mark = ","), "bp\n")
      }
    } else {
      # Show first 3 and last 1
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

    # Show feature type breakdown if type column exists
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

  # Metadata records
  n_records <- nrow(x$metadata)
  if (n_records > 0) {
    cat("Records:", n_records, "record")
    if (n_records != 1) cat("s")
    cat("\n")
  }

  invisible(x)
}

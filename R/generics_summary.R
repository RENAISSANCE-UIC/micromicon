#' Summary Method for genome_entity
#'
#' @description
#' Presenter for generating comprehensive summaries of genome_entity objects.
#' This is part of the Interface Adapters layer, responsible for output formatting.
#'
#' @param object A genome_entity object
#' @param ... Additional arguments (ignored)
#' @return A list with summary statistics (invisibly)
#' @export
summary.genome_entity <- function(object, ...) {
  # Validate object first
  validate_genome_entity(object)

  # Build summary statistics
  summary_list <- list(
    n_sequences = length(object$sequences$dna_raw),
    total_bp = sum(nchar(object$sequences$dna_raw)),
    n_features = nrow(object$features),
    n_records = nrow(object$metadata),
    seqnames = object$indices$seqnames,
    source = NA_character_,
    import_date = NA
  )

  # Add source info if available
  if ("source" %in% names(object$metadata) && length(object$metadata$source) > 0) {
    summary_list$source <- object$metadata$source[1]
  }

  if ("import_date" %in% names(object$metadata) && length(object$metadata$import_date) > 0) {
    summary_list$import_date <- object$metadata$import_date[1]
  }

  # Add feature type breakdown if available
  if (nrow(object$features) > 0 && "type" %in% names(object$features)) {
    summary_list$feature_types <- table(object$features$type)
  }

  # Add sequence-level info
  if (length(object$sequences$dna_raw) > 0) {
    summary_list$sequences_info <- data.frame(
      seqname = names(object$sequences$dna_raw),
      length_bp = nchar(object$sequences$dna_raw),
      stringsAsFactors = FALSE
    )
  }

  # Print formatted summary
  cat("genome_entity Summary\n")
  cat("=====================\n\n")

  # Basic info
  if (!is.null(summary_list$source)) {
    cat("Source:", summary_list$source, "\n")
  }
  cat("Sequences:", summary_list$n_sequences, "sequence")
  if (summary_list$n_sequences != 1) cat("s")
  cat(", ", format(summary_list$total_bp, big.mark = ","), " bp total\n", sep = "")
  cat("Features:", summary_list$n_features, "\n")
  cat("Records:", summary_list$n_records, "\n")

  # Feature types
  if (!is.null(summary_list$feature_types)) {
    cat("\nFeature Types:\n")
    print(summary_list$feature_types)
    cat("\n")
  }

  # Sequence info
  if (!is.null(summary_list$sequences_info) && nrow(summary_list$sequences_info) > 0) {
    cat("Sequences:\n")
    print(summary_list$sequences_info, row.names = FALSE)
    cat("\n")
  }

  invisible(summary_list)
}

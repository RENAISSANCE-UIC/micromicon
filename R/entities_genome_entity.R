#' Create a new genome_entity object
#'
#' @description
#' Internal constructor for genome_entity S3 class. This is the unified
#' internal representation for genomic data from both GenBank and GFF3/FASTA sources.
#'
#' This is a pure domain entity with no external dependencies. It represents the
#' core business concept of a "genome" independent of any particular file format
#' or framework.
#'
#' @param sequences List with components:
#'   - dna_raw (character): Named character vector of DNA sequences
#'   - dna_bio (DNAStringSet or NULL): Biostrings representation (optional)
#'   - indexed_fa (FaFile or NULL): Indexed FASTA file handle (optional)
#' @param features data.frame with genomic features (seqname, start, end, strand, type, ...)
#' @param metadata data.frame with sequence metadata (seqname, length, topology, molecule_type, ...)
#' @param indices List with components:
#'   - seqnames (character): Vector of sequence names
#'   - locus_tag_index (named integer): Map locus_tag to feature row
#'   - gene_index (named integer): Map gene name to feature row
#'
#' @return A genome_entity object
#' @keywords internal
#' @export
new_genome_entity <- function(sequences_list = list(
                                 dna_raw = character(),
                                 dna_bio = NULL,
                                 indexed_fa = NULL
                               ),
                               features_df = data.frame(),
                               metadata_df = data.frame(),
                               indices_list = list(
                                 seqnames = character(),
                                 locus_tag_index = integer(),
                                 gene_index = list()
                               )) {
  # Normalize inputs to expected structure
  if (is.null(sequences_list$dna_bio)) sequences_list$dna_bio <- NULL
  if (is.null(sequences_list$indexed_fa)) sequences_list$indexed_fa <- NULL

  structure(
    list(
      sequences = sequences_list,
      features = features_df,
      metadata = metadata_df,
      indices = indices_list
    ),
    class = "genome_entity"
  )
}

#' Validate a genome_entity object
#'
#' @description
#' Validates that a genome_entity object has the correct structure and
#' satisfies domain invariants. This is a pure validation function with
#' no side effects except stopping on error.
#'
#' Domain rules enforced:
#' - Required components present
#' - Data types correct
#' - Feature coordinates valid (within sequence bounds)
#' - Indices consistent with data
#'
#' @param x A genome_entity object
#' @return TRUE invisibly if valid, stops with error if invalid
#' @keywords internal
#' @export
validate_genome_entity <- function(x) {
  if (!inherits(x, "genome_entity")) {
    cli::cli_abort("Object is not a genome_entity")
  }

  # Check required top-level components
  required_components <- c("sequences", "features", "metadata", "indices")
  missing <- setdiff(required_components, names(x))
  if (length(missing) > 0) {
    cli::cli_abort("Missing required components: {paste(missing, collapse = ', ')}")
  }

  # Validate sequences
  if (!is.list(x$sequences)) {
    cli::cli_abort("sequences must be a list")
  }
  if (!"dna_raw" %in% names(x$sequences)) {
    cli::cli_abort("sequences must have 'dna_raw' component")
  }
  if (!is.character(x$sequences$dna_raw) && !is.null(x$sequences$dna_raw)) {
    cli::cli_abort("sequences$dna_raw must be a character vector or NULL")
  }

  # Validate features
  if (!is.data.frame(x$features)) {
    cli::cli_abort("features must be a data.frame")
  }

  # Validate metadata
  if (!is.data.frame(x$metadata)) {
    cli::cli_abort("metadata must be a data.frame")
  }

  # Validate indices
  if (!is.list(x$indices)) {
    cli::cli_abort("indices must be a list")
  }
  if (!"seqnames" %in% names(x$indices)) {
    cli::cli_abort("indices must have 'seqnames' component")
  }

  # Domain validation: Check feature coordinates if we have features and sequences
  if (nrow(x$features) > 0 && length(x$sequences$dna_raw) > 0) {
    # Get sequence lengths
    seq_lengths <- nchar(x$sequences$dna_raw)
    names(seq_lengths) <- names(x$sequences$dna_raw)

    # Check required feature columns
    required_feat_cols <- c("seqname", "start", "end")
    missing_cols <- setdiff(required_feat_cols, names(x$features))
    if (length(missing_cols) > 0) {
      cli::cli_abort("features missing required columns: {paste(missing_cols, collapse = ', ')}")
    }

    # Validate each feature
    for (i in seq_len(nrow(x$features))) {
      feat <- x$features[i, ]

      # Check seqname exists
      if (!feat$seqname %in% names(seq_lengths)) {
        cli::cli_abort("Feature {i} references unknown sequence: {feat$seqname}")
      }

      # Validate coordinates using domain validator
      tryCatch(
        validate_genomic_coordinates(
          start = feat$start,
          end = feat$end,
          seqname = feat$seqname,
          max_length = seq_lengths[[feat$seqname]]
        ),
        error = function(e) {
          cli::cli_abort("Feature {i} ({feat$seqname}:{feat$start}-{feat$end}): {e$message}")
        }
      )

      # Validate strand if present
      if ("strand" %in% names(feat)) {
        tryCatch(
          validate_strand(feat$strand),
          error = function(e) {
            cli::cli_abort("Feature {i} has invalid strand: {e$message}")
          }
        )
      }
    }
  }

  # Validate indices consistency
  if (length(x$sequences$dna_raw) > 0) {
    seq_names_from_seqs <- names(x$sequences$dna_raw)
    if (!all(x$indices$seqnames %in% seq_names_from_seqs)) {
      cli::cli_abort("indices$seqnames contains names not in sequences$dna_raw")
    }
  }

  invisible(TRUE)
}

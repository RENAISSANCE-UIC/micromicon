#' Extract Sequences by Coordinates (Internal Use Case)
#'
#' @description
#' Internal use case for extracting DNA sequences from specified genomic coordinates.
#' Supports extracting regions from one or more sequences with optional
#' reverse complement for minus strand.
#'
#' This is an internal function called by the controller extract_sequences_by_coords().
#' Users should use the controller function instead.
#'
#' @param entity A genome_entity object
#' @param seqname Character, sequence name(s)
#' @param start Integer, start position(s) (1-based, inclusive)
#' @param end Integer, end position(s) (1-based, inclusive)
#' @param options List of options:
#'   - strand: Character, strand ("+", "-", or NULL for forward strand)
#'   - translate: Logical, translate to protein (default FALSE)
#'   - names: Character vector of names for output sequences (optional)
#'
#' @return Character vector of sequences (named if names provided)
#' @keywords internal
execute_extract_sequences_by_coords <- function(entity, seqname, start, end,
                                                 options = list()) {
  # Validate entity
  validate_genome_entity(entity)

  # Parse options
  strand <- options$strand %||% "+"
  translate <- options$translate %||% FALSE
  names_out <- options$names %||% NULL

  # Vectorize inputs
  n <- max(length(seqname), length(start), length(end))
  seqname <- rep_len(seqname, n)
  start <- rep_len(start, n)
  end <- rep_len(end, n)
  strand <- rep_len(strand, n)

  # Validate coordinates
  for (i in seq_len(n)) {
    if (!seqname[i] %in% names(entity$sequences$dna_raw)) {
      cli::cli_abort("Sequence '{seqname[i]}' not found in genome_entity")
    }

    max_length <- nchar(entity$sequences$dna_raw[[seqname[i]]])
    validate_genomic_coordinates(start[i], end[i], seqname[i], max_length)
  }

  # Extract sequences
  sequences <- character(n)

  for (i in seq_len(n)) {
    full_seq <- entity$sequences$dna_raw[[seqname[i]]]

    # Extract region
    subseq <- substr(full_seq, start[i], end[i])

    # Handle strand
    if (!is.na(strand[i]) && strand[i] == "-") {
      subseq <- reverse_complement(subseq)
    }

    # Translate if requested
    if (translate) {
      subseq <- translate_dna(subseq)
    }

    sequences[i] <- subseq
  }

  # Add names
  if (!is.null(names_out)) {
    names(sequences) <- names_out
  } else {
    # Default names: seqname:start-end
    names(sequences) <- paste0(seqname, ":", start, "-", end)
  }

  sequences
}

# Helper: %||% operator
`%||%` <- function(x, y) if (is.null(x)) y else x

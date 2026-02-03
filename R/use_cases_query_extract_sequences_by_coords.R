#' Extract Sequences by Coordinates
#'
#' @description
#' Use case for extracting DNA sequences from specified genomic coordinates.
#' Supports extracting regions from one or more sequences with optional
#' reverse complement for minus strand.
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
#' @export
#' @examples
#' \dontrun{
#' # Extract single region
#' seq <- execute_extract_sequences_by_coords(genome, "chr1", 1, 100)
#'
#' # Extract multiple regions
#' seqs <- execute_extract_sequences_by_coords(
#'   genome,
#'   seqname = c("chr1", "chr1", "plasmid"),
#'   start = c(1, 500, 1),
#'   end = c(100, 600, 200)
#' )
#'
#' # Extract with reverse complement
#' seq <- execute_extract_sequences_by_coords(
#'   genome,
#'   "chr1", 1, 100,
#'   options = list(strand = "-")
#' )
#' }
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
      stop("Sequence '", seqname[i], "' not found in genome_entity", call. = FALSE)
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

#' Create GenBank Gateway
#'
#' @description
#' Creates a gateway for reading GenBank format files. This gateway implements
#' the repository pattern, abstracting file I/O and parsing details from the
#' business logic layer (use cases).
#'
#' Gateway Contract (for use cases):
#' - read(path): Returns list of records, each with:
#'   - metadata: list with fields (locus, accession, length_bp, mol_type, topology, etc.)
#'   - features: data.frame with columns (type, start, end, strand, gene, locus_tag, etc.)
#'   - sequence: character string
#'
#' @return Gateway object (list with read method)
#' @export
#' @examples
#' \dontrun{
#' gateway <- create_genbank_gateway()
#' records <- gateway$read("data.gbk")
#' }
create_genbank_gateway <- function() {
  list(
    read = function(path) {
      # Validate file exists
      if (!file.exists(path)) {
        cli::cli_abort("GenBank file not found: {path}")
      }

      # Read file
      lines <- readLines(path, warn = FALSE)

      # Normalize tabs to spaces (GenBank spec uses spaces)
      lines <- gsub("\t", "    ", lines, fixed = TRUE)

      # Split into records by // delimiter
      record_chunks <- split_records_by_slashes(lines)

      if (length(record_chunks) == 0) {
        cli::cli_abort("No valid GenBank records found in file: {path}")
      }

      # Parse each record
      records <- lapply(record_chunks, parse_gbk_record)

      # Normalize to gateway contract
      normalized <- lapply(records, normalize_gbk_record)

      normalized
    }
  )
}

#' Normalize GenBank Record to Gateway Contract
#'
#' @description
#' Transforms the parser output to the standard gateway contract format
#' expected by use cases.
#'
#' @param record Parsed GenBank record from parse_gbk_record()
#' @return Normalized record with metadata, features, sequence
#' @keywords internal
normalize_gbk_record <- function(record) {
  # Extract metadata
  meta <- record$metadata

  # Normalize LOCUS information
  locus_info <- meta$LOCUS %||% list()

  metadata_normalized <- list(
    locus = locus_info$locus_name %||% NA_character_,
    accession = meta$ACCESSION %||% NA_character_,
    version = meta$VERSION %||% NA_character_,
    length_bp = locus_info$length %||% NA_integer_,
    mol_type = locus_info$mol_type %||% NA_character_,
    topology = locus_info$topology %||% NA_character_,
    division = locus_info$division %||% NA_character_,
    date = locus_info$date %||% NA_character_,
    definition = meta$DEFINITION %||% NA_character_,
    organism = meta$ORGANISM %||% NA_character_,
    taxonomy = if (!is.null(meta$TAXONOMY)) {
      paste(meta$TAXONOMY, collapse = "; ")
    } else {
      NA_character_
    }
  )

  # Extract features (already a data.frame from parser)
  features <- record$features

  # Ensure features data.frame has expected columns
  if (!is.data.frame(features) || nrow(features) == 0) {
    # Create empty features data.frame with expected columns
    features <- data.frame(
      type = character(),
      start = integer(),
      end = integer(),
      strand = character(),
      gene = character(),
      locus_tag = character(),
      product = character(),
      stringsAsFactors = FALSE
    )
  }

  # Extract sequence
  sequence <- record$sequence %||% ""

  list(
    metadata = metadata_normalized,
    features = features,
    sequence = sequence
  )
}

# Helper: %||% operator (like rlang)
`%||%` <- function(x, y) if (is.null(x)) y else x

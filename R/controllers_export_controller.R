#' Write Sequences to FASTA File
#'
#' @description
#' Export sequences from a genome_entity or character vector to FASTA format.
#'
#' @param x A genome_entity object or named character vector
#' @param file Output file path
#' @param wrap_width Integer; wrap sequences at this width (default 80)
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisibly returns the file path
#' @export
#'
#' @examples
#' \dontrun{
#' genome <- read_genome("data.gbk")
#'
#' # Export all sequences
#' write_fasta(genome, "output.fasta")
#'
#' # Export specific sequences
#' seqs <- sequences(genome)
#' write_fasta(seqs[1:5], "subset.fasta")
#' }
write_fasta <- function(x, file, wrap_width = 80, ...) {
  # Handle different input types
  if (inherits(x, "genome_entity")) {
    # Extract sequences from genome_entity
    sequences_to_write <- x$sequences$dna_raw
  } else if (is.character(x)) {
    # Already a character vector
    sequences_to_write <- x
  } else {
    stop("x must be a genome_entity object or character vector", call. = FALSE)
  }

  # Create FASTA gateway
  gateway <- create_fasta_gateway(use_bioconductor = TRUE)

  # Write file
  gateway$write(sequences_to_write, file, wrap_width)

  invisible(file)
}

#' Write Features to GFF3 File
#'
#' @description
#' Export features from a genome_entity to GFF3 format.
#'
#' @param x A genome_entity object
#' @param file Output file path
#' @param source Character; source field for GFF3 (default "micromicon")
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisibly returns the file path
#' @export
#'
#' @examples
#' \dontrun{
#' genome <- read_genome("data.gbk")
#' write_gff3(genome, "output.gff3")
#' }
write_gff3 <- function(x, file, source = "micromicon", ...) {
  # Validate input
  if (!inherits(x, "genome_entity")) {
    stop("x must be a genome_entity object", call. = FALSE)
  }

  validate_genome_entity(x)

  # Get features
  feats <- x$features

  if (nrow(feats) == 0) {
    warning("No features to write", call. = FALSE)
    return(invisible(file))
  }

  # Open file for writing
  con <- file(file, open = "wt")
  on.exit(close(con))

  # Write GFF3 header
  writeLines("##gff-version 3", con)

  # Build GFF3 lines
  for (i in seq_len(nrow(feats))) {
    feat <- feats[i, ]

    # Build attributes string
    attrs <- build_gff3_attributes(feat)

    # Build GFF3 line (9 columns: seqname, source, type, start, end, score, strand, phase, attributes)
    gff_line <- paste(
      feat$seqname,
      source,
      feat$type,
      feat$start,
      feat$end,
      ".",  # score
      ifelse(is.na(feat$strand), ".", feat$strand),
      ".",  # phase
      attrs,
      sep = "\t"
    )

    writeLines(gff_line, con)
  }

  invisible(file)
}

#' Build GFF3 Attributes String
#'
#' @description
#' Helper function to build the attributes column (column 9) for GFF3 format.
#'
#' @param feat Single row from features data.frame
#' @return Character string with key=value pairs
#' @keywords internal
build_gff3_attributes <- function(feat) {
  attrs <- character()

  # Common attributes to include
  attr_fields <- c("ID", "Name", "gene", "locus_tag", "product", "protein_id", "Parent")

  for (field in attr_fields) {
    if (field %in% names(feat) && !is.na(feat[[field]]) && nzchar(feat[[field]])) {
      # URL encode value
      value <- utils::URLencode(feat[[field]], reserved = TRUE)
      attrs <- c(attrs, paste0(field, "=", value))
    }
  }

  if (length(attrs) == 0) {
    return(".")
  }

  paste(attrs, collapse = ";")
}

#' Write GenBank and FASTA (Combined Export)
#'
#' @description
#' Export genome to FASTA format. This is a convenience wrapper for write_fasta()
#' that maintains backward compatibility with the old API.
#'
#' @param x A genome_entity object
#' @param file Output file path
#' @param ... Additional arguments passed to write_fasta
#'
#' @return Invisibly returns the file path
#' @export
#'
#' @examples
#' \dontrun{
#' genome <- read_genome("data.gbk")
#' write_gbk_fasta(genome, "output.fasta")
#' }
write_gbk_fasta <- function(x, file, ...) {
  write_fasta(x, file, ...)
}

#' Write GenBank and GFF3 (Combined Export)
#'
#' @description
#' Export genome to GFF3 format. This is a convenience wrapper for write_gff3()
#' that maintains backward compatibility with the old API.
#'
#' @param x A genome_entity object
#' @param file Output file path
#' @param ... Additional arguments passed to write_gff3
#'
#' @return Invisibly returns the file path
#' @export
#'
#' @examples
#' \dontrun{
#' genome <- read_genome("data.gbk")
#' write_gbk_gff3(genome, "output.gff3")
#' }
write_gbk_gff3 <- function(x, file, ...) {
  write_gff3(x, file, ...)
}

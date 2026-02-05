#' Write Sequences to FASTA File
#'
#' @description
#' Export sequences from a genome_entity or character vector to FASTA format.
#'
#' **⚠️ METADATA LOSS**: If the genome was imported from GenBank, this export
#' will LOSE organism information, taxonomic lineage, references, comments,
#' and accession numbers. FASTA format only stores sequence IDs and sequences.
#'
#' **❌ NO REVERSE CONVERSION**: GFF3+FASTA cannot be converted back to GenBank.
#' No `write_genbank()` function exists - this is FORBIDDEN, not a missing feature.
#'
#' @details
#' ## What Gets Preserved
#' - Sequence IDs (from seqnames)
#' - DNA sequences
#'
#' ## What Gets LOST (from GenBank sources)
#' - Organism name and taxonomic lineage
#' - References and publication citations
#' - Comments and curator notes
#' - Accession numbers and version history
#' - Sequence topology (circular vs. linear)
#' - All feature annotations (genes, CDS, etc.)
#'
#' ## Why GenBank Export is FORBIDDEN
#' There is no reverse conversion (GFF3+FASTA → GenBank) because:
#' - GFF3+FASTA lacks the metadata required for valid GenBank files
#' - GenBank format requires organism, taxonomy, and reference information
#' - Creating incomplete GenBank files would violate NCBI format specifications
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
    # Use new export_genome generic for genome_entity objects
    export_genome(x, file, format = "fasta", wrap_width = wrap_width, ...)
  } else if (is.character(x)) {
    # Direct character vector export (backward compatibility)
    gateway <- create_fasta_gateway(use_bioconductor = TRUE)
    gateway$write(x, file, wrap_width)
    invisible(file)
  } else {
    cli::cli_abort("x must be a genome_entity object or character vector")
  }
}

#' Write Features to GFF3 File
#'
#' @description
#' Export features from a genome_entity to GFF3 format.
#'
#' **⚠️ METADATA LOSS**: If the genome was imported from GenBank, this export
#' will LOSE organism information, taxonomic lineage, references, comments,
#' and accession numbers. GFF3 format cannot represent GenBank metadata.
#'
#' **❌ NO REVERSE CONVERSION**: GFF3+FASTA cannot be converted back to GenBank.
#' No `write_genbank()` function exists - this is FORBIDDEN, not a missing feature.
#'
#' @details
#' ## What Gets Preserved
#' - Feature coordinates (start, end, strand)
#' - Feature types (gene, CDS, tRNA, etc.)
#' - Gene names, products, locus tags
#'
#' ## What Gets LOST (from GenBank sources)
#' - Organism name and taxonomic lineage
#' - References and publication citations
#' - Comments and curator notes
#' - Accession numbers and version history
#' - Sequence topology (circular vs. linear)
#'
#' ## Why GenBank Export is FORBIDDEN
#' There is no reverse conversion (GFF3+FASTA → GenBank) because:
#' - GFF3+FASTA lacks the metadata required for valid GenBank files
#' - GenBank format requires organism, taxonomy, and reference information
#' - Creating incomplete GenBank files would violate NCBI format specifications
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
    cli::cli_abort("x must be a genome_entity object")
  }

  # Use new export_genome generic
  export_genome(x, file, format = "gff3", source = source, ...)
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

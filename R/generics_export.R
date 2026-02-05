#' Export Genome to File
#'
#' @description
#' Generic function for exporting genome objects to various file formats.
#' Provides a unified interface with automatic format detection from file extensions.
#'
#' @details
#' ## Supported Formats
#' - `fasta` - FASTA format (sequences only)
#' - `gff3` - GFF3 format (annotations only)
#' - `auto` - Automatically detect format from file extension
#'
#' ## Metadata Loss Warnings
#' When exporting GenBank-sourced data to FASTA or GFF3, metadata will be lost:
#' - Organism information and taxonomic lineage
#' - References and publication citations
#' - Comments and curator notes
#' - Accession numbers and version history
#'
#' This is a ONE-WAY conversion. Keep your original GenBank file as the source of truth.
#'
#' @param x Genome object to export
#' @param file Output file path
#' @param format Export format ("auto", "fasta", "gff3"). Default "auto" detects from file extension.
#' @param ... Additional arguments passed to format-specific exporters
#'
#' @return Invisibly returns the file path
#' @export
#'
#' @examples
#' \dontrun{
#' genome <- read_genome("data.gbk")
#'
#' # Auto-detect format from extension
#' export_genome(genome, "output.fasta")      # FASTA
#' export_genome(genome, "output.gff3")       # GFF3
#'
#' # Explicit format specification
#' export_genome(genome, "output.txt", format = "fasta")
#' export_genome(genome, "annotations.txt", format = "gff3")
#'
#' # With additional options
#' export_genome(genome, "output.fasta", wrap_width = 60)
#' export_genome(genome, "output.gff3", source = "my_analysis")
#' }
export_genome <- function(x, file, format = c("auto", "fasta", "gff3"), ...) {
  UseMethod("export_genome")
}

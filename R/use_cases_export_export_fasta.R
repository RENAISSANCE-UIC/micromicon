#' Execute FASTA Export Use Case
#'
#' @description
#' Pure use case function for exporting genome sequences to FASTA format.
#' This is the core business logic extracted from the controller layer.
#'
#' @details
#' This use case:
#' 1. Validates the genome entity
#' 2. Checks for metadata loss warnings (GenBank sources)
#' 3. Extracts sequences from the genome entity
#' 4. Delegates to FASTA gateway for file writing
#'
#' ## Metadata Loss
#' When exporting GenBank-sourced genomes, this function warns about lost metadata:
#' - Organism information and taxonomic lineage
#' - References and citations
#' - Comments and curator notes
#' - Accession numbers
#' - Feature annotations
#'
#' @param entity A genome_entity object
#' @param file Output file path
#' @param wrap_width Integer; wrap sequences at this width (default 80)
#' @param warn_metadata Logical; whether to show metadata loss warnings (default TRUE)
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisibly returns the file path
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Internal use only - called by export_genome() generic
#' entity <- read_genome("data.gbk")
#' execute_export_fasta(entity, "output.fasta")
#' }
execute_export_fasta <- function(entity, file, wrap_width = 80,
                                 warn_metadata = TRUE, ...) {
  # Validate entity
  validate_genome_entity(entity)

  # Check for GenBank source and warn about metadata loss
  if (warn_metadata) {
    import_source <- attr(entity, "import_source")

    if (!is.null(import_source) && grepl("genbank", import_source, ignore.case = TRUE)) {
      # Check if warnings are enabled globally
      if (!isFALSE(getOption("micromicon.warn_export", default = TRUE))) {
        cli::cli_alert_warning("Exporting GenBank-sourced data to FASTA LOSES metadata")
        cli::cli_inform(c(
          "i" = "Lost: organism, taxonomy, references, comments, accession, features",
          "i" = "This conversion is ONE-WAY - no reverse conversion exists",
          "i" = "Keep your original GenBank file as source of truth",
          "i" = "Suppress: options(micromicon.warn_export = FALSE)"
        ))
      }
    }
  }

  # Extract sequences from genome_entity
  sequences_to_write <- entity$sequences$dna_raw

  # Create FASTA gateway
  gateway <- create_fasta_gateway(use_bioconductor = TRUE)

  # Write file
  gateway$write(sequences_to_write, file, wrap_width)

  invisible(file)
}

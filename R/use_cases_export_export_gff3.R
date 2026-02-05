#' Execute GFF3 Export Use Case
#'
#' @description
#' Pure use case function for exporting genome features to GFF3 format.
#' This is the core business logic extracted from the controller layer.
#'
#' @details
#' This use case:
#' 1. Validates the genome entity
#' 2. Checks for metadata loss warnings (GenBank sources)
#' 3. Extracts features from the genome entity
#' 4. Formats and writes GFF3 file with proper headers and attributes
#'
#' ## Metadata Loss
#' When exporting GenBank-sourced genomes, this function warns about lost metadata:
#' - Organism information and taxonomic lineage
#' - References and citations
#' - Comments and curator notes
#' - Accession numbers
#' - Sequence topology information
#'
#' @param entity A genome_entity object
#' @param file Output file path
#' @param source Character; source field for GFF3 (default "micromicon")
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
#' execute_export_gff3(entity, "output.gff3")
#' }
execute_export_gff3 <- function(entity, file, source = "micromicon",
                                warn_metadata = TRUE, ...) {
  # Validate entity
  validate_genome_entity(entity)

  # Check for GenBank source and warn about metadata loss
  if (warn_metadata) {
    import_source <- attr(entity, "import_source")

    if (!is.null(import_source) && grepl("genbank", import_source, ignore.case = TRUE)) {
      # Check if warnings are enabled globally
      if (!isFALSE(getOption("micromicon.warn_export", default = TRUE))) {
        cli::cli_alert_warning("Exporting GenBank-sourced data to GFF3 LOSES metadata")
        cli::cli_inform(c(
          "i" = "Lost: organism, taxonomy, references, comments, accession",
          "i" = "This conversion is ONE-WAY - no reverse conversion exists",
          "i" = "Keep your original GenBank file as source of truth",
          "i" = "Suppress: options(micromicon.warn_export = FALSE)"
        ))
      }
    }
  }

  # Get features
  feats <- entity$features

  if (nrow(feats) == 0) {
    cli::cli_warn("No features to write")
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

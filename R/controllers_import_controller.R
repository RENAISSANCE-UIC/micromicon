#' Read Genome Data from Various Formats
#'
#' @description
#' Unified entry point for reading genomic data from GenBank or GFF3+FASTA files.
#' Automatically detects format and returns a genome_entity object.
#'
#' This is a controller function that coordinates gateways and use cases.
#'
#' @param path Character string or named vector specifying input file(s):
#'   - Single file path with .gb/.gbk extension → GenBank format
#'   - Named vector c(gff = "path.gff3", fasta = "path.fasta") → GFF3+FASTA
#' @param gff Path to GFF3 file (alternative to using named vector in path)
#' @param fasta Path to FASTA file (alternative to using named vector in path)
#' @param format Character string to override auto-detection:
#'   "auto" (default), "genbank", or "gff3_fasta"
#' @param ... Additional arguments passed to use cases (e.g., auto_harmonize, verbose)
#'
#' @return A genome_entity object
#' @export
#'
#' @examples
#' \dontrun{
#' # GenBank format
#' genome1 <- read_genome("data/ecoli.gbk")
#'
#' # GFF3 + FASTA (named vector)
#' genome2 <- read_genome(c(gff = "data/ecoli.gff3", fasta = "data/ecoli.fna"))
#'
#' # GFF3 + FASTA (separate arguments)
#' genome3 <- read_genome(gff = "data/ecoli.gff3", fasta = "data/ecoli.fna")
#'
#' # Explicit format specification
#' genome4 <- read_genome("data/genome.txt", format = "genbank")
#' }
read_genome <- function(path = NULL, gff = NULL, fasta = NULL,
                        format = c("auto", "genbank", "gff3_fasta"), ...) {
  format <- match.arg(format)

  # Handle different input modes
  if (!is.null(gff) && !is.null(fasta)) {
    # Mode 1: gff= and fasta= arguments
    format <- "gff3_fasta"
    gff_path <- gff
    fasta_path <- fasta
  } else if (!is.null(path) && is.character(path) && length(path) == 2 &&
             !is.null(names(path)) && all(c("gff", "fasta") %in% names(path))) {
    # Mode 2: named vector c(gff = "...", fasta = "...")
    format <- "gff3_fasta"
    gff_path <- path["gff"]
    fasta_path <- path["fasta"]
  } else if (!is.null(path) && is.character(path) && length(path) == 1) {
    # Mode 3: single file path
    if (format == "auto") {
      # Auto-detect based on extension
      if (grepl("\\.(gb|gbk|genbank)$", path, ignore.case = TRUE)) {
        format <- "genbank"
      } else if (grepl("\\.gff3?$", path, ignore.case = TRUE)) {
        cli::cli_abort(c(
          "GFF3 file detected but FASTA file required.",
          "i" = "Provide both: read_genome(gff = '{path}', fasta = 'genome.fasta')"
        ))
      } else {
        cli::cli_abort(c(
          "Cannot auto-detect format from path: {path}",
          "i" = "Use format= argument or provide gff= and fasta= arguments"
        ))
      }
    }

    if (format == "genbank") {
      gbk_path <- path
    } else {
      cli::cli_abort("Single file path only works for GenBank format")
    }
  } else {
    cli::cli_abort(c(
      "Invalid arguments to read_genome().",
      "i" = "Provide either: path (GenBank file), or gff= and fasta= arguments, or c(gff = ..., fasta = ...)"
    ))
  }

  # Read based on detected format
  if (format == "genbank") {
    # Create GenBank gateway
    gateway <- create_genbank_gateway()

    # Call use case
    entity <- execute_import_genbank(gateway, gbk_path, options = list(...))
  } else if (format == "gff3_fasta") {
    # Create gateways
    gff_gateway <- create_gff_gateway(use_bioconductor = TRUE)
    fasta_gateway <- create_fasta_gateway(use_bioconductor = TRUE)

    # Call use case
    entity <- execute_import_gff_fasta(
      gff_gateway, fasta_gateway,
      gff_path, fasta_path,
      options = list(...)
    )
  }

  entity
}

#' Read GenBank File
#'
#' @description
#' Read a GenBank (.gb, .gbk) file. By default returns a genome_entity object.
#' Set return_entity=FALSE for backward compatibility with old list format.
#'
#' @param path Path to GenBank file
#' @param return_entity Logical; if TRUE (default), return genome_entity,
#'   if FALSE, return legacy list format
#' @param ... Additional arguments (currently unused)
#'
#' @return A genome_entity object (if return_entity=TRUE) or list of records (if FALSE)
#' @export
#'
#' @examples
#' \dontrun{
#' # New format (genome_entity)
#' genome <- read_gbk("data/ecoli.gbk")
#'
#' # Old format (for backward compatibility)
#' gbk_list <- read_gbk("data/ecoli.gbk", return_entity = FALSE)
#' }
read_gbk <- function(path, return_entity = TRUE, ...) {
  # Create gateway
  gateway <- create_genbank_gateway()

  if (return_entity) {
    # Return genome_entity
    entity <- execute_import_genbank(gateway, path, options = list(...))
    entity
  } else {
    # Return raw records (legacy format)
    records <- gateway$read(path)
    records
  }
}

#' Initialize Genome from GFF3 and FASTA Files
#'
#' @description
#' Initialize a genome object from GFF3 annotation and FASTA sequence files.
#' By default returns a genome_entity object. Set return_entity=FALSE for
#' backward compatibility with old format.
#'
#' @param gff_path Path to GFF3 file
#' @param fasta_path Path to FASTA file
#' @param return_entity Logical; if TRUE (default), return genome_entity,
#'   if FALSE, return legacy genome_obj format (not yet implemented)
#' @param auto_harmonize Logical; attempt to harmonize seqname mismatches
#' @param verbose Logical; if TRUE, print progress messages
#' @param ... Additional arguments passed to use case
#'
#' @return A genome_entity object (if return_entity=TRUE)
#' @export
#'
#' @examples
#' \dontrun{
#' # New format (genome_entity)
#' genome <- init_genome("data/ecoli.gff3", "data/ecoli.fna")
#'
#' # With options
#' genome <- init_genome(
#'   "data/ecoli.gff3",
#'   "data/ecoli.fna",
#'   auto_harmonize = TRUE,
#'   verbose = TRUE
#' )
#' }
init_genome <- function(gff_path, fasta_path,
                        return_entity = TRUE,
                        auto_harmonize = TRUE,
                        verbose = TRUE,
                        ...) {
  # Create gateways
  gff_gateway <- create_gff_gateway(use_bioconductor = TRUE)
  fasta_gateway <- create_fasta_gateway(use_bioconductor = TRUE)

  # Build options
  options <- list(
    auto_harmonize = auto_harmonize,
    verbose = verbose,
    ...
  )

  # Call use case
  entity <- execute_import_gff_fasta(
    gff_gateway, fasta_gateway,
    gff_path, fasta_path,
    options = options
  )

  if (return_entity) {
    entity
  } else {
    # Legacy format conversion (not yet implemented)
    cli::cli_warn("Legacy format (return_entity=FALSE) not yet implemented. Returning genome_entity.")
    entity
  }
}

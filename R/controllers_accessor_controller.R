#' Access Sequences from genome_entity
#'
#' @description
#' Extract sequences from a genome_entity object in various formats.
#'
#' @param x A genome_entity object
#' @param format Character; output format:
#'   - "character": Named character vector (default)
#'   - "DNAStringSet": Biostrings DNAStringSet (requires Biostrings)
#' @param ... Additional arguments (currently unused)
#'
#' @return Sequences in requested format
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' genome <- read_genome("data.gbk")
#'
#' # As character vector
#' seqs <- sequences(genome)
#'
#' # As DNAStringSet (if Biostrings available)
#' seqs_bio <- sequences(genome, format = "DNAStringSet")
#' }
sequences <- function(x, format = c("character", "DNAStringSet"), ...) {
  format <- match.arg(format)

  # Validate input
  if (!inherits(x, "genome_entity")) {
    cli::cli_abort("x must be a genome_entity object")
  }

  validate_genome_entity(x)

  if (format == "character") {
    # Return character vector
    return(x$sequences$dna_raw)
  } else if (format == "DNAStringSet") {
    # Check if Bioconductor gateway is available
    bio_gateway <- create_bioconductor_gateway()

    if (!bio_gateway$is_available()) {
      cli::cli_abort(c(
        "Biostrings package required for DNAStringSet format.",
        "i" = "Install with: BiocManager::install('Biostrings')"
      ))
    }

    # Convert to DNAStringSet
    if (!is.null(x$sequences$dna_bio)) {
      # Already cached
      return(x$sequences$dna_bio)
    } else {
      # Convert
      return(bio_gateway$to_dnastringset(x$sequences$dna_raw))
    }
  }
}

#' Access Features from genome_entity
#'
#' @description
#' Extract features from a genome_entity object as data.frame or GRanges.
#'
#' @param x A genome_entity object
#' @param format Character; output format:
#'   - "data.frame": data.frame with features (default)
#'   - "GRanges": GenomicRanges GRanges object (requires GenomicRanges)
#' @param type Character; filter by feature type (e.g., "gene", "CDS")
#' @param ... Additional arguments (currently unused)
#'
#' @return Features in requested format
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' genome <- read_genome("data.gbk")
#'
#' # All features as data.frame
#' feats <- features(genome)
#'
#' # Only genes
#' genes <- features(genome, type = "gene")
#'
#' # As GRanges (if GenomicRanges available)
#' gr <- features(genome, format = "GRanges")
#' }
features <- function(x, format = c("data.frame", "GRanges"), type = NULL, ...) {
  format <- match.arg(format)

  # Validate input
  if (!inherits(x, "genome_entity")) {
    cli::cli_abort("x must be a genome_entity object")
  }

  validate_genome_entity(x)

  # Get features
  feats <- x$features

  # Filter by type if requested
  if (!is.null(type) && "type" %in% names(feats)) {
    feats <- feats[feats$type == type, ]
  }

  if (format == "data.frame") {
    return(feats)
  } else if (format == "GRanges") {
    # Check if Bioconductor gateway is available
    bio_gateway <- create_bioconductor_gateway()

    if (!bio_gateway$is_available()) {
      cli::cli_abort(c(
        "GenomicRanges package required for GRanges format.",
        "i" = "Install with: BiocManager::install('GenomicRanges')"
      ))
    }

    # Convert to GRanges
    return(bio_gateway$to_granges(feats))
  }
}

#' Access Metadata from genome_entity
#'
#' @description
#' Extract metadata (sequence-level information) from a genome_entity object.
#'
#' @param x A genome_entity object
#' @param ... Additional arguments (currently unused)
#'
#' @return data.frame with metadata
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' genome <- read_genome("data.gbk")
#' meta <- metadata(genome)
#' print(meta)
#' }
metadata <- function(x, ...) {
  # Validate input
  if (!inherits(x, "genome_entity")) {
    cli::cli_abort("x must be a genome_entity object")
  }

  validate_genome_entity(x)

  x$metadata
}

#' Get Sequence Names (Scaffolds/Contigs/Chromosomes/Plasmids)
#'
#' @description
#' Convenience function to quickly retrieve the names of all sequences
#' (scaffolds, contigs, chromosomes, plasmids) in a genome object.
#'
#' This is a short alias for `get_genome_seqnames()` and is equivalent to
#' `unique(genome$features$seqname)` but more convenient.
#'
#' @param x A genome_entity object
#' @param ... Additional arguments (currently unused)
#'
#' @return Character vector of unique sequence names
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' gd <- read_genome("genome.gbk")
#'
#' # Get all sequence names
#' seqnames(gd)  # e.g., c("NC_000913.3", "plasmid_pXYZ")
#'
#' # Also available as get_genome_seqnames()
#' get_genome_seqnames(gd)
#' }
seqnames <- function(x, ...) {
  # Validate input
  if (!inherits(x, "genome_entity")) {
    cli::cli_abort("x must be a genome_entity object")
  }

  validate_genome_entity(x)

  x$indices$seqnames
}

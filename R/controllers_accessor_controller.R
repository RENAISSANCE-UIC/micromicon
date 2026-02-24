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


#' Access Features
#'
#' @param x Object to extract features from
#' @param ... Additional arguments
#' @export
get_features <- function(x, ...) {
  UseMethod("get_features")
}

#' @export
get_features.genome_entity <- function(x, format = c("tibble", "data.frame", "GRanges"), type = NULL, ...) {
  format <- match.arg(format)
  validate_genome_entity(x)

  feats <- x$features

  if (!is.null(type) && "type" %in% names(feats)) {
    feats <- feats[feats$type == type, ]
  }

  if (format == "tibble") {
    return(tibble::as_tibble(feats))
  } else if (format == "data.frame") {
    return(feats)
  } else {
    bio_gateway <- create_bioconductor_gateway()
    if (!bio_gateway$is_available()) {
      cli::cli_abort(c(
        "GenomicRanges required for GRanges format",
        "i" = "Install with: BiocManager::install('GenomicRanges')"
      ))
    }
    return(bio_gateway$to_granges(feats))
  }
}

#' @export
get_features.default <- function(x, ...) {
  cli::cli_abort("get_features() not implemented for class {.cls {class(x)[1]}}")
}


#' Access Sequences
#'
#' @param x Object to extract sequences from
#' @param ... Additional arguments
#' @export
sequences <- function(x, ...) {
  UseMethod("sequences")
}

#' @export
sequences.genome_entity <- function(x, format = c("character", "DNAStringSet"), ...) {
  format <- match.arg(format)
  validate_genome_entity(x)

  if (format == "character") {
    return(x$sequences$dna_raw)
  } else {
    bio_gateway <- create_bioconductor_gateway()
    if (!bio_gateway$is_available()) {
      cli::cli_abort(c(
        "Biostrings required for DNAStringSet format",
        "i" = "Install with: BiocManager::install('Biostrings')"
      ))
    }
    return(bio_gateway$to_dnastringset(x$sequences$dna_raw))
  }
}

#' @export
sequences.default <- function(x, ...) {
  cli::cli_abort("sequences() not implemented for class {.cls {class(x)[1]}}")
}


#' Access Genome Metadata
#'
#' @param x Object to extract metadata from
#' @param ... Additional arguments
#' @export
get_genome_metadata <- function(x, ...) {
  UseMethod("get_genome_metadata")
}

#' @export
get_genome_metadata.genome_entity <- function(x, ...) {
  validate_genome_entity(x)
  x$metadata
}

#' @export
get_genome_metadata.default <- function(x, ...) {
  cli::cli_abort("get_genome_metadata() not implemented for class {.cls {class(x)[1]}}")
}


#' Get Sequence Names (Scaffolds/Contigs/Chromosomes/Plasmids)
#'
#' @description
#' Convenience function to quickly retrieve the names of all sequences
#' (scaffolds, contigs, chromosomes, plasmids) in a genome object.
#'
#' This is equivalent to `unique(genome$features$seqname)` but more convenient
#' and standardized across genome objects.
#'
#' @param x A genome object (e.g., genome_entity)
#' @param ... Additional arguments passed to methods
#'
#' @return Character vector of unique sequence names
#'
#' @examples
#' \dontrun{
#' gd <- read_genome("genome.gbk")
#'
#' # Full function name
#' get_genome_seqnames(gd)
#'
#' # Short alias
#' seqnames(gd)
#' }
#'
#' @export
get_genome_seqnames <- function(x, ...) {
  UseMethod("get_genome_seqnames")
}

#' @export
get_genome_seqnames.genome_entity <- function(x, ...) {
  validate_genome_entity(x)
  x$indices$seqnames
}

#' @export
get_genome_seqnames.default <- function(x, ...) {
  cli::cli_abort("get_genome_seqnames() not implemented for class {.cls {class(x)[1]}}")
}

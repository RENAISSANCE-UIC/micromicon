#' Access Features
#'
#' @param x Object to extract features from
#' @param ... Additional arguments
#' @export
features <- function(x, ...) {
  UseMethod("features")
}

#' @export
features.genome_entity <- function(x, format = c("data.frame", "GRanges"), type = NULL, ...) {
  format <- match.arg(format)
  validate_genome_entity(x)

  feats <- x$features

  if (!is.null(type) && "type" %in% names(feats)) {
    feats <- feats[feats$type == type, ]
  }

  if (format == "data.frame") {
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
features.default <- function(x, ...) {
  cli::cli_abort("features() not implemented for class {.cls {class(x)[1]}}")
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
genome_metadata <- function(x, ...) {
  UseMethod("genome_metadata")
}

#' @export
genome_metadata.genome_entity <- function(x, ...) {
  validate_genome_entity(x)
  x$metadata
}

#' @export
genome_metadata.default <- function(x, ...) {
  cli::cli_abort("genome_metadata() not implemented for class {.cls {class(x)[1]}}")
}


#' Get Sequence Names
#'
#' @param x Object to extract sequence names from
#' @param ... Additional arguments
#' @export
genome_seqnames <- function(x, ...) {
  UseMethod("genome_seqnames")
}

#' @export
genome_seqnames.genome_entity <- function(x, ...) {
  validate_genome_entity(x)
  x$indices$seqnames
}

#' @export
genome_seqnames.default <- function(x, ...) {
  cli::cli_abort("genome_seqnames() not implemented for class {.cls {class(x)[1]}}")
}

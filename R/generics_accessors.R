#' Retrieve annotated features from a genome object
#'
#' @description
#' Returns the complete feature table of a loaded genome -- every annotated
#' element (genes, CDS, rRNAs, tRNAs, misc features) with its coordinates,
#' strand, type, and available identifiers.
#'
#' This is the canonical first stop after \code{\link{read_genome}()}: take
#' stock of what the annotation contains before querying individual elements.
#' The returned tibble has one row per feature. Columns present depend on
#' the source format, but the following are always available:
#'
#' | Column | Type | Description |
#' |--------|------|-------------|
#' | `seqname` | character | Contig or chromosome name |
#' | `start`, `end` | integer | 1-based, inclusive coordinates |
#' | `strand` | character | `"+"`, `"-"`, or `"*"` |
#' | `type` | character | Feature type (`CDS`, `gene`, `rRNA`, `tRNA`, ...) |
#' | `gene` | character | Gene symbol, when present in the annotation |
#' | `locus_tag` | character | Systematic locus identifier |
#' | `product` | character | Functional description |
#'
#' @param x A \code{genome_entity} returned by \code{\link{read_genome}()}.
#'   Works equally on \code{genome_entity_gd} objects.
#' @param format Output format. One of:
#'   \describe{
#'     \item{`"tibble"`}{A tidy tibble (default). Best for interactive exploration
#'       and downstream \code{dplyr} pipelines.}
#'     \item{`"data.frame"`}{A plain base-R \code{data.frame}.}
#'     \item{`"GRanges"`}{A \code{GenomicRanges::GRanges} object for Bioconductor
#'       workflows. Requires \code{BiocManager::install("GenomicRanges")}.}
#'   }
#' @param type Character. Restrict output to a single feature type, e.g.
#'   \code{"CDS"}, \code{"rRNA"}, \code{"tRNA"}, or \code{"gene"}.
#'   \code{NULL} (default) returns every feature type.
#' @param seqname Character. Restrict output to features on a specific contig
#'   or chromosome. \code{NULL} (default) returns features from all contigs.
#' @param gene Character. Restrict output to features with this gene symbol.
#'   \code{NULL} (default) does not filter by gene.
#' @param locus_tag Character. Restrict output to features with this locus tag.
#'   \code{NULL} (default) does not filter by locus_tag.
#' @param ... Reserved for future arguments; currently ignored.
#'
#' @return
#' A tibble (or data.frame / GRanges, per \code{format}) with one row per
#' genomic feature. The seven columns listed in the Description are always
#' present; annotation-specific extras (\code{protein_id}, \code{note}, etc.)
#' appear when the source file supplies them.
#'
#' @seealso
#' \code{\link{search_features}()} to filter by name, type, or coordinates;
#' \code{\link{get_roi_features}()} for features overlapping a specific window;
#' \code{\link{gene_info}()} for structured metadata on a single named gene.
#'
#' @examples
#' \dontrun{
#' ref <- read_genome("reference.gbk")
#'
#' # Census: what feature types does this annotation contain?
#' feats <- get_features(ref)
#' table(feats$type)
#'
#' # Protein-coding genes only
#' cds <- get_features(ref, type = "CDS")
#' nrow(cds)
#'
#' # Quick annotation overview
#' cds[, c("gene", "locus_tag", "product")]
#'
#' # Strand balance
#' table(cds$strand)
#'
#' # Ribosomal RNA loci -- useful for coverage QC
#' get_features(ref, type = "rRNA")
#'
#' # Filter by gene symbol
#' get_features(ref, gene = "motA")
#'
#' # Filter by contig and gene
#' get_features(ref, seqname = "U00096", gene = "motA")
#'
#' # Filter by locus tag
#' get_features(ref, locus_tag = "b1891")
#'
#' # Hand off to a Bioconductor workflow
#' get_features(ref, format = "GRanges")
#' }
#'
#' @export
get_features <- function(x, ...) {
  UseMethod("get_features")
}

#' @rdname get_features
#' @export
get_features.genome_entity <- function(x, format = c("tibble", "data.frame", "GRanges"), type = NULL, seqname = NULL, gene = NULL, locus_tag = NULL, ...) {
  format <- match.arg(format)
  validate_genome_entity(x)

  feats <- x$features

  if (!is.null(type) && "type" %in% names(feats)) {
    feats <- feats[feats$type == type, ]
  }

  if (!is.null(seqname) && "seqname" %in% names(feats)) {
    feats <- feats[!is.na(feats$seqname) & feats$seqname == seqname, ]
  }

  if (!is.null(gene) && "gene" %in% names(feats)) {
    feats <- feats[!is.na(feats$gene) & feats$gene == gene, ]
  }

  if (!is.null(locus_tag) && "locus_tag" %in% names(feats)) {
    feats <- feats[!is.na(feats$locus_tag) & feats$locus_tag == locus_tag, ]
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


#' Get Contig Sequences
#'
#' @param x Object to extract sequences from
#' @param ... Additional arguments
#' @export
get_contig_sequences <- function(x, ...) {
  UseMethod("get_contig_sequences")
}

#' @export
get_contig_sequences.genome_entity <- function(x, format = c("character", "DNAStringSet"), ...) {
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
get_contig_sequences.default <- function(x, ...) {
  cli::cli_abort("get_contig_sequences() not implemented for class {.cls {class(x)[1]}}")
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


#' Get Contig Names
#'
#' @description
#' Retrieve the names of all contigs, chromosomes, scaffolds, or plasmids
#' in a genome object.
#'
#' @param x A genome object (e.g., genome_entity)
#' @param ... Additional arguments passed to methods
#'
#' @return Character vector of contig/sequence names
#'
#' @examples
#' \dontrun{
#' gd <- read_genome("genome.gbk")
#' get_contig_names(gd)
#' }
#'
#' @export
get_contig_names <- function(x, ...) {
  UseMethod("get_contig_names")
}

#' @export
get_contig_names.genome_entity <- function(x, ...) {
  validate_genome_entity(x)
  x$indices$seqnames
}

#' @export
get_contig_names.default <- function(x, ...) {
  cli::cli_abort("get_contig_names() not implemented for class {.cls {class(x)[1]}}")
}

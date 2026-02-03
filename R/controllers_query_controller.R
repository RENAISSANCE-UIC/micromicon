#' Extract Sequences by Feature Name
#'
#' @description
#' Extract DNA sequences for features matching a pattern. Searches features
#' by gene name, locus_tag, or product and extracts their sequences.
#'
#' This is a controller function that calls the extract_sequences_by_name use case.
#'
#' @param x A genome_entity object
#' @param pattern Character string or regex pattern to match feature names
#' @param fields Character vector of fields to search (default: c("gene", "locus_tag", "product"))
#' @param type Character; filter by feature type (e.g., "CDS", "gene")
#' @param ignore_case Logical; case-insensitive matching (default TRUE)
#' @param translate Logical; translate CDS to protein (default FALSE)
#' @param ... Additional arguments (currently unused)
#'
#' @return Named character vector of sequences
#' @export
#'
#' @examples
#' \dontrun{
#' genome <- read_genome("data.gbk")
#'
#' # Extract sequences for genes matching "dnaA"
#' seqs <- extract_sequences_by_name(genome, "dnaA")
#'
#' # Extract CDS sequences for locus tags starting with "EC"
#' seqs <- extract_sequences_by_name(
#'   genome,
#'   "^EC",
#'   fields = "locus_tag",
#'   type = "CDS"
#' )
#'
#' # Translate to protein
#' proteins <- extract_sequences_by_name(genome, "dnaA", type = "CDS", translate = TRUE)
#' }
extract_sequences_by_name <- function(x, pattern,
                                      fields = c("gene", "locus_tag", "product"),
                                      type = NULL,
                                      ignore_case = TRUE,
                                      translate = FALSE,
                                      ...) {
  # Validate input
  if (!inherits(x, "genome_entity")) {
    stop("x must be a genome_entity object", call. = FALSE)
  }

  # Build options
  options <- list(
    ignore_case = ignore_case,
    type_filter = type,
    translate = translate
  )

  # Call use case
  execute_extract_sequences_by_name(x, pattern, fields, options)
}

#' Extract Sequences by Genomic Coordinates
#'
#' @description
#' Extract DNA sequences from specified genomic coordinates. Supports extracting
#' regions from one or more sequences with optional reverse complement.
#'
#' This is a controller function that calls the extract_sequences_by_coords use case.
#'
#' @param x A genome_entity object
#' @param seqname Character; sequence name(s)
#' @param start Integer; start position(s) (1-based, inclusive)
#' @param end Integer; end position(s) (1-based, inclusive)
#' @param strand Character; strand ("+", "-", or NULL for forward strand)
#' @param translate Logical; translate to protein (default FALSE)
#' @param names Character vector of names for output sequences (optional)
#' @param ... Additional arguments (currently unused)
#'
#' @return Character vector of sequences (named if names provided)
#' @export
#'
#' @examples
#' \dontrun{
#' genome <- read_genome("data.gbk")
#'
#' # Extract single region
#' seq <- extract_sequences_by_coords(genome, "chr1", 1, 100)
#'
#' # Extract multiple regions
#' seqs <- extract_sequences_by_coords(
#'   genome,
#'   seqname = c("chr1", "chr1", "plasmid"),
#'   start = c(1, 500, 1),
#'   end = c(100, 600, 200)
#' )
#'
#' # Extract with reverse complement
#' seq <- extract_sequences_by_coords(
#'   genome,
#'   "chr1", 1, 100,
#'   strand = "-"
#' )
#' }
extract_sequences_by_coords <- function(x, seqname, start, end,
                                        strand = "+",
                                        translate = FALSE,
                                        names = NULL,
                                        ...) {
  # Validate input
  if (!inherits(x, "genome_entity")) {
    stop("x must be a genome_entity object", call. = FALSE)
  }

  # Build options
  options <- list(
    strand = strand,
    translate = translate,
    names = names
  )

  # Call use case
  execute_extract_sequences_by_coords(x, seqname, start, end, options)
}

#' Search Features
#'
#' @description
#' Search for features matching criteria (type, name, coordinates).
#'
#' @param x A genome_entity object
#' @param type Character; filter by feature type (e.g., "gene", "CDS")
#' @param pattern Character; regex pattern to match in gene, locus_tag, or product
#' @param seqname Character; filter by sequence name
#' @param start Integer; filter features starting after this position
#' @param end Integer; filter features ending before this position
#' @param strand Character; filter by strand ("+", "-")
#' @param ... Additional arguments (currently unused)
#'
#' @return data.frame of matching features
#' @export
#'
#' @examples
#' \dontrun{
#' genome <- read_genome("data.gbk")
#'
#' # Find all genes
#' genes <- search_features(genome, type = "gene")
#'
#' # Find CDS on chr1
#' cds_chr1 <- search_features(genome, type = "CDS", seqname = "chr1")
#'
#' # Find features matching pattern
#' dna_features <- search_features(genome, pattern = "dna")
#' }
search_features <- function(x, type = NULL, pattern = NULL,
                           seqname = NULL, start = NULL, end = NULL,
                           strand = NULL, ...) {
  # Validate input
  if (!inherits(x, "genome_entity")) {
    stop("x must be a genome_entity object", call. = FALSE)
  }

  validate_genome_entity(x)

  # Get features
  feats <- x$features

  # Apply filters
  if (!is.null(type) && "type" %in% names(feats)) {
    feats <- feats[feats$type == type, ]
  }

  if (!is.null(pattern)) {
    # Search in gene, locus_tag, product fields
    matches <- rep(FALSE, nrow(feats))

    for (field in c("gene", "locus_tag", "product")) {
      if (field %in% names(feats)) {
        field_matches <- grepl(pattern, feats[[field]], ignore.case = TRUE)
        matches <- matches | field_matches
      }
    }

    feats <- feats[matches, ]
  }

  if (!is.null(seqname) && "seqname" %in% names(feats)) {
    feats <- feats[feats$seqname == seqname, ]
  }

  if (!is.null(start) && "start" %in% names(feats)) {
    feats <- feats[feats$start >= start, ]
  }

  if (!is.null(end) && "end" %in% names(feats)) {
    feats <- feats[feats$end <= end, ]
  }

  if (!is.null(strand) && "strand" %in% names(feats)) {
    feats <- feats[feats$strand == strand, ]
  }

  feats
}

#' Domain Validators for Genomic Data
#'
#' @description
#' Pure validation functions that enforce business rules for genomic entities.
#' These validators are framework-independent and contain no I/O or external dependencies.
#'
#' @name validators
NULL

#' Validate Genomic Coordinates
#'
#' @description
#' Validates that genomic coordinates satisfy domain rules:
#' - Start position must be positive (1-based)
#' - Start must be <= end
#' - Coordinates must be within sequence bounds (if max_length provided)
#'
#' @param start Integer start position (1-based)
#' @param end Integer end position (1-based, inclusive)
#' @param seqname Character sequence name (for error messages)
#' @param max_length Integer maximum sequence length (optional)
#'
#' @return TRUE if valid, stops with error if invalid
#' @export
#' @examples
#' # Valid coordinates
#' validate_genomic_coordinates(1, 100, "chr1")
#' validate_genomic_coordinates(50, 150, "chr1", max_length = 200)
#'
#' \dontrun{
#' # Invalid: negative start
#' validate_genomic_coordinates(-1, 100, "chr1")  # Error
#'
#' # Invalid: start > end
#' validate_genomic_coordinates(100, 50, "chr1")  # Error
#'
#' # Invalid: exceeds bounds
#' validate_genomic_coordinates(1, 200, "chr1", max_length = 100)  # Error
#' }
validate_genomic_coordinates <- function(start, end, seqname = NULL, max_length = NULL) {
  # Rule 1: Start must be positive (1-based coordinates)
  if (any(start < 1, na.rm = TRUE)) {
    msg <- "Genomic coordinates must be positive (1-based indexing)"
    if (!is.null(seqname)) {
      msg <- paste0(msg, " for sequence '", seqname, "'")
    }
    stop(msg, call. = FALSE)
  }

  # Rule 2: Start must be <= end
  if (any(start > end, na.rm = TRUE)) {
    msg <- "Start coordinate must be <= end coordinate"
    if (!is.null(seqname)) {
      msg <- paste0(msg, " for sequence '", seqname, "'")
    }
    stop(msg, call. = FALSE)
  }

  # Rule 3: Coordinates must be within sequence bounds (if provided)
  if (!is.null(max_length)) {
    if (any(end > max_length, na.rm = TRUE)) {
      msg <- paste0("Coordinates exceed sequence length (", max_length, " bp)")
      if (!is.null(seqname)) {
        msg <- paste0(msg, " for sequence '", seqname, "'")
      }
      stop(msg, call. = FALSE)
    }
  }

  TRUE
}

#' Validate Strand
#'
#' @description
#' Validates that strand value is one of the allowed genomic strand values.
#' Allowed values: "+", "-", ".", "*", NA
#'
#' @param strand Character strand value
#'
#' @return TRUE if valid, stops with error if invalid
#' @export
#' @examples
#' # Valid strands
#' validate_strand("+")
#' validate_strand("-")
#' validate_strand(".")  # Unknown
#' validate_strand("*")  # Either strand
#' validate_strand(NA)   # Not applicable
#'
#' \dontrun{
#' # Invalid strand
#' validate_strand("F")  # Error
#' validate_strand("forward")  # Error
#' }
validate_strand <- function(strand) {
  valid_strands <- c("+", "-", ".", "*", NA_character_)

  # Check each strand value
  invalid <- !strand %in% valid_strands & !is.na(strand)

  if (any(invalid)) {
    invalid_values <- unique(strand[invalid])
    stop(
      "Invalid strand value(s): ", paste(invalid_values, collapse = ", "),
      ". Allowed: +, -, ., *, NA",
      call. = FALSE
    )
  }

  TRUE
}

#' Validate DNA Sequence
#'
#' @description
#' Validates that a DNA sequence string contains only valid IUPAC nucleotide codes.
#' Allows both standard codes (A, C, G, T, U) and ambiguity codes (R, Y, S, W, K, M, etc.).
#'
#' Valid IUPAC nucleotide codes:
#' - A, C, G, T, U: Standard bases
#' - R (A or G), Y (C or T), S (G or C), W (A or T), K (G or T), M (A or C)
#' - B (C, G, or T), D (A, G, or T), H (A, C, or T), V (A, C, or G)
#' - N: Any base
#' - -: Gap
#'
#' @param seq_char Character string of DNA sequence
#' @param allow_gaps Logical; allow gap characters ("-")? Default TRUE
#'
#' @return TRUE if valid, stops with error if invalid
#' @export
#' @examples
#' # Valid sequences
#' validate_dna_sequence("ATCG")
#' validate_dna_sequence("ATCGNNNATCG")  # With ambiguity codes
#' validate_dna_sequence("ATCG--ATCG")   # With gaps
#'
#' \dontrun{
#' # Invalid sequence
#' validate_dna_sequence("ATCGX")  # X not allowed
#' validate_dna_sequence("ATCG123")  # Numbers not allowed
#' }
validate_dna_sequence <- function(seq_char, allow_gaps = TRUE) {
  # IUPAC nucleotide codes
  valid_codes <- c(
    "A", "C", "G", "T", "U",           # Standard bases
    "R", "Y", "S", "W", "K", "M",      # Purines/pyrimidines
    "B", "D", "H", "V",                # 3-way ambiguity
    "N"                                 # Any base
  )

  if (allow_gaps) {
    valid_codes <- c(valid_codes, "-")
  }

  # Convert to uppercase and split into characters
  seq_upper <- toupper(strsplit(seq_char, "")[[1]])

  # Check for invalid codes
  invalid <- !seq_upper %in% valid_codes
  if (any(invalid)) {
    invalid_codes <- unique(seq_upper[invalid])
    stop(
      "Invalid nucleotide code(s): ", paste(invalid_codes, collapse = ", "),
      ". Allowed: ", paste(valid_codes, collapse = ", "),
      call. = FALSE
    )
  }

  TRUE
}

#' Validate Feature Type
#'
#' @description
#' Validates that a feature type is a recognized genomic feature type.
#' This uses a permissive list that includes common Sequence Ontology (SO) terms.
#'
#' @param type Character feature type
#' @param allow_custom Logical; allow custom (unrecognized) types? Default TRUE
#'
#' @return TRUE if valid, emits warning if unrecognized and allow_custom=TRUE
#' @export
#' @examples
#' # Valid types
#' validate_feature_type("gene")
#' validate_feature_type("CDS")
#' validate_feature_type("mRNA")
#'
#' # Custom types (with warning)
#' validate_feature_type("my_custom_feature")
validate_feature_type <- function(type, allow_custom = TRUE) {
  # Common Sequence Ontology (SO) feature types
  recognized_types <- c(
    "gene", "CDS", "mRNA", "tRNA", "rRNA", "ncRNA",
    "exon", "intron", "UTR", "five_prime_UTR", "three_prime_UTR",
    "promoter", "enhancer", "terminator",
    "repeat_region", "mobile_element",
    "misc_feature", "misc_binding", "misc_RNA",
    "region", "source"
  )

  if (!type %in% recognized_types) {
    if (allow_custom) {
      warning(
        "Unrecognized feature type: '", type, "'. ",
        "Using anyway (allow_custom=TRUE). ",
        "Common types: ", paste(head(recognized_types, 10), collapse = ", "),
        call. = FALSE
      )
    } else {
      stop(
        "Unrecognized feature type: '", type, "'. ",
        "Allowed types: ", paste(recognized_types, collapse = ", "),
        call. = FALSE
      )
    }
  }

  TRUE
}

#' Validate Locus Tag Format
#'
#' @description
#' Validates that a locus tag follows common naming conventions.
#' Typically: alphanumeric with underscores, no spaces, reasonable length.
#'
#' @param locus_tag Character locus tag
#' @param max_length Integer maximum length (default 50)
#'
#' @return TRUE if valid, stops with error if invalid
#' @export
#' @examples
#' # Valid locus tags
#' validate_locus_tag("GENE001")
#' validate_locus_tag("ABC_12345")
#' validate_locus_tag("locus_tag_v2")
#'
#' \dontrun{
#' # Invalid locus tags
#' validate_locus_tag("GENE 001")  # Space not allowed
#' validate_locus_tag("GENE-001")  # Hyphen typically not allowed
#' }
validate_locus_tag <- function(locus_tag, max_length = 50) {
  # Check length
  if (nchar(locus_tag) > max_length) {
    stop(
      "Locus tag too long (", nchar(locus_tag), " characters). ",
      "Maximum: ", max_length,
      call. = FALSE
    )
  }

  # Check format: alphanumeric with underscores only
  if (!grepl("^[A-Za-z0-9_]+$", locus_tag)) {
    stop(
      "Invalid locus tag format: '", locus_tag, "'. ",
      "Must contain only letters, numbers, and underscores.",
      call. = FALSE
    )
  }

  TRUE
}

#' Validate Sequence Topology
#'
#' @description
#' Validates that sequence topology is one of the allowed values.
#' Allowed: "linear", "circular", NA
#'
#' @param topology Character topology value
#'
#' @return TRUE if valid, stops with error if invalid
#' @export
#' @examples
#' validate_sequence_topology("linear")
#' validate_sequence_topology("circular")
#' validate_sequence_topology(NA)
#'
#' \dontrun{
#' validate_sequence_topology("supercoiled")  # Error
#' }
validate_sequence_topology <- function(topology) {
  valid_topologies <- c("linear", "circular", NA_character_)

  if (!topology %in% valid_topologies & !is.na(topology)) {
    stop(
      "Invalid topology: '", topology, "'. ",
      "Allowed: linear, circular, NA",
      call. = FALSE
    )
  }

  TRUE
}

#' Create GFF3 Gateway
#'
#' @description
#' Creates a gateway for reading GFF3 format annotation files. This gateway
#' provides optional Bioconductor integration via rtracklayer for robust parsing.
#'
#' Gateway Contract (for use cases):
#' - read(path): Returns data.frame with columns:
#'   seqname, start, end, strand, type, and optional annotation columns
#'   (gene, locus_tag, product, etc.)
#'
#' @param use_bioconductor Logical; if TRUE and rtracklayer is available, use it.
#'   Default TRUE.
#'
#' @return Gateway object (list with read method)
#' @export
#' @examples
#' \dontrun{
#' gateway <- create_gff_gateway()
#' features <- gateway$read("annotation.gff3")
#' }
create_gff_gateway <- function(use_bioconductor = TRUE) {
  # Check if rtracklayer is available
  has_rtracklayer <- use_bioconductor &&
                     requireNamespace("rtracklayer", quietly = TRUE)

  list(
    read = function(path) {
      # Validate file exists
      if (!file.exists(path)) {
        stop("GFF3 file not found: ", path, call. = FALSE)
      }

      if (has_rtracklayer) {
        # Use rtracklayer for reading (robust, handles edge cases)
        features <- read_gff_with_rtracklayer(path)
      } else {
        # Use simple parser
        features <- parse_gff_simple(path)
      }

      features
    }
  )
}

#' Read GFF3 with rtracklayer
#'
#' @description
#' Uses Bioconductor's rtracklayer to read GFF3 files into a data.frame.
#' Converts GRanges to data.frame format expected by use cases.
#'
#' @param path Path to GFF3 file
#' @return data.frame with features
#' @keywords internal
read_gff_with_rtracklayer <- function(path) {
  # Import as GRanges
  gr <- rtracklayer::import(path, format = "gff3")

  # Convert to data.frame
  df <- as.data.frame(gr, stringsAsFactors = FALSE)

  # Rename columns to match contract
  if ("seqnames" %in% names(df)) {
    names(df)[names(df) == "seqnames"] <- "seqname"
  }

  # Ensure required columns exist
  required_cols <- c("seqname", "start", "end", "strand", "type")
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    stop("GFF3 file missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # Convert strand to character
  df$strand <- as.character(df$strand)

  # Drop GRanges-specific columns
  cols_to_drop <- c("width", "seqlevels", "seqlengths", "isCircular", "genome")
  df <- df[, !names(df) %in% cols_to_drop, drop = FALSE]

  df
}

#' Simple GFF3 Parser
#'
#' @description
#' Simple parser for GFF3 files when rtracklayer is not available.
#' Handles basic GFF3 format without advanced features.
#'
#' @param path Path to GFF3 file
#' @return data.frame with features
#' @keywords internal
parse_gff_simple <- function(path) {
  lines <- readLines(path, warn = FALSE)

  # Remove comment lines and blank lines
  is_comment <- grepl("^\\s*#", lines)
  is_blank <- !nzchar(trimws(lines))
  lines <- lines[!(is_comment | is_blank)]

  if (length(lines) == 0) {
    stop("No feature lines found in GFF3 file: ", path, call. = FALSE)
  }

  # Split by tabs
  parts <- strsplit(lines, "\t", fixed = TRUE)
  lens <- lengths(parts)

  # Keep only lines with at least 9 fields (GFF3 spec)
  valid <- lens >= 9
  if (sum(valid) == 0) {
    stop("No valid GFF3 lines found (need 9 tab-separated fields)", call. = FALSE)
  }

  parts <- parts[valid]

  # Build data.frame
  n_features <- length(parts)

  df <- data.frame(
    seqname = character(n_features),
    source = character(n_features),
    type = character(n_features),
    start = integer(n_features),
    end = integer(n_features),
    score = character(n_features),
    strand = character(n_features),
    phase = character(n_features),
    attributes = character(n_features),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_features)) {
    fields <- parts[[i]]
    df$seqname[i] <- fields[1]
    df$source[i] <- fields[2]
    df$type[i] <- fields[3]
    df$start[i] <- as.integer(fields[4])
    df$end[i] <- as.integer(fields[5])
    df$score[i] <- fields[6]
    df$strand[i] <- fields[7]
    df$phase[i] <- fields[8]
    df$attributes[i] <- fields[9]
  }

  # Parse attributes column (key=value pairs)
  df <- parse_gff_attributes(df)

  # Remove attributes column (now parsed into separate columns)
  df$attributes <- NULL

  df
}

#' Parse GFF3 Attributes
#'
#' @description
#' Parses the GFF3 attributes column (column 9) into separate columns.
#' Common attributes: ID, Name, gene, locus_tag, product, etc.
#'
#' @param df data.frame with attributes column
#' @return data.frame with parsed attributes as columns
#' @keywords internal
parse_gff_attributes <- function(df) {
  # Common attribute keys to extract
  common_attrs <- c("ID", "Name", "gene", "locus_tag", "product", "protein_id", "Parent")

  # Initialize columns
  for (attr in common_attrs) {
    df[[attr]] <- NA_character_
  }

  # Parse each row's attributes
  for (i in seq_len(nrow(df))) {
    attr_string <- df$attributes[i]

    if (is.na(attr_string) || attr_string == ".") {
      next
    }

    # Split by semicolon
    pairs <- strsplit(attr_string, ";", fixed = TRUE)[[1]]

    for (pair in pairs) {
      # Split by first = sign
      kv <- strsplit(pair, "=", fixed = TRUE)[[1]]

      if (length(kv) == 2) {
        key <- trimws(kv[1])
        value <- trimws(kv[2])

        # URL decode value (e.g., %20 -> space)
        value <- utils::URLdecode(value)

        if (key %in% common_attrs) {
          df[[key]][i] <- value
        }
      }
    }
  }

  df
}

#' Extract Sequences by Feature Name (Internal Use Case)
#'
#' @description
#' Internal use case for extracting DNA sequences for features matching a pattern.
#' Searches features by gene name, locus_tag, or product and extracts
#' their sequences from the genome.
#'
#' This is an internal function called by the controller extract_sequences_by_name().
#' Users should use the controller function instead.
#'
#' @param entity A genome_entity object
#' @param pattern Character string or regex pattern to match feature names
#' @param fields Character vector of fields to search (default: c("gene", "locus_tag", "product"))
#' @param options List of options:
#'   - ignore_case: Logical, case-insensitive matching (default TRUE)
#'   - type_filter: Character, filter by feature type (e.g., "CDS", "gene")
#'   - translate: Logical, translate CDS to protein (default FALSE)
#'
#' @return Named character vector of sequences
#' @keywords internal
execute_extract_sequences_by_name <- function(entity, pattern,
                                               fields = c("gene", "locus_tag", "product"),
                                               options = list()) {
  # Validate entity
  validate_genome_entity(entity)

  # Parse options
  ignore_case <- options$ignore_case %||% TRUE
  type_filter <- options$type_filter %||% NULL
  translate <- options$translate %||% FALSE

  # Get features
  features_df <- entity$features

  if (nrow(features_df) == 0) {
    cli::cli_warn("No features in genome_entity")
    return(character())
  }

  # Search for matching features
  matches <- rep(FALSE, nrow(features_df))

  for (field in fields) {
    if (field %in% names(features_df)) {
      field_values <- features_df[[field]]
      field_matches <- grepl(pattern, field_values, ignore.case = ignore_case)
      matches <- matches | field_matches
    }
  }

  # Apply type filter if specified
  if (!is.null(type_filter) && "type" %in% names(features_df)) {
    matches <- matches & (features_df$type == type_filter)
  }

  matching_features <- features_df[matches, ]

  if (nrow(matching_features) == 0) {
    cli::cli_inform("No features matched pattern: {pattern}")
    return(character())
  }

  # Extract sequences for matching features
  sequences <- character(nrow(matching_features))
  names(sequences) <- paste0(
    matching_features$seqname, ":",
    matching_features$start, "-",
    matching_features$end,
    ifelse(!is.na(matching_features$gene), paste0(" ", matching_features$gene), "")
  )

  for (i in seq_len(nrow(matching_features))) {
    feat <- matching_features[i, ]

    # Get sequence for this feature's seqname
    if (!feat$seqname %in% names(entity$sequences$dna_raw)) {
      cli::cli_warn("Seqname '{feat$seqname}' not found in sequences")
      sequences[i] <- NA_character_
      next
    }

    full_seq <- entity$sequences$dna_raw[[feat$seqname]]

    # Extract region
    start_pos <- feat$start
    end_pos <- feat$end

    # Validate coordinates
    if (is.na(start_pos) || is.na(end_pos)) {
      cli::cli_warn("Feature has missing coordinates")
      sequences[i] <- NA_character_
      next
    }

    if (start_pos < 1 || end_pos > nchar(full_seq)) {
      cli::cli_warn("Feature coordinates out of bounds")
      sequences[i] <- NA_character_
      next
    }

    subseq <- substr(full_seq, start_pos, end_pos)

    # Handle strand (reverse complement if on minus strand)
    if (!is.na(feat$strand) && feat$strand == "-") {
      subseq <- reverse_complement(subseq)
    }

    # Translate if requested (only for CDS)
    if (translate && !is.na(feat$type) && feat$type == "CDS") {
      subseq <- translate_dna(subseq, frame = 1, genetic_code = "11", .internal = TRUE)
    }

    sequences[i] <- subseq
  }

  # Remove NA sequences
  sequences <- sequences[!is.na(sequences)]

  sequences
}

# Helper: %||% operator
`%||%` <- function(x, y) if (is.null(x)) y else x

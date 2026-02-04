#' Extract Sequences by Feature Name
#'
#' @description
#' Use case for extracting DNA sequences for features matching a pattern.
#' Searches features by gene name, locus_tag, or product and extracts
#' their sequences from the genome.
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
#' @export
#' @examples
#' \dontrun{
#' # Extract sequences for genes matching "dnaA"
#' seqs <- execute_extract_sequences_by_name(genome, "dnaA")
#'
#' # Extract CDS sequences for locus tags starting with "EC"
#' seqs <- execute_extract_sequences_by_name(
#'   genome,
#'   "^EC",
#'   fields = "locus_tag",
#'   options = list(type_filter = "CDS")
#' )
#' }
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
      subseq <- translate_dna(subseq)
    }

    sequences[i] <- subseq
  }

  # Remove NA sequences
  sequences <- sequences[!is.na(sequences)]

  sequences
}

#' Reverse Complement DNA Sequence
#'
#' @param seq Character string of DNA sequence
#' @return Reverse complement
#' @keywords internal
reverse_complement <- function(seq) {
  # Complement lookup
  complement_map <- c(
    A = "T", T = "A", G = "C", C = "G",
    a = "t", t = "a", g = "c", c = "g",
    N = "N", n = "n",
    R = "Y", Y = "R", S = "S", W = "W", K = "M", M = "K",
    B = "V", V = "B", D = "H", H = "D"
  )

  # Split into characters
  chars <- strsplit(seq, "")[[1]]

  # Complement
  comp_chars <- complement_map[chars]

  # Handle unknown characters
  comp_chars[is.na(comp_chars)] <- "N"

  # Reverse
  paste(rev(comp_chars), collapse = "")
}

#' Translate DNA to Protein
#'
#' @param dna Character string of DNA sequence
#' @return Protein sequence (single-letter amino acids)
#' @keywords internal
translate_dna <- function(dna) {
  # Standard genetic code
  codon_table <- c(
    TTT = "F", TTC = "F", TTA = "L", TTG = "L",
    TCT = "S", TCC = "S", TCA = "S", TCG = "S",
    TAT = "Y", TAC = "Y", TAA = "*", TAG = "*",
    TGT = "C", TGC = "C", TGA = "*", TGG = "W",
    CTT = "L", CTC = "L", CTA = "L", CTG = "L",
    CCT = "P", CCC = "P", CCA = "P", CCG = "P",
    CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
    CGT = "R", CGC = "R", CGA = "R", CGG = "R",
    ATT = "I", ATC = "I", ATA = "I", ATG = "M",
    ACT = "T", ACC = "T", ACA = "T", ACG = "T",
    AAT = "N", AAC = "N", AAA = "K", AAG = "K",
    AGT = "S", AGC = "S", AGA = "R", AGG = "R",
    GTT = "V", GTC = "V", GTA = "V", GTG = "V",
    GCT = "A", GCC = "A", GCA = "A", GCG = "A",
    GAT = "D", GAC = "D", GAA = "E", GAG = "E",
    GGT = "G", GGC = "G", GGA = "G", GGG = "G"
  )

  # Convert to uppercase
  dna <- toupper(dna)

  # Extract codons (groups of 3)
  n_codons <- floor(nchar(dna) / 3)
  if (n_codons == 0) return("")

  codons <- character(n_codons)
  for (i in seq_len(n_codons)) {
    start <- (i - 1) * 3 + 1
    codons[i] <- substr(dna, start, start + 2)
  }

  # Translate
  aa <- codon_table[codons]

  # Handle unknown codons
  aa[is.na(aa)] <- "X"

  paste(aa, collapse = "")
}

# Helper: %||% operator
`%||%` <- function(x, y) if (is.null(x)) y else x

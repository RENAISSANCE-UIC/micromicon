#' Import GFF3 + FASTA Use Case
#'
#' @description
#' Business logic for importing GFF3 annotation and FASTA sequence data
#' into a genome_entity. This use case orchestrates reading from two
#' different gateways and harmonizing the data.
#'
#' Gateway Contracts:
#' - GFF gateway must provide read(path) returning data.frame with columns:
#'   seqname, start, end, strand, type, and optional (gene, locus_tag, product, etc.)
#' - FASTA gateway must provide read(path) returning named character vector
#'   where names are sequence IDs and values are DNA sequences
#'
#' @param gff_gateway Gateway object for GFF3 files
#' @param fasta_gateway Gateway object for FASTA files
#' @param gff_path Path to GFF3 file
#' @param fasta_path Path to FASTA file
#' @param options List of options:
#'   - auto_harmonize: Logical, attempt to harmonize seqname mismatches (default TRUE)
#'   - verbose: Logical, print progress messages (default TRUE)
#'
#' @return A validated genome_entity object
#' @export
#' @examples
#' \dontrun{
#' # With real gateways
#' gff_gw <- create_gff_gateway(get_gff_parser())
#' fasta_gw <- create_fasta_gateway(get_fasta_parser())
#' entity <- execute_import_gff_fasta(gff_gw, fasta_gw, "anno.gff3", "genome.fasta")
#'
#' # With mock gateways (for testing)
#' mock_gff <- list(
#'   read = function(path) {
#'     data.frame(
#'       seqname = "chr1",
#'       start = 1,
#'       end = 100,
#'       strand = "+",
#'       type = "gene",
#'       gene = "geneA",
#'       stringsAsFactors = FALSE
#'     )
#'   }
#' )
#' mock_fasta <- list(
#'   read = function(path) {
#'     c(chr1 = "ATCGATCG")
#'   }
#' )
#' entity <- execute_import_gff_fasta(mock_gff, mock_fasta, "test.gff3", "test.fasta")
#' }
execute_import_gff_fasta <- function(gff_gateway, fasta_gateway,
                                     gff_path, fasta_path,
                                     options = list()) {
  # Parse options
  auto_harmonize <- options$auto_harmonize %||% TRUE
  verbose <- options$verbose %||% TRUE

  inform <- function(...) if (verbose) message(...)

  # Step 1: Read GFF3 via gateway
  inform("Reading GFF3 file...")
  features_df <- gff_gateway$read(gff_path)

  if (!is.data.frame(features_df)) {
    stop("GFF gateway did not return a data.frame", call. = FALSE)
  }

  # Ensure required columns
  required_cols <- c("seqname", "start", "end", "strand", "type")
  missing_cols <- setdiff(required_cols, names(features_df))
  if (length(missing_cols) > 0) {
    stop("GFF data missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # Step 2: Read FASTA via gateway
  inform("Reading FASTA file...")
  dna_raw <- fasta_gateway$read(fasta_path)

  if (!is.character(dna_raw) || is.null(names(dna_raw))) {
    stop("FASTA gateway did not return a named character vector", call. = FALSE)
  }

  # Step 3: Harmonize seqnames if needed
  gff_seqnames <- unique(features_df$seqname)
  fasta_seqnames <- names(dna_raw)

  if (!all(gff_seqnames %in% fasta_seqnames)) {
    if (auto_harmonize) {
      inform("Seqname mismatch detected. Attempting to harmonize...")

      # Try simple prefix removal (e.g., "chr1" -> "1", "lcl|chr1" -> "chr1")
      harmonized <- harmonize_seqnames_simple(gff_seqnames, fasta_seqnames, verbose = verbose)

      if (!is.null(harmonized)) {
        # Apply mapping
        features_df$seqname <- harmonized$mapping[features_df$seqname]
        inform("Seqname harmonization successful")
      } else {
        warning(
          "Could not harmonize seqnames. ",
          "GFF seqnames: ", paste(head(gff_seqnames, 5), collapse = ", "),
          ". FASTA seqnames: ", paste(head(fasta_seqnames, 5), collapse = ", "),
          ". Some features may be orphaned.",
          call. = FALSE
        )
      }
    } else {
      warning(
        "Seqname mismatch detected but auto_harmonize=FALSE. ",
        "Some features may not match sequences.",
        call. = FALSE
      )
    }
  }

  # Step 4: Filter features to only those with matching sequences
  features_df <- features_df[features_df$seqname %in% names(dna_raw), ]

  if (nrow(features_df) == 0) {
    warning("No features matched FASTA sequences after harmonization", call. = FALSE)
  }

  # Step 5: Build metadata from FASTA
  metadata_df <- data.frame(
    seqname = names(dna_raw),
    length_bp = nchar(dna_raw),
    topology = "linear",  # Default assumption
    mol_type = "DNA",     # Default assumption
    stringsAsFactors = FALSE
  )

  # Step 6: Build indices
  indices_list <- execute_build_indices(features_df, metadata_df)

  # Step 7: Create entity
  entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = dna_raw,
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = features_df,
    metadata_df = metadata_df,
    indices_list = indices_list
  )

  # Step 8: Validate
  validate_genome_entity(entity)

  # Step 9: Return
  entity
}

#' Simple Seqname Harmonization
#'
#' @description
#' Attempts to harmonize seqnames between GFF and FASTA by trying common
#' transformations (prefix removal, case changes, prefix addition).
#'
#' @param gff_seqnames Character vector of GFF seqnames
#' @param fasta_seqnames Character vector of FASTA seqnames
#' @param verbose Logical, print messages?
#'
#' @return List with mapping (named vector) or NULL if harmonization failed
#' @keywords internal
harmonize_seqnames_simple <- function(gff_seqnames, fasta_seqnames, verbose = TRUE) {
  inform <- function(...) if (verbose) message(...)

  # Strategy 1: Remove common prefixes from GFF names
  prefixes_to_try <- c("chr", "lcl\\|", "CHR", "Chr")

  for (prefix in prefixes_to_try) {
    cleaned_gff <- sub(paste0("^", prefix), "", gff_seqnames, ignore.case = TRUE)

    if (all(cleaned_gff %in% fasta_seqnames)) {
      mapping <- setNames(cleaned_gff, gff_seqnames)
      inform("  Found mapping by removing prefix '", prefix, "'")
      return(list(mapping = mapping))
    }
  }

  # Strategy 2: Add "chr" prefix to GFF names
  prefixed_gff <- paste0("chr", gff_seqnames)
  if (all(prefixed_gff %in% fasta_seqnames)) {
    mapping <- setNames(prefixed_gff, gff_seqnames)
    inform("  Found mapping by adding 'chr' prefix")
    return(list(mapping = mapping))
  }

  # Strategy 3: Case-insensitive matching
  gff_lower <- tolower(gff_seqnames)
  fasta_lower <- tolower(fasta_seqnames)

  matches <- match(gff_lower, fasta_lower)
  if (all(!is.na(matches))) {
    mapping <- setNames(fasta_seqnames[matches], gff_seqnames)
    inform("  Found mapping by case-insensitive matching")
    return(list(mapping = mapping))
  }

  # Failed
  NULL
}

# Helper: %||% operator (like rlang)
`%||%` <- function(x, y) if (is.null(x)) y else x

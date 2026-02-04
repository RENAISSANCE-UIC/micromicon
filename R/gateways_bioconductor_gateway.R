#' Create Bioconductor Gateway
#'
#' @description
#' Creates a gateway for converting between domain objects (data.frames, character vectors)
#' and Bioconductor S4 objects (GRanges, DNAStringSet). This provides optional
#' Bioconductor integration while keeping the core domain logic framework-independent.
#'
#' Gateway Methods:
#' - is_available(): Check if Bioconductor packages are available
#' - to_granges(features_df): Convert features data.frame to GRanges
#' - from_granges(granges): Convert GRanges to features data.frame
#' - to_dnastringset(sequences_char): Convert character vector to DNAStringSet
#' - from_dnastringset(dna_bio): Convert DNAStringSet to character vector
#'
#' @return Gateway object (list with conversion methods)
#' @export
#' @examples
#' \dontrun{
#' gateway <- create_bioconductor_gateway()
#'
#' if (gateway$is_available()) {
#'   # Convert to Bioconductor objects
#'   granges <- gateway$to_granges(features_df)
#'   dna <- gateway$to_dnastringset(sequences)
#' }
#' }
create_bioconductor_gateway <- function() {
  # Check package availability once
  has_genomic_ranges <- requireNamespace("GenomicRanges", quietly = TRUE)
  has_iranges <- requireNamespace("IRanges", quietly = TRUE)
  has_biostrings <- requireNamespace("Biostrings", quietly = TRUE)
  has_s4vectors <- requireNamespace("S4Vectors", quietly = TRUE)

  is_available <- has_genomic_ranges && has_iranges && has_biostrings && has_s4vectors

  list(
    is_available = function() {
      is_available
    },

    to_granges = function(features_df) {
      if (!is_available) {
        cli::cli_abort(c(
          "Bioconductor packages required for GRanges conversion.",
          "i" = "Install with: BiocManager::install(c('GenomicRanges', 'IRanges', 'S4Vectors'))"
        ))
      }

      if (!is.data.frame(features_df) || nrow(features_df) == 0) {
        # Return empty GRanges
        return(GenomicRanges::GRanges())
      }

      # Ensure required columns exist
      if (!all(c("seqname", "start", "end") %in% names(features_df))) {
        cli::cli_abort("features_df must have columns: seqname, start, end")
      }

      # Filter out features with NA coordinates
      valid_coords <- !is.na(features_df$start) & !is.na(features_df$end) & !is.na(features_df$seqname)

      if (!all(valid_coords)) {
        n_invalid <- sum(!valid_coords)
        cli::cli_warn(c(
          "Removing {n_invalid} feature{?s} with missing coordinates (NA start/end/seqname).",
          "i" = "These features will not be included in the GRanges object."
        ))
        features_df <- features_df[valid_coords, ]
      }

      # Handle empty result after filtering
      if (nrow(features_df) == 0) {
        return(GenomicRanges::GRanges())
      }

      # Build GRanges
      gr <- GenomicRanges::GRanges(
        seqnames = features_df$seqname,
        ranges = IRanges::IRanges(
          start = features_df$start,
          end = features_df$end
        )
      )

      # Add strand if present
      if ("strand" %in% names(features_df)) {
        GenomicRanges::strand(gr) <- features_df$strand
      }

      # Add metadata columns (all columns except coordinate columns and GRanges reserved names)
      coord_cols <- c("seqname", "start", "end", "strand")
      # GRanges reserved slot names that cannot be used as metadata column names
      reserved_names <- c("seqnames", "ranges", "strand", "seqlevels", "seqlengths",
                         "isCircular", "start", "end", "width", "element")
      meta_cols <- setdiff(names(features_df), c(coord_cols, reserved_names))

      if (length(meta_cols) > 0) {
        mcols_df <- features_df[, meta_cols, drop = FALSE]
        S4Vectors::mcols(gr) <- mcols_df
      }

      gr
    },

    from_granges = function(granges) {
      if (!is_available) {
        cli::cli_abort("Bioconductor packages required for GRanges conversion")
      }

      if (!inherits(granges, "GRanges")) {
        cli::cli_abort("Input must be a GRanges object")
      }

      if (length(granges) == 0) {
        # Return empty data.frame
        return(data.frame(
          seqname = character(),
          start = integer(),
          end = integer(),
          strand = character(),
          stringsAsFactors = FALSE
        ))
      }

      # Extract coordinates
      df <- data.frame(
        seqname = as.character(GenomicRanges::seqnames(granges)),
        start = GenomicRanges::start(granges),
        end = GenomicRanges::end(granges),
        strand = as.character(GenomicRanges::strand(granges)),
        stringsAsFactors = FALSE
      )

      # Add metadata columns
      mcols_data <- S4Vectors::mcols(granges)
      if (!is.null(mcols_data) && ncol(mcols_data) > 0) {
        mcols_df <- as.data.frame(mcols_data, stringsAsFactors = FALSE)
        df <- cbind(df, mcols_df)
      }

      df
    },

    to_dnastringset = function(sequences_char) {
      if (!is_available) {
        cli::cli_abort(c(
          "Bioconductor packages required for DNAStringSet conversion.",
          "i" = "Install with: BiocManager::install('Biostrings')"
        ))
      }

      if (!is.character(sequences_char)) {
        cli::cli_abort("Input must be a character vector")
      }

      Biostrings::DNAStringSet(sequences_char)
    },

    from_dnastringset = function(dna_bio) {
      if (!is_available) {
        cli::cli_abort("Bioconductor packages required for DNAStringSet conversion")
      }

      if (!inherits(dna_bio, "DNAStringSet")) {
        cli::cli_abort("Input must be a DNAStringSet object")
      }

      # Convert to character vector
      sequences <- as.character(dna_bio)
      names(sequences) <- names(dna_bio)

      sequences
    },

    # Additional utility: Create indexed FASTA file
    index_fasta = function(fasta_path) {
      if (!requireNamespace("Rsamtools", quietly = TRUE)) {
        cli::cli_abort(c(
          "Rsamtools package required for FASTA indexing.",
          "i" = "Install with: BiocManager::install('Rsamtools')"
        ))
      }

      if (!file.exists(fasta_path)) {
        cli::cli_abort("FASTA file not found: {fasta_path}")
      }

      # Create .fai index
      Rsamtools::indexFa(fasta_path)

      # Return FaFile object
      Rsamtools::FaFile(fasta_path)
    }
  )
}

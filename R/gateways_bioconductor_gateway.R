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
        stop(
          "Bioconductor packages required for GRanges conversion. ",
          "Install with: BiocManager::install(c('GenomicRanges', 'IRanges', 'S4Vectors'))",
          call. = FALSE
        )
      }

      if (!is.data.frame(features_df) || nrow(features_df) == 0) {
        # Return empty GRanges
        return(GenomicRanges::GRanges())
      }

      # Ensure required columns exist
      if (!all(c("seqname", "start", "end") %in% names(features_df))) {
        stop("features_df must have columns: seqname, start, end", call. = FALSE)
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

      # Add metadata columns (all columns except seqname, start, end, strand)
      coord_cols <- c("seqname", "start", "end", "strand")
      meta_cols <- setdiff(names(features_df), coord_cols)

      if (length(meta_cols) > 0) {
        mcols_df <- features_df[, meta_cols, drop = FALSE]
        S4Vectors::mcols(gr) <- mcols_df
      }

      gr
    },

    from_granges = function(granges) {
      if (!is_available) {
        stop(
          "Bioconductor packages required for GRanges conversion",
          call. = FALSE
        )
      }

      if (!inherits(granges, "GRanges")) {
        stop("Input must be a GRanges object", call. = FALSE)
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
        stop(
          "Bioconductor packages required for DNAStringSet conversion. ",
          "Install with: BiocManager::install('Biostrings')",
          call. = FALSE
        )
      }

      if (!is.character(sequences_char)) {
        stop("Input must be a character vector", call. = FALSE)
      }

      Biostrings::DNAStringSet(sequences_char)
    },

    from_dnastringset = function(dna_bio) {
      if (!is_available) {
        stop(
          "Bioconductor packages required for DNAStringSet conversion",
          call. = FALSE
        )
      }

      if (!inherits(dna_bio, "DNAStringSet")) {
        stop("Input must be a DNAStringSet object", call. = FALSE)
      }

      # Convert to character vector
      sequences <- as.character(dna_bio)
      names(sequences) <- names(dna_bio)

      sequences
    },

    # Additional utility: Create indexed FASTA file
    index_fasta = function(fasta_path) {
      if (!requireNamespace("Rsamtools", quietly = TRUE)) {
        stop(
          "Rsamtools package required for FASTA indexing. ",
          "Install with: BiocManager::install('Rsamtools')",
          call. = FALSE
        )
      }

      if (!file.exists(fasta_path)) {
        stop("FASTA file not found: ", fasta_path, call. = FALSE)
      }

      # Create .fai index
      Rsamtools::indexFa(fasta_path)

      # Return FaFile object
      Rsamtools::FaFile(fasta_path)
    }
  )
}

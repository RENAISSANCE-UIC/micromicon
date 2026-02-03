#' Create FASTA Gateway
#'
#' @description
#' Creates a gateway for reading and writing FASTA format files. This gateway
#' provides optional Bioconductor integration for enhanced performance and features.
#'
#' Gateway Contract (for use cases):
#' - read(path): Returns named character vector where names are sequence IDs
#'   and values are DNA sequences
#' - write(sequences, path): Writes sequences to FASTA file
#'
#' @param use_bioconductor Logical; if TRUE and Biostrings is available, use it
#'   for reading/writing. Default TRUE.
#'
#' @return Gateway object (list with read and write methods)
#' @export
#' @examples
#' \dontrun{
#' gateway <- create_fasta_gateway()
#' sequences <- gateway$read("genome.fasta")
#' gateway$write(sequences, "output.fasta")
#' }
create_fasta_gateway <- function(use_bioconductor = TRUE) {
  # Check if Bioconductor is available
  has_biostrings <- use_bioconductor &&
                    requireNamespace("Biostrings", quietly = TRUE)

  list(
    read = function(path) {
      # Validate file exists
      if (!file.exists(path)) {
        stop("FASTA file not found: ", path, call. = FALSE)
      }

      if (has_biostrings) {
        # Use Biostrings for reading (faster, more robust)
        dna <- Biostrings::readDNAStringSet(path)
        # Convert to named character vector
        sequences <- as.character(dna)
        names(sequences) <- names(dna)
      } else {
        # Use simple parser
        sequences <- parse_fasta_simple(path)
      }

      sequences
    },

    write = function(sequences, path, wrap_width = 80) {
      # Validate sequences
      if (!is.character(sequences)) {
        stop("Sequences must be a character vector", call. = FALSE)
      }

      if (has_biostrings) {
        # Use Biostrings for writing
        dna <- Biostrings::DNAStringSet(sequences)
        Biostrings::writeXStringSet(dna, path, width = wrap_width)
      } else {
        # Use simple writer
        write_fasta_simple(sequences, path, wrap_width)
      }

      invisible(path)
    }
  )
}

#' Simple FASTA Parser
#'
#' @description
#' Simple parser for FASTA files when Biostrings is not available.
#' Handles basic FASTA format without advanced features.
#'
#' @param path Path to FASTA file
#' @return Named character vector of sequences
#' @keywords internal
parse_fasta_simple <- function(path) {
  lines <- readLines(path, warn = FALSE)

  # Find header lines (start with >)
  header_idx <- which(grepl("^>", lines))

  if (length(header_idx) == 0) {
    stop("No FASTA headers found in file: ", path, call. = FALSE)
  }

  # Extract sequences
  n_seqs <- length(header_idx)
  sequences <- character(n_seqs)
  names_out <- character(n_seqs)

  for (i in seq_len(n_seqs)) {
    # Extract header (remove > and take first word as ID)
    header_line <- lines[header_idx[i]]
    header_clean <- sub("^>", "", header_line)
    # ID is first word (before space)
    seq_id <- strsplit(header_clean, "\\s+")[[1]][1]
    names_out[i] <- seq_id

    # Extract sequence lines
    start_line <- header_idx[i] + 1
    end_line <- if (i < n_seqs) {
      header_idx[i + 1] - 1
    } else {
      length(lines)
    }

    if (start_line <= end_line) {
      seq_lines <- lines[start_line:end_line]
      # Remove whitespace and combine
      seq_lines <- gsub("\\s", "", seq_lines)
      sequences[i] <- paste(seq_lines, collapse = "")
    } else {
      sequences[i] <- ""
    }
  }

  names(sequences) <- names_out
  sequences
}

#' Simple FASTA Writer
#'
#' @description
#' Simple writer for FASTA files when Biostrings is not available.
#'
#' @param sequences Named character vector of sequences
#' @param path Output file path
#' @param wrap_width Integer; wrap sequences at this width (default 80)
#' @keywords internal
write_fasta_simple <- function(sequences, path, wrap_width = 80) {
  # Open file for writing
  con <- file(path, open = "wt")
  on.exit(close(con))

  for (i in seq_along(sequences)) {
    # Write header
    header <- paste0(">", names(sequences)[i])
    writeLines(header, con)

    # Write sequence (wrapped)
    seq <- sequences[i]
    if (nchar(seq) <= wrap_width) {
      writeLines(seq, con)
    } else {
      # Wrap sequence
      seq_wrapped <- wrap_sequence(seq, wrap_width)
      writeLines(seq_wrapped, con)
    }
  }
}

#' Wrap Sequence to Fixed Width
#'
#' @param seq Character string
#' @param width Integer; wrap width
#' @return Character vector of wrapped lines
#' @keywords internal
wrap_sequence <- function(seq, width) {
  n_chars <- nchar(seq)
  n_lines <- ceiling(n_chars / width)

  lines <- character(n_lines)
  for (i in seq_len(n_lines)) {
    start_pos <- (i - 1) * width + 1
    end_pos <- min(i * width, n_chars)
    lines[i] <- substr(seq, start_pos, end_pos)
  }

  lines
}

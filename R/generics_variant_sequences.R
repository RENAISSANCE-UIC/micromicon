#' Get Reference DNA Sequence for Variants
#'
#' @description
#' Extracts the \code{dna_ref} column from a consequence table produced by
#' \code{\link{annotate_variants}()}, returning a named character vector keyed
#' by gene label.  Only rows with a non-missing reference sequence are returned
#' by default.
#'
#' @param consequence_tbl A data frame from \code{annotate_variants()}.
#' @param gene Optional character string. When supplied, only rows whose
#'   \code{gene} column contains this pattern (partial match, case-insensitive,
#'   strand arrows ignored) are returned. Mutually exclusive with
#'   \code{position}.
#' @param position Numeric scalar or length-2 numeric vector. When a scalar,
#'   keeps only the row whose \code{position} equals that value exactly. When a
#'   length-2 vector \code{c(start, end)}, keeps rows whose \code{position}
#'   falls within that range (inclusive). Mutually exclusive with \code{gene}.
#'   Errors if no matching rows are found. Supply a plain number
#'   (e.g. \code{4438305}); the comma-formatted display value in the table
#'   (e.g. \code{"4,438,305"}) is handled automatically.
#' @param na.rm Logical; if \code{TRUE} (default) rows with \code{NA}
#'   \code{dna_ref} are silently dropped.
#' @param ... Reserved for future arguments.
#'
#' @return A named character vector.  Names are taken from the \code{gene}
#'   column; values are the reference DNA sequences.
#'
#' @seealso \code{\link{get_variant_dna}()}, \code{\link{get_reference_aa}()},
#'   \code{\link{get_variant_aa}()}, \code{\link{annotate_variants}()}
#'
#' @examples
#' \dontrun{
#' ct <- annotate_variants(gd, predict_variants(gd))
#' get_reference_dna(ct)
#' get_reference_dna(ct, gene = "rpoB")
#' get_reference_dna(ct, position = 1834467)
#' }
#'
#' @export
get_reference_dna <- function(consequence_tbl, gene = NULL, position = NULL, na.rm = TRUE, ...) {
  .get_variant_seq_col(consequence_tbl, col = "dna_ref", gene = gene, position = position, na.rm = na.rm)
}


#' Get Variant DNA Sequence for Variants
#'
#' @description
#' Extracts the \code{dna_alt} column from a consequence table produced by
#' \code{\link{annotate_variants}()}.  For coding mutations this is the full
#' mutated CDS; for intergenic mutations it is the flanking window with the
#' mutation applied.  Structural variants return \code{NA}.
#'
#' @inheritParams get_reference_dna
#'
#' @return A named character vector.  Names are taken from the \code{gene}
#'   column; values are the variant DNA sequences.
#'
#' @seealso \code{\link{get_reference_dna}()}, \code{\link{get_reference_aa}()},
#'   \code{\link{get_variant_aa}()}, \code{\link{annotate_variants}()}
#'
#' @examples
#' \dontrun{
#' ct <- annotate_variants(gd, predict_variants(gd))
#' get_variant_dna(ct)
#' get_variant_dna(ct, gene = "rpoB")
#' get_variant_dna(ct, position = 1834467)
#' }
#'
#' @export
get_variant_dna <- function(consequence_tbl, gene = NULL, position = NULL, na.rm = TRUE, ...) {
  .get_variant_seq_col(consequence_tbl, col = "dna_alt", gene = gene, position = position, na.rm = na.rm)
}


#' Get Reference Amino Acid Sequence for Variants
#'
#' @description
#' Extracts the \code{aa_ref} column from a consequence table produced by
#' \code{\link{annotate_variants}()}.  Only available for coding-region
#' mutations; intergenic and structural variants return \code{NA}.
#'
#' @inheritParams get_reference_dna
#'
#' @return A named character vector.  Names are taken from the \code{gene}
#'   column; values are the reference amino acid sequences.
#'
#' @seealso \code{\link{get_reference_dna}()}, \code{\link{get_variant_dna}()},
#'   \code{\link{get_variant_aa}()}, \code{\link{annotate_variants}()}
#'
#' @examples
#' \dontrun{
#' ct <- annotate_variants(gd, predict_variants(gd))
#' get_reference_aa(ct)
#' get_reference_aa(ct, gene = "rpoB")
#' get_reference_aa(ct, position = 1834467)
#' }
#'
#' @export
get_reference_aa <- function(consequence_tbl, gene = NULL, position = NULL, na.rm = TRUE, ...) {
  .get_variant_seq_col(consequence_tbl, col = "aa_ref", gene = gene, position = position, na.rm = na.rm)
}


#' Get Variant Amino Acid Sequence for Variants
#'
#' @description
#' Extracts the \code{aa_alt} column from a consequence table produced by
#' \code{\link{annotate_variants}()}.  For coding mutations this is the full
#' translated variant protein (truncated at the first stop codon for
#' nonsense variants).  Intergenic and structural variants return \code{NA}.
#'
#' @inheritParams get_reference_dna
#'
#' @return A named character vector.  Names are taken from the \code{gene}
#'   column; values are the variant amino acid sequences.
#'
#' @seealso \code{\link{get_reference_dna}()}, \code{\link{get_variant_dna}()},
#'   \code{\link{get_reference_aa}()}, \code{\link{annotate_variants}()}
#'
#' @examples
#' \dontrun{
#' ct <- annotate_variants(gd, predict_variants(gd))
#' get_variant_aa(ct)
#' get_variant_aa(ct, gene = "rpoB")
#' get_variant_aa(ct, position = 1834467)
#' }
#'
#' @export
get_variant_aa <- function(consequence_tbl, gene = NULL, position = NULL, na.rm = TRUE, ...) {
  .get_variant_seq_col(consequence_tbl, col = "aa_alt", gene = gene, position = position, na.rm = na.rm)
}


# -- internal helper ----------------------------------------------------------

#' @keywords internal
.get_variant_seq_col <- function(consequence_tbl, col, gene = NULL, position = NULL, na.rm = TRUE) {
  if (!is.data.frame(consequence_tbl)) {
    cli::cli_abort(c(
      "{.arg consequence_tbl} must be a data frame.",
      "i" = "Supply the output of {.fn annotate_variants}."
    ))
  }
  if (!col %in% names(consequence_tbl)) {
    cli::cli_abort(c(
      "Column {.field {col}} not found in {.arg consequence_tbl}.",
      "i" = "Run {.fn annotate_variants} first to compute variant sequences."
    ))
  }
  if (!"gene" %in% names(consequence_tbl)) {
    cli::cli_abort("{.arg consequence_tbl} must have a {.field gene} column.")
  }

  # --- mutual exclusivity ------------------------------------------------------
  if (!is.null(gene) && !is.null(position)) {
    cli::cli_abort(c(
      "x" = "Supply {.arg gene} or {.arg position}, not both.",
      "i" = "They are alternative ways to identify a locus."
    ))
  }

  tbl <- consequence_tbl

  # --- gene filter -------------------------------------------------------------
  if (!is.null(gene)) {
    gene_clean <- gsub("[\u2190\u2192]", "", tbl$gene)
    gene_clean <- trimws(gene_clean)
    tbl <- tbl[grepl(gene, gene_clean, ignore.case = TRUE), , drop = FALSE]
  }

  # --- position filter ---------------------------------------------------------
  if (!is.null(position)) {
    if (is.character(position)) {
      parsed <- suppressWarnings(as.numeric(gsub(",", "", position)))
      if (!anyNA(parsed) && length(parsed) == 1L) {
        cli::cli_abort(c(
          "x" = "{.arg position} must be numeric, not a character string.",
          "i" = "Did you mean {.code position = {parsed}}?"
        ))
      } else if (!anyNA(parsed) && length(parsed) == 2L) {
        cli::cli_abort(c(
          "x" = "{.arg position} must be numeric, not a character string.",
          "i" = "Did you mean {.code position = c({parsed[1]}, {parsed[2]})}?"
        ))
      } else {
        cli::cli_abort(
          "{.arg position} must be a numeric scalar or a length-2 numeric vector."
        )
      }
    }

    if (!is.numeric(position) && !is.integer(position)) {
      cli::cli_abort(
        "{.arg position} must be a numeric scalar or a length-2 numeric vector."
      )
    }

    pos_numeric <- suppressWarnings(as.integer(gsub(",", "", as.character(tbl$position))))

    if (length(position) == 1L) {
      tbl <- tbl[!is.na(pos_numeric) & pos_numeric == position, , drop = FALSE]
      if (nrow(tbl) == 0L) {
        cli::cli_abort("No variants found at position {position}.")
      }
    } else if (length(position) == 2L) {
      if (position[1] > position[2]) {
        cli::cli_abort(c(
          "x" = "Range {.arg position} must be {.code c(start, end)} with start \u2264 end.",
          "i" = "Did you mean {.code position = c({position[2]}, {position[1]})}?"
        ))
      }
      tbl <- tbl[!is.na(pos_numeric) & pos_numeric >= position[1] & pos_numeric <= position[2], , drop = FALSE]
      if (nrow(tbl) == 0L) {
        cli::cli_abort(
          "No variants found in the range {position[1]}\u2013{position[2]}."
        )
      }
    } else {
      cli::cli_abort(c(
        "x" = "{.arg position} must be a scalar (exact) or length-2 vector (range).",
        "i" = "Got a vector of length {length(position)}."
      ))
    }
  }

  seqs  <- tbl[[col]]
  names <- tbl$gene

  if (na.rm) {
    keep  <- !is.na(seqs)
    seqs  <- seqs[keep]
    names <- names[keep]
  }

  stats::setNames(seqs, names)
}

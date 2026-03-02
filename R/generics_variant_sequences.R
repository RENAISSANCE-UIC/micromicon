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
#'   strand arrows ignored) are returned.
#' @param na.rm Logical; if \code{TRUE} (default) rows with \code{NA}
#'   \code{dna_ref} are silently dropped.
#' @param ... Reserved for future arguments.
#'
#' @return A named character vector.  Names are taken from the \code{gene}
#'   column; values are the reference DNA sequences.
#'
#' @seealso \code{\link{get_alternate_dna}()}, \code{\link{get_reference_aa}()},
#'   \code{\link{get_alternate_aa}()}, \code{\link{annotate_variants}()}
#'
#' @examples
#' \dontrun{
#' ct <- annotate_variants(gd, predict_variants(gd))
#' get_reference_dna(ct)
#' get_reference_dna(ct, gene = "rpoB")
#' }
#'
#' @export
get_reference_dna <- function(consequence_tbl, gene = NULL, na.rm = TRUE, ...) {
  .get_variant_seq_col(consequence_tbl, col = "dna_ref", gene = gene, na.rm = na.rm)
}


#' Get Alternate DNA Sequence for Variants
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
#'   column; values are the alternate DNA sequences.
#'
#' @seealso \code{\link{get_reference_dna}()}, \code{\link{get_reference_aa}()},
#'   \code{\link{get_alternate_aa}()}, \code{\link{annotate_variants}()}
#'
#' @examples
#' \dontrun{
#' ct <- annotate_variants(gd, predict_variants(gd))
#' get_alternate_dna(ct)
#' get_alternate_dna(ct, gene = "rpoB")
#' }
#'
#' @export
get_alternate_dna <- function(consequence_tbl, gene = NULL, na.rm = TRUE, ...) {
  .get_variant_seq_col(consequence_tbl, col = "dna_alt", gene = gene, na.rm = na.rm)
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
#' @seealso \code{\link{get_reference_dna}()}, \code{\link{get_alternate_dna}()},
#'   \code{\link{get_alternate_aa}()}, \code{\link{annotate_variants}()}
#'
#' @examples
#' \dontrun{
#' ct <- annotate_variants(gd, predict_variants(gd))
#' get_reference_aa(ct)
#' get_reference_aa(ct, gene = "rpoB")
#' }
#'
#' @export
get_reference_aa <- function(consequence_tbl, gene = NULL, na.rm = TRUE, ...) {
  .get_variant_seq_col(consequence_tbl, col = "aa_ref", gene = gene, na.rm = na.rm)
}


#' Get Alternate Amino Acid Sequence for Variants
#'
#' @description
#' Extracts the \code{aa_alt} column from a consequence table produced by
#' \code{\link{annotate_variants}()}.  For coding mutations this is the full
#' translated alternate protein (truncated at the first stop codon for
#' nonsense variants).  Intergenic and structural variants return \code{NA}.
#'
#' @inheritParams get_reference_dna
#'
#' @return A named character vector.  Names are taken from the \code{gene}
#'   column; values are the alternate amino acid sequences.
#'
#' @seealso \code{\link{get_reference_dna}()}, \code{\link{get_alternate_dna}()},
#'   \code{\link{get_reference_aa}()}, \code{\link{annotate_variants}()}
#'
#' @examples
#' \dontrun{
#' ct <- annotate_variants(gd, predict_variants(gd))
#' get_alternate_aa(ct)
#' get_alternate_aa(ct, gene = "rpoB")
#' }
#'
#' @export
get_alternate_aa <- function(consequence_tbl, gene = NULL, na.rm = TRUE, ...) {
  .get_variant_seq_col(consequence_tbl, col = "aa_alt", gene = gene, na.rm = na.rm)
}


# ── internal helper ──────────────────────────────────────────────────────────

#' @keywords internal
.get_variant_seq_col <- function(consequence_tbl, col, gene = NULL, na.rm = TRUE) {
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

  tbl <- consequence_tbl

  if (!is.null(gene)) {
    gene_clean <- gsub("[\u2190\u2192]", "", tbl$gene)
    gene_clean <- trimws(gene_clean)
    tbl <- tbl[grepl(gene, gene_clean, ignore.case = TRUE), , drop = FALSE]
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

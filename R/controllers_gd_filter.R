#' Filter a Variants Table
#'
#' @description
#' Convenience function for filtering the output of \code{\link{predict_variants}()}
#' or \code{\link{annotate_variants}()}. Handles display-formatted columns
#' transparently so you do not need to know how frequencies are stored or
#' that gene names carry strand-direction arrows.
#'
#' All arguments are optional and stack (AND logic): only rows satisfying
#' every supplied criterion are returned.
#'
#' @param tbl A data frame returned by \code{predict_variants()} or
#'   \code{annotate_variants()}.
#' @param gene Character scalar; case-insensitive partial match against the
#'   \code{gene} column. Strand arrows (\code{→}/\code{←}) are stripped before
#'   matching, so \code{gene = "mdtA"} matches \code{"mdtA →"} and
#'   \code{gene = "mdt"} matches any \emph{mdt*} gene.
#' @param min_freq Numeric (0–1); retain rows whose allele frequency is at
#'   least this value. The display string (e.g. \code{"85.1\%"}) is converted
#'   automatically. Structural variants reported at 100\% always pass this
#'   filter.
#' @param type Character vector; keep only rows whose \code{type} matches one
#'   of the supplied values (case-insensitive). Common values:
#'   \code{"SNP"}, \code{"DEL"}, \code{"INS"}, \code{"SUB"}, \code{"MOB"},
#'   \code{"AMP"}, \code{"CON"}, \code{"INV"}.
#' @param consequence Character vector; keep only rows whose \code{consequence}
#'   matches one of the supplied values (case-insensitive). Requires
#'   \code{annotate_variants()} to have been run first. Common values:
#'   \code{"synonymous"}, \code{"missense"}, \code{"nonsense"},
#'   \code{"frameshift"}, \code{"inframe_deletion"},
#'   \code{"inframe_insertion"}.
#' @param seq_id Character scalar; keep only rows on this contig or chromosome
#'   (exact match).
#'
#' @return A data frame with the same columns as \code{tbl}, containing only
#'   the rows that satisfy all supplied criteria.
#'
#' @examples
#' \dontrun{
#' variants <- predict_variants(gd)
#'
#' # Any mdt* gene above 80% frequency
#' filter_variants(variants, gene = "mdt", min_freq = 0.80)
#'
#' # Only SNPs and DELs
#' filter_variants(variants, type = c("SNP", "DEL"))
#'
#' # Missense mutations in marR
#' annotated <- annotate_variants(gd, variants)
#' filter_variants(annotated, gene = "marR", consequence = "missense")
#' }
#'
#' @export
filter_variants <- function(tbl,
                            gene        = NULL,
                            min_freq    = NULL,
                            type        = NULL,
                            consequence = NULL,
                            seq_id      = NULL) {
  stopifnot(is.data.frame(tbl))

  keep <- rep(TRUE, nrow(tbl))

  # --- gene: strip strand arrows, then partial case-insensitive match ----------
  if (!is.null(gene)) {
    stopifnot(is.character(gene), length(gene) == 1L)
    gene_clean <- gsub("[\u2190\u2192]", "", tbl$gene)
    gene_clean <- trimws(gene_clean)
    keep <- keep & grepl(gene, gene_clean, ignore.case = TRUE)
  }

  # --- min_freq: parse "85.1%" -> 0.851 ----------------------------------------
  if (!is.null(min_freq)) {
    stopifnot(is.numeric(min_freq), length(min_freq) == 1L,
              min_freq >= 0, min_freq <= 1)
    freq_num <- suppressWarnings(
      as.numeric(sub("%$", "", tbl$freq)) / 100
    )
    keep <- keep & !is.na(freq_num) & (freq_num >= min_freq)
  }

  # --- type: case-insensitive exact match against a set ------------------------
  if (!is.null(type)) {
    stopifnot(is.character(type))
    keep <- keep & (toupper(tbl$type) %in% toupper(type))
  }

  # --- consequence: annotate_variants() output only ----------------------------
  if (!is.null(consequence)) {
    stopifnot(is.character(consequence))
    if (!"consequence" %in% names(tbl)) {
      cli::cli_abort(c(
        "x" = "Column {.field consequence} not found in {.arg tbl}.",
        "i" = "Run {.fn annotate_variants} first to compute molecular consequences."
      ))
    }
    keep <- keep & !is.na(tbl$consequence) &
      (tolower(tbl$consequence) %in% tolower(consequence))
  }

  # --- seq_id: exact match -----------------------------------------------------
  if (!is.null(seq_id)) {
    stopifnot(is.character(seq_id), length(seq_id) == 1L)
    keep <- keep & !is.na(tbl$seq_id) & (tbl$seq_id == seq_id)
  }

  tbl[keep, , drop = FALSE]
}

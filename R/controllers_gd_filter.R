#' Filter a Variants Table
#'
#' @description
#' Convenience function for filtering the output of \code{\link{predict_variants}()}
#' or \code{\link{annotate_variants}()}. Handles display-formatted columns
#' transparently so you do not need to know how frequencies are stored,
#' that positions are comma-formatted, or that gene names carry strand-direction
#' arrows.
#'
#' All arguments are optional and stack (AND logic): only rows satisfying
#' every supplied criterion are returned.
#'
#' @param tbl A data frame or tibble returned by \code{predict_variants()} or
#'   \code{annotate_variants()}.
#' @param gene Character scalar; case-insensitive partial match against the
#'   \code{gene} column. Strand arrows (\code{→}/\code{←}) are stripped before
#'   matching, so \code{gene = "mdtA"} matches \code{"mdtA →"} and
#'   \code{gene = "mdt"} matches any \emph{mdt*} gene. Mutually exclusive
#'   with \code{position}.
#' @param position Numeric scalar or length-2 numeric vector. When a scalar,
#'   keeps only rows whose \code{position} equals that value exactly. When a
#'   length-2 vector \code{c(start, end)}, keeps rows whose \code{position}
#'   falls within that range (inclusive). Mutually exclusive with \code{gene}.
#'   Errors if no matching rows are found. Supply a plain number
#'   (e.g. \code{4438305}); the comma-formatted display value in the table
#'   (e.g. \code{"4,438,305"}) is handled automatically.
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
#'   \code{"synonymous"}, \code{"nonsynonymous"}, \code{"nonsense"},
#'   \code{"frameshift"}, \code{"inframe_deletion"},
#'   \code{"inframe_insertion"}, \code{"inactivating"}.
#' @param contig Character scalar; keep only rows on this contig or chromosome
#'   (exact match).
#'
#' @return A tibble with the same columns as \code{tbl}, containing only
#'   the rows that satisfy all supplied criteria.
#'
#' @examples
#' \dontrun{
#' gd       <- predict_variants(gd)
#' variants <- get_predicted_variants_table(gd)
#'
#' # Any mdt* gene above 80% frequency
#' filter_variants(variants, gene = "mdt", min_freq = 0.80)
#'
#' # Only SNPs and DELs
#' filter_variants(variants, type = c("SNP", "DEL"))
#'
#' # Variant at an exact position
#' filter_variants(variants, position = 1834467)
#'
#' # All variants within a genomic window
#' filter_variants(variants, position = c(1830000, 1840000))
#'
#' # Nonsynonymous mutations in marR (requires annotate_variants first)
#' annotated <- annotate_variants(gd)
#' filter_variants(annotated, gene = "marR", consequence = "nonsynonymous")
#' }
#'
#' @export
filter_variants <- function(tbl,
                            gene        = NULL,
                            position    = NULL,
                            min_freq    = NULL,
                            type        = NULL,
                            consequence = NULL,
                            contig      = NULL) {
  stopifnot(is.data.frame(tbl))

  # --- mutual exclusivity ------------------------------------------------------
  if (!is.null(gene) && !is.null(position)) {
    cli::cli_abort(c(
      "x" = "Supply {.arg gene} or {.arg position}, not both.",
      "i" = "They are alternative ways to identify a locus."
    ))
  }

  keep <- rep(TRUE, nrow(tbl))

  # --- gene: strip strand arrows, then partial case-insensitive match ----------
  if (!is.null(gene)) {
    stopifnot(is.character(gene), length(gene) == 1L)
    gene_clean <- gsub("[\u2190\u2192]", "", tbl$gene)
    gene_clean <- trimws(gene_clean)
    keep <- keep & grepl(gene, gene_clean, ignore.case = TRUE)
  }

  # --- position: exact or range ------------------------------------------------
  if (!is.null(position)) {
    # Friendly error for comma-formatted strings like "1,834,467"
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
      keep <- keep & !is.na(pos_numeric) & (pos_numeric == position)
    } else if (length(position) == 2L) {
      if (position[1] > position[2]) {
        cli::cli_abort(c(
          "x" = "Range {.arg position} must be {.code c(start, end)} with start \u2264 end.",
          "i" = "Did you mean {.code position = c({position[2]}, {position[1]})}?"
        ))
      }
      keep <- keep & !is.na(pos_numeric) &
        (pos_numeric >= position[1]) & (pos_numeric <= position[2])
    } else {
      cli::cli_abort(c(
        "x" = "{.arg position} must be a scalar (exact) or length-2 vector (range).",
        "i" = "Got a vector of length {length(position)}."
      ))
    }

    if (!any(keep)) {
      if (length(position) == 1L) {
        cli::cli_abort("No variants found at position {position}.")
      } else {
        cli::cli_abort(
          "No variants found in the range {position[1]}\u2013{position[2]}."
        )
      }
    }
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

  # --- contig: exact match -----------------------------------------------------
  if (!is.null(contig)) {
    stopifnot(is.character(contig), length(contig) == 1L)
    keep <- keep & !is.na(tbl$seq_id) & (tbl$seq_id == contig)
  }

  tibble::as_tibble(tbl[keep, , drop = FALSE])
}

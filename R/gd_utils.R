
#' Assert length-one logical and return a scalar
#' @keywords internal
gd_as_bool1 <- function(x, what = "value") {
  if (length(x) != 1L || (!isTRUE(is.na(x)) && !is.logical(x))) {
    cli::cli_abort("{.arg {what}} must be a length-one logical")
  }
  as.logical(x)
}


#' Bounds check for 1-indexed string positions
#' @keywords internal
gd_check_pos1 <- function(pos, n, what = "position") {
  if (!is.numeric(pos) || length(pos) != 1L || is.na(pos) || pos < 1L || pos > n) {
    cli::cli_abort("{.arg {what}} out of range [1, {n}] (got {pos})")
  }
  invisible(pos)
}

# ---- DNA utilities ----

#' Reverse-complement a DNA string
#' Falls back to a project-level reverse_complement() if present.
#' @keywords internal
gd_revcomp <- function(dna) {
  if (exists("reverse_complement", mode = "function")) {
    return(reverse_complement(dna))
  }
  if (!is.character(dna) || length(dna) != 1L) {
    cli::cli_abort("{.arg dna} must be a length-one character string")
  }
  s <- strsplit(dna, "", fixed = TRUE)[[1]]
  s <- chartr("ACGTNacgtn", "TGCANtgcan", s)
  paste0(rev(s), collapse = "")
}

#' Complement a base (vectorized; accepts A/C/G/T/N)
#' @keywords internal
gd_complement_base <- function(b) {
  if (!is.character(b)) cli::cli_abort("{.arg b} must be character")
  out <- chartr("ACGTNacgtn", "TGCANtgcan", b)
  bad <- is.na(out) | !nzchar(out)
  if (any(bad)) {
    cli::cli_abort("Invalid base(s) for complement: {paste(unique(b[bad]), collapse = ', ')}")
  }
  toupper(out)
}

#' Replace a single character at a 1-based position
#' @keywords internal
gd_replace_substring1 <- function(s, pos, newchar) {
  if (!is.character(s) || length(s) != 1L) cli::cli_abort("{.arg s} must be length-one character")
  n <- nchar(s, type = "chars", allowNA = FALSE)
  gd_check_pos1(pos, n, "pos")
  pre  <- if (pos > 1L) substr(s, 1L, pos - 1L) else ""
  post <- if (pos < n) substr(s, pos + 1L, n) else ""
  paste0(pre, newchar, post)
}

#' Validate DNA (ACGTN only), returns uppercase no-whitespace
#' @keywords internal
gd_validate_dna <- function(s, allow_n = TRUE) {
  if (!is.character(s) || length(s) != 1L) cli::cli_abort("{.arg s} must be length-one character")
  s2 <- toupper(gsub("\\s+", "", s, perl = TRUE))
  pat <- if (allow_n) "^[ACGTN]*$" else "^[ACGT]*$"
  if (!grepl(pat, s2, perl = TRUE)) cli::cli_abort("Non-DNA characters detected in sequence")
  s2
}


# ---- AA / codon utilities ----

#' Minimal DNA -> AA translation (standard code, '*' for stop, 'X' for ambiguity)
#' Assumes length multiple of 3.
#' @keywords internal
gd_translate_cds <- function(dna) {
  dna <- gd_validate_dna(dna, allow_n = TRUE)
  n <- nchar(dna)
  if (n %% 3L != 0L) cli::cli_abort("CDS length must be divisible by 3 (got {n})")
  codons <- substring(dna, seq(1L, n - 2L, by = 3L), seq(3L, n, by = 3L))
  tbl <- c(
    TTT="F", TTC="F", TTA="L", TTG="L", TCT="S", TCC="S", TCA="S", TCG="S",
    TAT="Y", TAC="Y", TAA="*", TAG="*", TGT="C", TGC="C", TGA="*", TGG="W",
    CTT="L", CTC="L", CTA="L", CTG="L", CCT="P", CCC="P", CCA="P", CCG="P",
    CAT="H", CAC="H", CAA="Q", CAG="Q", CGT="R", CGC="R", CGA="R", CGG="R",
    ATT="I", ATC="I", ATA="I", ATG="M", ACT="T", ACC="T", ACA="T", ACG="T",
    AAT="N", AAC="N", AAA="K", AAG="K", AGT="S", AGC="S", AGA="R", AGG="R",
    GTT="V", GTC="V", GTA="V", GTG="V", GCT="A", GCC="A", GCA="A", GCG="A",
    GAT="D", GAC="D", GAA="E", GAG="E", GGT="G", GGC="G", GGA="G", GGG="G"
  )
  is_amb <- grepl("N", codons, fixed = TRUE)
  aa <- tbl[codons]
  aa[is.na(aa)] <- "X"
  aa[is_amb]    <- "X"
  paste0(unname(aa), collapse = "")
}

#' Codon window and AA index from CDS position (1-based)
#' @keywords internal
gd_codon_window <- function(cds_pos) {
  stopifnot(length(cds_pos) == 1L, is.finite(cds_pos), cds_pos >= 1L)
  codon_start <- ((cds_pos - 1L) %/% 3L) * 3L + 1L
  aa_index    <- ((cds_pos - 1L) %/% 3L) + 1L
  list(codon_start = codon_start, aa_index = aa_index)
}

# ---- light formatting ----

#' Format integer with thousands separator
#' @keywords internal
gd_fmt_int <- function(x) {
  ifelse(is.na(x), NA_character_, formatC(as.integer(x), big.mark = ",", format = "d"))
}

#' Format frequency (0 to 1) as percentage with 1 decimal
#' @keywords internal
gd_fmt_pct <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.1f%%", 100 * as.numeric(x)))
}

#' Normalize arrows / slashes spacing in gene/product fields
#' @keywords internal
gd_arrowize <- function(x) {
  if (is.null(x)) return(x)
  x <- gsub("\\s+", " ", x, perl = TRUE)
  x <- gsub("/", " / ", x, perl = TRUE)
  trimws(x)
}



# # strand-aware base to apply in CDS space
# new_base_cds <- if (strand == "+") new_base else gd_complement_base(new_base)
# 
# # splice into CDS
# gene_dna_mut <- gd_replace_substring1(gene_dna_ref, cds_pos, new_base_cds)
# 
# # codon/AA bookkeeping
# w <- gd_codon_window(cds_pos)
# codon_ref <- substr(gene_dna_ref, w$codon_start, w$codon_start + 2L)
# codon_new <- substr(gene_dna_mut, w$codon_start, w$codon_start + 2L)
# gene_aa_mut <- gd_translate_cds(gene_dna_mut)
# aa_ref_i <- substr(gene_aa_ref, w$aa_index, w$aa_index)
# aa_new_i <- substr(gene_aa_mut, w$aa_index, w$aa_index)

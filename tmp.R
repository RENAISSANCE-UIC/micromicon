

# --- 0) Small local helpers (session-only; no package edits) ------------------

.translate_cds_dna <- function(dna) {
  # minimal, standard code; '*' for stop; 'X' for ambiguity
  dna <- toupper(gsub("\\s+", "", dna, perl = TRUE))
  n <- nchar(dna)
  if (n %% 3L != 0L) stop("CDS length must be divisible by 3; got ", n, call. = FALSE)
  codons <- substring(dna, seq(1L, n-2L, by = 3L), seq(3L, n, by = 3L))
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
  amb <- grepl("N", codons, fixed = TRUE)
  aa <- tbl[codons]; aa[is.na(aa)] <- "X"; aa[amb] <- "X"
  paste0(unname(aa), collapse = "")
}

.cds_window <- function(cds_pos) {
  stopifnot(length(cds_pos) == 1L, is.finite(cds_pos), cds_pos >= 1L)
  codon_start <- ((cds_pos - 1L) %/% 3L) * 3L + 1L
  aa_index    <- ((cds_pos - 1L) %/% 3L) + 1L
  list(codon_start = codon_start, aa_index = aa_index)
}

# --- 1) Find the RA evidence at 1:1,834,467 ----------------------------------

ev_tbl <- gd_events_table(gd, kinds = "evidence", expand_tags = FALSE)
ra <- subset(ev_tbl, type == "RA" & seq_id == 1L & position == 1834467L)

if (nrow(ra) != 1L) stop("Expected exactly one RA at 1:1,834,467; found ", nrow(ra), call. = FALSE)

# Keep declared bases (RA carries these explicitly)
ref_base <- ra$ra_ref_base
alt_base <- ra$ra_new_base

# Optional: corroborate via tags if present (parser preserved them)
# (These helpers exist per the handoff; safe to call and ignore if missing.)
ra <- try(tag_get_first(ra, "aa_ref_seq"), silent = TRUE)
ra <- try(tag_get_first(ra, "aa_new_seq"), silent = TRUE)
ra <- try(tag_get_first(ra, "codon_ref_seq"), silent = TRUE)
ra <- try(tag_get_first(ra, "codon_new_seq"), silent = TRUE)
ra <- try(tag_get_first(ra, "codon_number"),  silent = TRUE)

# --- 2) Confirm the target gene and strand -----------------------------------
search_features(gd, type = "gene")
feat <- search_features(gd, type = "CDS", pattern = "^acrB_1$")
if (!nrow(feat)) stop("Gene 'acrB_1' not found")
if (nrow(feat) > 1L) cli::cli_warn("Multiple 'acrB_1' features; using the first.")
strand <- feat$strand[1L]  # "+" or "-"

# --- 3) Map genomic position -> CDS position ---------------------------------

cds_map <- map_genomic_to_cds(gd, "acrB_1", 1834467L)
cds_pos <- if (is.list(cds_map)) cds_map$cds_pos %||% cds_map$cds_position else cds_map
stopifnot(length(cds_pos) == 1L, is.finite(cds_pos), !is.na(cds_pos))

# --- 4) Pull reference gene DNA/AA, splice the SNP in CDS space --------------

gene_dna_ref <- get_gene_dna(gd, "acrB_1")
gene_aa_ref  <- get_gene_aa(gd,  "acrB_1")

# Sanity: does the CDS base match the expected ref (accounting for strand)?
ref_base_cds <- substr(gene_dna_ref, cds_pos, cds_pos)
ref_expect   <- if (strand == "+") ref_base else gd_complement_base(ref_base)
if (!identical(toupper(ref_base_cds), toupper(ref_expect))) {
  cli::cli_warn(c(
    "!" = "CDS base @ {cds_pos} is {ref_base_cds}, but RA ref is {ref_base} on {strand} strand.",
    " " = "Coordinates may be offset in annotation; proceed but verify."
  ))
}

# Apply the alternate base in CDS coordinates (complement if minus strand)
alt_cds      <- if (strand == "+") alt_base else gd_complement_base(alt_base)
gene_dna_mut <- gd_replace_substring1(gene_dna_ref, cds_pos, alt_cds)

# --- 5) Codon/AA bookkeeping + translation -----------------------------------

w <- .cds_window(cds_pos)
codon_ref <- substr(gene_dna_ref, w$codon_start, w$codon_start + 2L)
codon_new <- substr(gene_dna_mut, w$codon_start, w$codon_start + 2L)

gene_aa_mut <- .translate_cds_dna(gene_dna_mut)
aa_ref_i    <- substr(gene_aa_ref, w$aa_index, w$aa_index)
aa_new_i    <- substr(gene_aa_mut, w$aa_index, w$aa_index)

cli::cli_inform(c(
  "v" = "acrB_1 {strand}-strand | genomic 1:{formatC(1834467L, big.mark=',')} {ref_base}->{alt_base}",
  "v" = "CDS pos {w$aa_index * 3 - (3 - ((cds_pos - 1L) %% 3L) - 1L)} (in codon starting {w$codon_start})",
  "v" = "Codon {toupper(codon_ref)} -> {toupper(codon_new)}",
  "v" = "AA {aa_ref_i}{w$aa_index}{aa_new_i}"
))

# Optional: RA frequency if present
if (!is.null(ra$ev_frequency)) {
  cli::cli_inform(c("i" = "Evidence frequency ≈ {sprintf('%.1f%%', 100*as.numeric(ra$ev_frequency))}"))
}


# File: gateways_reference.R  (Mode B)
# API: create_reference_gateway(entity, logger = NULL)
# Exposes: nt_window, translate, revcomp
# Zero collisions with Mode A.


# Reference gateway for genome_entity objects used by gd utilities.
# - Uses official get_roi_dna(x, chrom, start, end, strand = "+")
# - translate_dna() is required and remains literal/start-agnostic
# - reverse complement: local fallback, auto-switch to exported fn if you add it

gd_create_reference_gateway <- function(entity, logger = NULL) {
  stopifnot(inherits(entity, "genome_entity"))
  
  has_translate <- exists("translate_dna", mode = "function")
  # Prefer a future exported function named reverse_complement()
  has_revcomp_exported <- exists("reverse_complement", mode = "function")
  
  nt_window <- function(seq_id, start, end, strict = TRUE) {
    if (start > end) stop("nt_window: start > end")
    get_roi_dna(entity, seq_id, start, end, strand = "+")
  }
  
  translate <- function(dna, transl_table = 11L) {
    if (!has_translate)
      stop("translate_dna() not found; required by gd utilities.")
    translate_dna(dna, codon_table = transl_table)
  }
  
  # Local fallback (safe and small) — auto-uses exported reverse_complement() when you export it
  .compl_map <- function(x) chartr("ACGTacgtnN", "TGCAtgcanN", x)
  .revcomp_local <- function(s) {
    ssplit <- strsplit(s, "", fixed = TRUE)[[1]]
    paste0(rev(.compl_map(ssplit)), collapse = "")
  }
  revcomp <- function(dna) {
    if (has_revcomp_exported) return(get("reverse_complement")(dna))
    .revcomp_local(dna)
  }
  
  list(
    nt_window = nt_window,
    translate = translate,
    revcomp   = revcomp
  )
}




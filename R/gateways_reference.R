
# File: gateways_reference.R  (Mode B)
# API: create_reference_gateway(entity, logger = NULL)
# Exposes: nt_window, translate, revcomp
# Zero collisions with Mode A.


gd_create_reference_gateway <- function(entity, logger = NULL) {
  stopifnot(inherits(entity, "genome_entity"))
  
  has_entity_roi <- isTRUE(!is.null(entity$get_roi_dna)) && is.function(entity$get_roi_dna)
  has_global_roi <- exists("get_roi_dna", mode = "function")
  has_rc_entity  <- isTRUE(!is.null(entity$revcomp)) && is.function(entity$revcomp)
  has_rc_global  <- exists("revcomp", mode = "function")
  has_translate  <- exists("translate_dna", mode = "function")
  
  seq_lookup <- NULL
  if (!has_entity_roi && !has_global_roi && !is.null(entity$sequences) && is.list(entity$sequences)) {
    seq_lookup <- entity$sequences
  }
  
  .compl_map <- function(x) chartr("ACGTacgtnN", "TGCAtgcanN", x)
  .revcomp_local <- function(s) {
    ssplit <- strsplit(s, "", fixed = TRUE)[[1]]
    paste0(rev(.compl_map(ssplit)), collapse = "")
  }
  
  nt_window <- function(seq_id, start, end, strict = TRUE) {
    if (start > end) stop("nt_window: start > end")
    if (has_entity_roi) {
      return(entity$get_roi_dna(seq_id, start, end, strand = "+"))
    } else if (has_global_roi) {
      return(get("get_roi_dna")(seq_id, start, end, strand = "+"))
    } else if (!is.null(seq_lookup) && !is.null(seq_lookup[[seq_id]])) {
      s <- seq_lookup[[seq_id]]
      if (strict && (start < 1L || end > nchar(s))) {
        stop(sprintf("nt_window out of bounds for %s: [%d,%d] (len=%d)", seq_id, start, end, nchar(s)))
      }
      substr(s, start, end)
    } else {
      stop("nt_window: supply get_roi_dna() or entity$sequences[[seq_id]].")
    }
  }
  
  translate <- function(dna, transl_table = 11L) {
    if (!has_translate) stop("translate_dna() not found; required by gd utilities.")
    translate_dna(dna, codon_table = transl_table)  # literal translator per your rule
  }
  
  revcomp <- function(dna) {
    if (has_rc_entity) return(entity$revcomp(dna))
    if (has_rc_global) return(get("revcomp")(dna))
    .revcomp_local(dna)
  }
  
  list(
    nt_window = nt_window,
    translate = translate,
    revcomp   = revcomp
  )
}



gd_create_reference_gateway <- function(entity, logger = NULL) {
  stopifnot(inherits(entity, "genome_entity"))
  has_translate  <- exists("translate_dna", mode = "function")
  has_revcomp_ex <- exists("reverse_complement", mode = "function")
  
  nt_window <- function(seq_id, start, end, strict = TRUE) {
    if (start > end) cli::cli_abort("nt_window: start > end")
    get_roi_dna(entity, seq_id, start, end, strand = "+")
  }
  
  translate <- function(dna, transl_table = 11L) {
    gd_translate_literal(dna, transl_table = transl_table, frame = 1L)
  }
  
  .compl_map <- function(x) chartr("ACGTacgtnN", "TGCAtgcanN", x)
  .revcomp_local <- function(s) {
    ssplit <- strsplit(s, "", fixed = TRUE)[[1]]
    paste0(rev(.compl_map(ssplit)), collapse = "")
  }
  revcomp <- function(dna) {
    if (has_revcomp_ex) return(get("reverse_complement")(dna))
    .revcomp_local(dna)
  }
  
  list(nt_window = nt_window, translate = translate, revcomp = revcomp)
}

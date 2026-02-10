
gd_create_annotation_gateway <- function(entity, logger = NULL) {
  stopifnot(inherits(entity, "genome_entity"))
  feat <- entity$features
  
  cds_by_position <- function(seq_id, pos) {
    feat[feat$seqname == seq_id &
           feat$type == "CDS" &
           feat$start <= pos & feat$end >= pos, , drop = FALSE]
  }
  
  codon_context <- function(ref_gw, cds_row, pos) {
    strand <- as.character(cds_row$strand)[1]
    start  <- as.integer(cds_row$start[1])
    end    <- as.integer(cds_row$end[1])
    seq_id <- as.character(cds_row$seqname[1])
    tt     <- suppressWarnings(as.integer(cds_row$transl_table[1])); if (is.na(tt)) tt <- 11L
    
    if (strand == "+") {
      idx0 <- pos - start
      codon_index <- floor(idx0 / 3L) + 1L
      codon_pos   <- (idx0 %% 3L) + 1L
      codon_start <- start + 3L * (codon_index - 1L)
      codon_end   <- codon_start + 2L
      codon_seq   <- ref_gw$nt_window(seq_id, codon_start, codon_end, strict = TRUE)
      aa_ref      <- ref_gw$translate(codon_seq, transl_table = tt)
    } else if (strand == "-") {
      idx0 <- end - pos
      codon_index <- floor(idx0 / 3L) + 1L
      codon_pos   <- (idx0 %% 3L) + 1L
      codon_end   <- end - 3L * (codon_index - 1L)
      codon_start <- codon_end - 2L
      codon_seq_fw <- ref_gw$nt_window(seq_id, codon_start, codon_end, strict = TRUE)
      codon_seq <- ref_gw$revcomp(codon_seq_fw)
      aa_ref    <- ref_gw$translate(codon_seq, transl_table = tt)
    } else {
      cli::cli_abort("Unknown strand value: {.val {strand}}")
    }
    
    list(
      seq_id = seq_id, strand = strand, transl_table = tt,
      codon_index = codon_index, codon_pos = codon_pos,
      codon_start = codon_start, codon_end = codon_end,
      codon_seq_ref = codon_seq, aa_ref = aa_ref
    )
  }
  
  list(
    cds_by_position = cds_by_position,
    codon_context   = codon_context
  )
}

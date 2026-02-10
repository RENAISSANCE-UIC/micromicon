gd_create_annotation_gateway <- function(entity, logger = NULL) {
  stopifnot(inherits(entity, "genome_entity"))
  feat <- entity$features
  
  cds_by_position <- function(seq_id, pos) {
    gd_filter_cds_covering(feat, seq_id, pos)
  }
  
  codon_context <- function(ref_gw, cds_row, pos) {
    strand <- as.character(cds_row$strand[1])
    start  <- gd_parse_int(cds_row$start[1])
    end    <- gd_parse_int(cds_row$end[1])
    seq_id <- as.character(cds_row$seqname[1])
    tt     <- gd_get_transl_table(cds_row, default = 11L)
    
    # Optional heads-up for multi-segment (we still approximate contiguity here)
    if (gd_col_has(cds_row, "location_type") && !identical(as.character(cds_row$location_type[1]), "single")) {
      cli::cli_warn("Multi-segment CDS at {seq_id}:{start}-{end}; using contiguous approximation for codon mapping.")
    }
    
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
  
  list(cds_by_position = cds_by_position, codon_context = codon_context)
}

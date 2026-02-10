

gd_mb_coalesce <- function(...) {
  args <- list(...)
  for (a in args) if (!is.null(a) && length(a) > 0 && !is.na(a)[1]) return(a)
  NA_character_
}

gd_complement_base <- function(b) {
  switch(toupper(b),
         "A" = "T", "T" = "A", "G" = "C", "C" = "G", "N" = "N",
         stop(sprintf("Invalid base for complement: %s", b))
  )
}

gd_replace_substring1 <- function(s, pos, newchar) {
  if (pos < 1L || pos > nchar(s)) stop("gd_replace_substring1: pos out of range")
  paste0(substr(s, 1L, pos - 1L), newchar, substr(s, pos + 1L, nchar(s)))
}


gd_summarize_snp <- function(gd_row, entity, gw_annot = NULL, 
                             gw_ref = NULL, logger = NULL) {
  stopifnot(inherits(entity, "genome_entity"))
  if (is.null(gw_annot)) gw_annot <- gd_create_annotation_gateway(entity, logger = logger)
  if (is.null(gw_ref))    gw_ref   <- gd_create_reference_gateway(entity,   logger = logger)
  
  seq_id <- as.character(gd_row$chrom %||% gd_row$seq_id %||% gd_row$contig)
  pos    <- as.integer(gd_row$pos)
  alt    <- as.character(gd_row$alt)
  ref    <- as.character(gd_row$ref)
  
  # Locate CDS
  cds_rows <- gw_annot$cds_by_position(seq_id, pos)
  if (nrow(cds_rows) == 0L) {
    return(data.frame(
      id = gd_row$id, type = gd_row$type, chrom = seq_id, pos = pos,
      gene = NA_character_, strand = NA_character_,
      dna_before = NA_character_, dna_after = NA_character_,
      aa_before  = NA_character_,  aa_after  = NA_character_,
      effect = "intergenic",
      ref = ref, alt = alt, length = 1L,
      cds_pos = NA_integer_, aa_pos = NA_integer_,
      messages = I(list(character()))
    ))
  }
  
  # Multiple CDS overlap — warn, pick first deterministically
  msg <- character()
  if (nrow(cds_rows) > 1L) {
    msg <- c(msg, sprintf("Multiple CDS overlap %s:%d; using first.", seq_id, pos))
  }
  row1   <- cds_rows[1, , drop = FALSE]
  ctx    <- gw_annot$codon_context(gw_ref, row1, pos)
  strand <- ctx$strand
  tt     <- ctx$transl_table
  
  # Validate reference base if provided
  ref_at_pos <- gw_ref$nt_window(seq_id, pos, pos, strict = TRUE)
  if (!is.na(ref) && nzchar(ref) && toupper(ref) != toupper(ref_at_pos)) {
    stop(sprintf("Ref base '%s' does not match reference '%s' at %s:%d", ref, ref_at_pos, seq_id, pos))
  }
  
  # Genomic alt -> gene-orientation alt if antisense
  alt_gene <- if (strand == "+") toupper(alt) else gd_complement_base(alt)
  
  codon_alt     <- gd_replace_substring1(ctx$codon_seq_ref, ctx$codon_pos, alt_gene)
  aa_ref_local  <- ctx$aa_ref
  aa_alt_local  <- gw_ref$translate(codon_alt, transl_table = tt)
  
  effect <- if (aa_ref_local == aa_alt_local) {
    "synonymous"
  } else if (aa_alt_local == "*") {
    "nonsense"
  } else {
    "missense"
  }
  
  gene_name <- gd_mb_coalesce(row1$gene, row1$locus_tag, row1$Name, row1$ID)
  
  # cds_pos: base index in CDS, aa_pos: codon index
  cds_pos <- {
    if (strand == "+") ctx$codon_start - as.integer(row1$start[1]) + ctx$codon_pos
    else               as.integer(row1$end[1]) - ctx$codon_end + ctx$codon_pos
  }
  aa_pos  <- ctx$codon_index
  
  data.frame(
    id = gd_row$id, type = gd_row$type, chrom = seq_id, pos = pos,
    gene = gene_name, strand = strand,
    dna_before = ctx$codon_seq_ref, dna_after = codon_alt,
    aa_before  = aa_ref_local,      aa_after  = aa_alt_local,
    effect = effect,
    ref = ref_at_pos, alt = alt, length = 1L,
    cds_pos = cds_pos, aa_pos = aa_pos,
    messages = I(list(msg))
  )
}


gd_summarize_effects <- function(gd, entity, logger = NULL) {
  stopifnot(inherits(entity, "genome_entity"),
            is.data.frame(gd))  # ok whether or not class 'genome_entity_gd' is set
  
  gw_annot <- gd_create_annotation_gateway(entity, logger)
  gw_ref   <- gd_create_reference_gateway(entity,   logger)
  
  out <- lapply(seq_len(nrow(gd)), function(i) {
    row  <- gd[i, , drop = FALSE]
    type <- as.character(row$type)
    if (identical(type, "SNP")) {
      gd_summarize_snp(row, entity, gw_annot, gw_ref, logger)
    } else {
      data.frame(
        id = row$id, type = type, chrom = row$chrom, pos = row$pos,
        gene = NA_character_, strand = NA_character_,
        dna_before = NA_character_, dna_after = NA_character_,
        aa_before  = NA_character_, aa_after  = NA_character_,
        effect = "undetermined",
        ref = row$ref, alt = row$alt, length = row$length %||% NA_integer_,
        cds_pos = NA_integer_, aa_pos = NA_integer_,
        messages = I(list(sprintf("Type '%s' not yet implemented", type)))
      )
    }
  })
  do.call(rbind, out)
}


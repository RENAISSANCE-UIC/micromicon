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
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
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
  
  seq_id <- as.character(gd_row$chrom %||% gd_row$seq_id %||% gd_row$contig)
  pos    <- as.integer(gd_row$pos)
  alt    <- as.character(gd_row$alt)
  ref    <- as.character(gd_row$ref)
  
  cds_rows <- gw_annot$cds_by_position(seq_id, pos)
  if (nrow(cds_rows) == 0L) {
    return(data.frame(
      id = gd_row$id, type = gd_row$type, chrom = seq_id, pos = pos,
      gene = NA_character_, strand = NA_character_,
      dna_before = NA_character_, dna_after = NA_character_,
      aa_before  = NA_character_, aa_after  = NA_character_,
      effect = "intergenic",
      ref = ref, alt = alt, length = 1L,
      cds_pos = NA_integer_, aa_pos = NA_integer_,
      messages = I(list(character()))
    ))
  }
  
  msg <- character()
  if (nrow(cds_rows) > 1L) {
    msg <- c(msg, sprintf("Multiple CDS overlap %s:%d; using first.", seq_id, pos))
  }
  row1   <- cds_rows[1, , drop = FALSE]
  ctx    <- gw_annot$codon_context(gw_ref, row1, pos)
  strand <- ctx$strand
  tt     <- ctx$transl_table
  
  ref_at_pos <- gw_ref$nt_window(seq_id, pos, pos, strict = TRUE)
  if (!is.na(ref) && nzchar(ref) && toupper(ref) != toupper(ref_at_pos)) {
    stop(sprintf("Ref base '%s' does not match reference '%s' at %s:%d", ref, ref_at_pos, seq_id, pos))
  }
  
  alt_gene     <- if (strand == "+") toupper(alt) else gd_complement_base(alt)
  codon_alt    <- gd_replace_substring1(ctx$codon_seq_ref, ctx$codon_pos, alt_gene)
  aa_ref_local <- ctx$aa_ref
  aa_alt_local <- gw_ref$translate(codon_alt, transl_table = tt)
  
  effect <- if (aa_ref_local == aa_alt_local) "synonymous"
  else if (aa_alt_local == "*")     "nonsense"
  else                               "missense"
  
  gene_id <- gd_resolve_gene_id(row1)
  
  cds_pos <- {
    start <- as.integer(row1$start[1]); end <- as.integer(row1$end[1])
    if (strand == "+") ctx$codon_start - start + ctx$codon_pos
    else               end - ctx$codon_end + ctx$codon_pos
  }
  aa_pos <- ctx$codon_index
  
  data.frame(
    id = gd_row$id, type = gd_row$type, chrom = seq_id, pos = pos,
    gene = gene_id, strand = strand,
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

gd_get_transl_table <- function(cds_row, default = 11L) {
  if (gd_col_has(cds_row, "transl_table")) {
    tt <- suppressWarnings(as.integer(cds_row$transl_table[1]))
    if (!is.na(tt)) return(tt)
  }
  qv <- suppressWarnings(as.integer(gd_get_qual(cds_row, "transl_table", NA_character_)))
  if (!is.na(qv)) return(qv)
  as.integer(default)
}

gd_resolve_gene_id <- function(cds_row) {
  v <- if (gd_col_has(cds_row, "locus_tag")) gd_get1(cds_row, "locus_tag") else gd_get_qual(cds_row, "locus_tag")
  if (!is.na(v) && nzchar(v)) return(v)
  v <- if (gd_col_has(cds_row, "gene"))      gd_get1(cds_row, "gene")      else gd_get_qual(cds_row, "gene")
  if (!is.na(v) && nzchar(v)) return(v)
  if (gd_col_has(cds_row, "ID"))   { v <- gd_get1(cds_row, "ID");   if (!is.na(v) && nzchar(v)) return(v) }
  if (gd_col_has(cds_row, "Name")) { v <- gd_get1(cds_row, "Name"); if (!is.na(v) && nzchar(v)) return(v) }
  seq_id <- gd_get1(cds_row, "seqname"); s <- gd_parse_int(cds_row$start[1]); e <- gd_parse_int(cds_row$end[1])
  paste0(seq_id, ":", s, "-", e)
}

gd_parse_int <- function(x) {
  # x is character or factor; always coerce to character first
  xs <- as.character(x)
  
  # remove all non-numeric characters (except minus sign just in case)
  xs_clean <- gsub("[^0-9-]", "", xs)
  
  # Convert to integer; return NA if empty or non-convertible
  out <- suppressWarnings(as.integer(xs_clean))
  out
}

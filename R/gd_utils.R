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

gd_translate_literal <- function(dna, transl_table = 11L, frame = 1L) {
  translate_dna(dna, frame = frame, genetic_code = as.character(transl_table), .internal = TRUE)
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
  xs <- as.character(x)
  xs_clean <- gsub("[^0-9-]", "", xs)
  out <- suppressWarnings(as.integer(xs_clean))
  out
}

gd_nt_window <- function(entity, seq_id, start, end, strict = TRUE) {
  if (start > end) stop("gd_nt_window: start > end")
  get_roi_dna(entity, seq_id, start, end, strand = "+")
}

gd_translate_literal <- function(dna, transl_table = 11L, frame = 1L) {
  # Your translate_dna signature uses 'genetic_code'
  translate_dna(dna, frame = frame, genetic_code = as.character(transl_table), .internal = TRUE)
}

gd_revcomp <- function(dna) {
  if (exists("reverse_complement", mode = "function")) return(reverse_complement(dna))
  map <- function(x) chartr("ACGTacgtnN", "TGCAtgcanN", x)
  s   <- strsplit(dna, "", fixed = TRUE)[[1]]
  paste0(rev(map(s)), collapse = "")
}

gd_pos_in_row <- function(cds_row, pos) {
  # Fast path: single interval
  if (!gd_col_has(cds_row, "location_type") || identical(as.character(cds_row$location_type[1]), "single")) {
    s <- gd_parse_int(cds_row$start[1]); e <- gd_parse_int(cds_row$end[1])
    return(!is.na(s) && !is.na(e) && pos >= s && pos <= e)
  }
  
  # Multi-segment via 'ranges' (PGAP)
  if (gd_col_has(cds_row, "ranges")) {
    rr <- cds_row$ranges[[1]]
    if (is.data.frame(rr) || is.matrix(rr)) {
      starts <- gd_parse_int(rr[, 1]); ends <- gd_parse_int(rr[, 2])
      hit <- (pos >= starts) & (pos <= ends)
      return(any(hit[!is.na(hit)]))
    } else if (is.list(rr)) {
      ok <- vapply(rr, function(x) {
        if (length(x) < 2) return(FALSE)
        a <- gd_parse_int(x[[1]]); b <- gd_parse_int(x[[2]])
        !is.na(a) && !is.na(b) && pos >= a && pos <= b
      }, logical(1))
      return(any(ok))
    }
  }
  
  # Fallback: parse 'location_string'
  if (gd_col_has(cds_row, "location_string")) {
    loc <- cds_row$location_string[1]
    loc <- gsub("complement\\(|join\\(|order\\(|\\)", "", loc)
    tokens <- unlist(strsplit(loc, ",", fixed = TRUE))
    ok <- vapply(tokens, function(tok) {
      m <- regexec("([0-9<>?]+)\\.\\.([0-9<>?]+)", tok)
      r <- regmatches(tok, m)[[1]]
      if (length(r) != 3) return(FALSE)
      a <- gd_parse_int(r[2]); b <- gd_parse_int(r[3])
      !is.na(a) && !is.na(b) && pos >= a && pos <= b
    }, logical(1))
    return(any(ok))
  }
  
  # Last resort
  s <- gd_parse_int(cds_row$start[1]); e <- gd_parse_int(cds_row$end[1])
  !is.na(s) && !is.na(e) && pos >= s && pos <= e
}

gd_cds_by_position <- function(entity, seq_id, pos) {
  feat <- entity$features
  cand <- feat[feat$seqname == seq_id & feat$type == "CDS", , drop = FALSE]
  keep <- vapply(seq_len(nrow(cand)), function(i) gd_pos_in_row(cand[i, , drop = FALSE], pos), logical(1))
  cand[keep, , drop = FALSE]
}

gd_codon_context <- function(entity, cds_row, pos) {
  strand <- as.character(cds_row$strand[1])
  start  <- gd_parse_int(cds_row$start[1])
  end    <- gd_parse_int(cds_row$end[1])
  seq_id <- as.character(cds_row$seqname[1])
  tt     <- gd_get_transl_table(cds_row, default = 11L)
  
  # Heads-up if PGAP multi-segment; still using contiguous approximation
  if (gd_col_has(cds_row, "location_type") && !identical(as.character(cds_row$location_type[1]), "single")) {
    cli::cli_warn("Multi-segment CDS at {seq_id}:{start}-{end}; using contiguous approximation for codon mapping.")
  }
  
  if (strand == "+") {
    idx0 <- pos - start
    codon_index <- floor(idx0 / 3L) + 1L
    codon_pos   <- (idx0 %% 3L) + 1L
    codon_start <- start + 3L * (codon_index - 1L)
    codon_end   <- codon_start + 2L
    codon_seq   <- gd_nt_window(entity, seq_id, codon_start, codon_end, strict = TRUE)
    aa_ref      <- gd_translate_literal(codon_seq, transl_table = tt, frame = 1L)
  } else if (strand == "-") {
    idx0 <- end - pos
    codon_index <- floor(idx0 / 3L) + 1L
    codon_pos   <- (idx0 %% 3L) + 1L
    codon_end   <- end - 3L * (codon_index - 1L)
    codon_start <- codon_end - 2L
    codon_fw    <- gd_nt_window(entity, seq_id, codon_start, codon_end, strict = TRUE)
    codon_seq   <- gd_revcomp(codon_fw)
    aa_ref      <- gd_translate_literal(codon_seq, transl_table = tt, frame = 1L)
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

gd_summarize_snp_row <- function(gd, i, logger = NULL) {
  stopifnot(inherits(gd, "genome_entity_gd"))
  entity <- gd$entity
  row    <- gd$events[[i]]
  if (!identical(row$type, "SNP")) {
    return(data.frame(
      id = row$id, type = row$type, chrom = row$seq_id, pos = row$position,
      gene = NA_character_, strand = NA_character_,
      dna_before = NA_character_, dna_after = NA_character_,
      aa_before  = NA_character_, aa_after  = NA_character_,
      effect = "undetermined",
      ref = row$snp_ref_base %||% NA_character_, alt = row$snp_alt_base %||% row$col6, length = 1L,
      cds_pos = NA_integer_, aa_pos = NA_integer_,
      messages = I(list("Not a SNP row"))
    ))
  }
  
  seq_id <- as.character(row$seq_id)
  pos    <- as.integer(row$position)
  alt    <- as.character(row$snp_alt_base %||% row$col6)
  ref    <- as.character(row$snp_ref_base %||% NA_character_)
  
  cds_rows <- gd_cds_by_position(entity, seq_id, pos)
  if (nrow(cds_rows) == 0L) {
    return(data.frame(
      id = row$id, type = row$type, chrom = seq_id, pos = pos,
      gene = NA_character_, strand = NA_character_,
      dna_before = NA_character_, dna_after = NA_character_,
      aa_before  = NA_character_, aa_after  = NA_character_,
      effect = "intergenic",
      ref = ref, alt = alt, length = 1L,
      cds_pos = NA_integer_, aa_pos = NA_integer_,
      messages = I(list(character()))
    ))
  }
  
  cds_row <- cds_rows[1, , drop = FALSE]
  ctx     <- gd_codon_context(entity, cds_row, pos)
  
  # Validate reference base if provided
  ref_at_pos <- gd_nt_window(entity, seq_id, pos, pos, strict = TRUE)
  if (!is.na(ref) && nzchar(ref) && toupper(ref) != toupper(ref_at_pos)) {
    stop(sprintf("Ref base '%s' does not match reference '%s' at %s:%d", ref, ref_at_pos, seq_id, pos))
  }
  
  alt_gene     <- if (ctx$strand == "+") toupper(alt) else gd_complement_base(alt)
  codon_alt    <- gd_replace_substring1(ctx$codon_seq_ref, ctx$codon_pos, alt_gene)
  aa_ref_local <- ctx$aa_ref
  aa_alt_local <- gd_translate_literal(codon_alt, transl_table = ctx$transl_table, frame = 1L)
  
  effect <- if (aa_ref_local == aa_alt_local) "synonymous"
  else if (aa_alt_local == "*")     "nonsense"
  else                               "missense"
  
  gene_id <- gd_resolve_gene_id(cds_row)
  
  cds_pos <- {
    start <- gd_parse_int(cds_row$start[1]); end <- gd_parse_int(cds_row$end[1])
    if (ctx$strand == "+") ctx$codon_start - start + ctx$codon_pos
    else                   end - ctx$codon_end + ctx$codon_pos
  }
  aa_pos <- ctx$codon_index
  
  data.frame(
    id = row$id, type = row$type, chrom = seq_id, pos = pos,
    gene = gene_id, strand = ctx$strand,
    dna_before = ctx$codon_seq_ref, dna_after = codon_alt,
    aa_before  = aa_ref_local,      aa_after  = aa_alt_local,
    effect = effect,
    ref = ref_at_pos, alt = alt, length = 1L,
    cds_pos = cds_pos, aa_pos = aa_pos,
    messages = I(list(character()))
  )
}

gd_summarize_effects <- function(gd, logger = NULL) {
  stopifnot(inherits(gd, "genome_entity_gd"))
  out <- lapply(seq_along(gd$events), function(i) {
    row <- gd$events[[i]]
    if (identical(row$type, "SNP")) gd_summarize_snp_row(gd, i, logger = logger) else
      data.frame(
        id = row$id, type = row$type, chrom = row$seq_id, pos = row$position,
        gene = NA_character_, strand = NA_character_,
        dna_before = NA_character_, dna_after = NA_character_,
        aa_before  = NA_character_, aa_after  = NA_character_,
        effect = "undetermined",
        ref = row$snp_ref_base %||% NA_character_,
        alt = (row$snp_alt_base %||% row$col6),
        length = if (!is.null(row$del_size)) row$del_size else NA_integer_,
        cds_pos = NA_integer_, aa_pos = NA_integer_,
        messages = I(list(sprintf("Type '%s' not yet implemented", row$type)))
      )
  })
  do.call(rbind, out)
}


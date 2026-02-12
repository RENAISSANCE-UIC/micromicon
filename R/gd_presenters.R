# Testing RA parsing in preparation for full table
predicted_mutations_ra <- function(gd, min_freq = 0) {
  # ---------- small utils ----------
  
  scalar_or <- function(a, b) {
    if (is.null(a) || length(a) == 0) return(b)
    if (length(a) > 1) a <- a[1]
    if (is.na(a)) return(b)
    a
  }
  
  # Extract first tag value (always returns a length-1 vector)
  tag_first <- function(tags_list, key) {
    vapply(tags_list, function(x) {
      if (is.null(x) || is.null(x[[key]]) || length(x[[key]]) == 0)
        return(NA_character_)
      x[[key]][1]
    }, character(1))
  }
  
  # Extract concatenated multi-value tags (length-1 vector)
  tag_concat <- function(tags_list, key, sep = " | ") {
    vapply(tags_list, function(x) {
      if (is.null(x) || is.null(x[[key]]) || length(x[[key]]) == 0)
        return(NA_character_)
      paste(x[[key]], collapse = sep)
    }, character(1))
  }
  
  # Arrow and slash normalization (breseq uses <br>; we use clean ASCII)
  arrowize <- function(x) {
    ifelse(
      is.na(x),
      x,
      gsub("/", " / ",
           gsub("\n", " | ", x, fixed = TRUE),
           fixed = TRUE
      )
    )
  }
  
  # Numeric + percentage formatting
  fmt_pos  <- function(x) formatC(x, big.mark = ",", format = "d")
  fmt_freq <- function(x) ifelse(is.na(x), NA, sprintf("%.1f%%", 100 * as.numeric(x)))
  
  # Base predicate (SAFE — doesn’t throw on length >1)
  is_base <- function(x) {
    if (is.null(x) || length(x) == 0) return(FALSE)
    x <- x[1]
    !is.na(x) && x %in% c("A","C","G","T")
  }
  
  
  # ---------- flatten evidence (RA) ----------
  ev <- gd_events_table(gd, kinds = "evidence", expand_tags = FALSE)
  ra <- subset(ev, type == "RA")
  n  <- nrow(ra)
  if (n == 0L) {
    return(data.frame(evidence=character(), `seq id`=character(), position=character(),
                      mutation=character(), freq=character(), annotation=character(),
                      gene=character(), description=character(), check.names=FALSE))
  }
  
  # ---------- build a map from RA evidence id -> primary mutation ----------
  # We walk the raw event list to access 'parent_ids' that your parser kept.
  # Pick the first mutation whose parent_ids contains the RA id (most RA have exactly one).
  events <- gd$events
  # Preindex mutations by event index to avoid repeated scans
  mut_idx <- which(vapply(events, function(e) identical(e$kind, "mutation"), logical(1)))
  # For rapid lookup, create a list: for each mutation, a fast set of parent_ids
  mut_parents <- lapply(events[mut_idx], function(m) m$parent_ids %||% integer(0))
  
  ra_id <- ra$id
  # For each RA id, find a mutation index whose parent_ids contains it
  map_mut_i <- vapply(seq_along(ra_id), function(i) {
    rid <- ra_id[i]
    j <- which(vapply(mut_parents, function(p) length(p) && rid %in% p, logical(1)))
    if (length(j)) mut_idx[j[1]] else NA_integer_
  }, integer(1))
  
  # ---------- pull mutation fields we might need for formatting ----------
  # For convenience, extract a small view for found mutations, NA otherwise.
  getm <- function(i, nm, default = NA) {
    if (is.na(map_mut_i[i])) return(default)
    events[[ map_mut_i[i] ]][[nm]] %||% default
  }
  
  # SNP label from RA (if RA truly carries base->base)
  ref <- ra$ra_ref_base
  alt <- ra$ra_new_base
  is_snp_vec <- mapply(function(a, b) is_base(a) && is_base(b), ref, alt)
  snp_label  <- ifelse(is_snp_vec, paste0(ref, "→", alt), NA_character_)
  
  # Mutation-backed labels for non-SNP RA
  # Gather type and per-type fields from the linked mutation (if any)
  mut_type    <- vapply(seq_len(n), function(i) getm(i, "type", NA_character_), character(1))
  del_size    <- vapply(seq_len(n), function(i) getm(i, "del_size", NA_integer_), integer(1))
  ins_seq     <- vapply(seq_len(n), function(i) getm(i, "ins_seq", NA_character_), character(1))
  ins_size    <- vapply(seq_len(n), function(i) getm(i, "ins_size", NA_integer_), integer(1))
  sub_size    <- vapply(seq_len(n), function(i) getm(i, "sub_size", NA_integer_), integer(1))
  sub_new_seq <- vapply(seq_len(n), function(i) getm(i, "sub_new_seq", NA_character_), character(1))
  amp_copies  <- vapply(seq_len(n), function(i) getm(i, "amp_new_copy_number", NA_integer_), integer(1))
  inv_size    <- vapply(seq_len(n), function(i) getm(i, "inv_size", NA_integer_), integer(1))
  
  fmt_bp_del <- function(s) ifelse(is.na(s), "DEL", ifelse(s == 1L, "Δ1 bp", paste0("Δ", format(s, big.mark=","), " bp")))
  fmt_bp_ins <- function(s, seq) {
    if (!is.na(s)) {
      if (s == 1L) "+1 bp" else paste0("+", format(s, big.mark=","), " bp")
    } else if (!is.na(seq) && nchar(seq) == 1L) {
      paste0("+", seq)
    } else if (!is.na(seq) && nchar(seq) > 1L) {
      paste0("+", nchar(seq), " bp")
    } else "INS"
  }
  fmt_sub <- function(sz, newseq) {
    if (is.na(sz)) "SUB"
    else if (sz == 1L && !is.na(newseq) && nzchar(newseq)) paste0("1 bp→", newseq)
    else paste0("SUB(", sz, " bp)")
  }
  fmt_amp <- function(cn) ifelse(is.na(cn), "AMP", paste0("AMP×", cn))
  fmt_inv <- function(sz) ifelse(is.na(sz), "INV", paste0("INV(", format(sz, big.mark=","), " bp)"))
  
  # Format mutation-backed label
  mut_backfill <- vapply(seq_len(n), function(i) {
    t <- mut_type[i]
    if (is.na(t)) return(NA_character_)
    switch(t,
           "DEL" = fmt_bp_del(del_size[i]),
           "INS" = fmt_bp_ins(ins_size[i], ins_seq[i]),
           "SUB" = fmt_sub(sub_size[i], sub_new_seq[i]),
           "AMP" = fmt_amp(amp_copies[i]),
           "INV" = fmt_inv(inv_size[i]),
           # SNP here (rare path for RA non-SNP) — fall back to RA snp_label if available
           "SNP" = snp_label[i] %||% "SNP",
           # other (MOB/CON/etc.) — just print the type
           t
    )
  }, character(1))
  
  # If neither SNP nor mutation backfill is available, try RA 'prediction' tag last.
  ra_pred <- tag_first(ra$tags, "prediction")
  
  mutation <- ifelse(!is.na(snp_label), snp_label,
                     ifelse(!is.na(mut_backfill), mut_backfill,
                            ifelse(!is.na(ra_pred) & nzchar(ra_pred), ra_pred, "RA")
                     ))
  
  # ---------- annotation / gene / description (tags) ----------
  aa_ref_seq  <- tag_first (ra$tags, "aa_ref_seq")
  aa_new_seq  <- tag_first (ra$tags, "aa_new_seq")
  cod_ref_seq <- tag_first (ra$tags, "codon_ref_seq")
  cod_new_seq <- tag_first (ra$tags, "codon_new_seq")
  cod_number  <- tag_first (ra$tags, "codon_number")
  gene_name   <- tag_concat(ra$tags, "gene_name")
  gene_prod   <- tag_concat(ra$tags, "gene_product")
  gene_pos    <- tag_concat(ra$tags, "gene_position")
  snp_type    <- tag_first (ra$tags, "snp_type")
  locus_tag   <- tag_concat(ra$tags, "locus_tag")
  
  have_aa <- !is.na(aa_ref_seq) & !is.na(aa_new_seq) &
    !is.na(cod_ref_seq) & !is.na(cod_new_seq) & !is.na(cod_number)
  
  annotation <- ifelse(
    have_aa,
    paste0(aa_ref_seq, cod_number, aa_new_seq, " (", cod_ref_seq, "→", cod_new_seq, ")"),
    ifelse(!is.na(gene_pos) & nzchar(gene_pos), gene_pos, snp_type)
  )
  
  gene <- ifelse(!is.na(gene_name) & nzchar(gene_name), gene_name, locus_tag)
  gene <- arrowize(gene)
  description <- arrowize(gene_prod)
  
  # ---------- position (include insert position like "860,824:1") ----------
  pos_str <- fmt_pos(ra$position)
  ins_pos <- ra$ra_insert_position
  pos_str <- ifelse(!is.na(ins_pos) & ins_pos > 0L, paste0(pos_str, ":", ins_pos), pos_str)
  
  # ---------- frequency filter ----------
  keep <- rep(TRUE, n)
  if (is.finite(min_freq) && min_freq > 0) {
    keep <- !is.na(ra$ev_frequency) & (ra$ev_frequency > min_freq)
  }
  
  # ---------- build output ----------
  out <- data.frame(
    evidence    = ra$type[keep],
    `seq id`    = ra$seq_id[keep],
    position    = pos_str[keep],
    mutation    = mutation[keep],
    freq        = fmt_freq(ra$ev_frequency[keep]),
    annotation  = annotation[keep],
    gene        = gene[keep],
    description = description[keep],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Stable order: seq id, numeric position
  numpos <- suppressWarnings(as.integer(gsub("[:,].*$", "", gsub(",", "", out$position))))
  o <- order(out$`seq id`, numpos)
  out[o, , drop = FALSE]
}

#' Build a breseq-style "Predicted mutations" table from genome_entity_gd.
#'
#' Columns: evidence, seq id, position, mutation, freq, annotation, gene, description
#'
#' @param gd genome_entity_gd parsed by your v0.2.7 parser
#' @param min_freq numeric in [0,1], applied to RA evidence (e.g., 0.20 for >20%)
#' @param include_structural logical; include MC/JC-backed structural predictions
#' @param join how to present multi-valued tags: "slash" (default), "pipe", or "newline"
#' @return data.frame with breseq-like predicted mutations
predicted_mutations <- function(gd, min_freq = 0, include_structural = TRUE,
                                join = c("slash","pipe","newline")) {
  
  join <- match.arg(join)
  
  # ---------------------------
  # Unborked scalar/tag helpers
  # ---------------------------
  
  scalar_or <- function(a, b) {
     # If it's a list, do NOT truncate—either return the whole list or fallback
     if (is.list(a)) return(if (is.null(a)) b else a)
     if (is.null(a) || length(a) == 0) return(b)
     if (length(a) > 1) a <- a[1]
     if (is.na(a)) return(b)
     a
  }
  tag_first <- function(tags_list, key) {
    vapply(tags_list, function(x) {
      if (is.null(x) || is.null(x[[key]]) || length(x[[key]]) == 0)
        return(NA_character_)
      x[[key]][1]
    }, character(1))
  }
  tag_concat <- function(tags_list, key, sep = " | ") {
    vapply(tags_list, function(x) {
      if (is.null(x) || is.null(x[[key]]) || length(x[[key]]) == 0)
        return(NA_character_)
      paste(x[[key]], collapse = sep)
    }, character(1))
  }
  is_base <- function(x) {
    if (is.null(x) || length(x) == 0) return(FALSE)
    x <- x[1]
    !is.na(x) && x %in% c("A","C","G","T")
  }
  arrowize <- function(x) {
    s <- ifelse(is.na(x), x,
                gsub("/", " / ", gsub("\n", " | ", x, fixed = TRUE), fixed = TRUE))
    if (join == "slash")   gsub("\\|", " / ", s, fixed = TRUE)
    else if (join == "newline") gsub("\\|", "\n", s, fixed = TRUE)
    else s  # "pipe"
  }
  fmt_pos  <- function(x) formatC(x, big.mark = ",", format = "d")
  fmt_freq <- function(x) ifelse(is.na(x), NA, sprintf("%.1f%%", 100 * as.numeric(x)))
  
  # ======================================
  # Phase A — RA (polymorphism) predictions
  # ======================================
  ev <- gd_events_table(gd, kinds = "evidence", expand_tags = FALSE)
  ra <- subset(ev, type == "RA")
  n_ra <- nrow(ra)
  ra_out <- data.frame()
  
  if (n_ra > 0L) {
    # Pull tag fields directly from list-column (length-safe)
    ra_pred      <- tag_first (ra$tags, "prediction")
    aa_ref_seq   <- tag_first (ra$tags, "aa_ref_seq")
    aa_new_seq   <- tag_first (ra$tags, "aa_new_seq")
    cod_ref_seq  <- tag_first (ra$tags, "codon_ref_seq")
    cod_new_seq  <- tag_first (ra$tags, "codon_new_seq")
    cod_number   <- tag_first (ra$tags, "codon_number")
    gene_name    <- tag_concat(ra$tags, "gene_name")
    gene_prod    <- tag_concat(ra$tags, "gene_product")
    gene_pos     <- tag_concat(ra$tags, "gene_position")
    snp_type     <- tag_first (ra$tags, "snp_type")
    locus_tag    <- tag_concat(ra$tags, "locus_tag")
    
    # SNP vs non-SNP RA
    ref <- ra$ra_ref_base; alt <- ra$ra_new_base
    is_snp_vec <- mapply(function(a,b) is_base(a) && is_base(b), ref, alt)
    snp_label  <- ifelse(is_snp_vec, paste0(ref, "→", alt), NA_character_)
    
    # Backfill mutation label for non-SNP RA via linked mutation (parent_ids)
    # Build id -> idx map
    ev_ids <- vapply(gd$events, function(e) scalar_or(e$id, NA_integer_), integer(1))
    idx_by_id <- setNames(seq_along(ev_ids), as.character(ev_ids))
    
    # Map each RA -> first mutation whose parent_ids contains this RA id
    map_mut_i <- vapply(seq_len(n_ra), function(i) {
      rid <- ra$id[i]
      hit <- NA_integer_
      for (j in seq_along(gd$events)) {
        e <- gd$events[[j]]
        if (!identical(e$kind, "mutation")) next
        parents <- scalar_or(e$parent_ids, integer(0))
        if (length(parents) && rid %in% parents) { hit <- j; break }
      }
      hit
    }, integer(1))
    
    getm <- function(i, nm, default = NA) {
      j <- map_mut_i[i]
      if (is.na(j)) return(default)
      scalar_or(gd$events[[j]][[nm]], default)
    }
    
    mut_type    <- vapply(seq_len(n_ra), function(i) getm(i, "type", NA_character_), character(1))
    del_size    <- vapply(seq_len(n_ra), function(i) getm(i, "del_size", NA_integer_), integer(1))
    ins_seq     <- vapply(seq_len(n_ra), function(i) getm(i, "ins_seq", NA_character_), character(1))
    ins_size    <- vapply(seq_len(n_ra), function(i) getm(i, "ins_size", NA_integer_), integer(1))
    sub_size    <- vapply(seq_len(n_ra), function(i) getm(i, "sub_size", NA_integer_), integer(1))
    sub_new_seq <- vapply(seq_len(n_ra), function(i) getm(i, "sub_new_seq", NA_character_), character(1))
    amp_copies  <- vapply(seq_len(n_ra), function(i) getm(i, "amp_new_copy_number", NA_integer_), integer(1))
    inv_size    <- vapply(seq_len(n_ra), function(i) getm(i, "inv_size", NA_integer_), integer(1))
    
    fmt_bp_del <- function(s) ifelse(is.na(s), "DEL", ifelse(s == 1L, "Δ1 bp", paste0("Δ", format(s, big.mark=","), " bp")))
    fmt_bp_ins <- function(s, seq) {
      if (!is.na(s)) {
        if (s == 1L) "+1 bp" else paste0("+", format(s, big.mark=","), " bp")
      } else if (!is.na(seq) && nzchar(seq) && nchar(seq) == 1L) {
        paste0("+", seq)
      } else if (!is.na(seq) && nzchar(seq) && nchar(seq) > 1L) {
        paste0("+", nchar(seq), " bp")
      } else "INS"
    }
    fmt_sub <- function(sz, newseq) {
      if (is.na(sz)) "SUB"
      else if (sz == 1L && !is.na(newseq) && nzchar(newseq)) paste0("1 bp→", newseq)
      else paste0("SUB(", sz, " bp)")
    }
    fmt_amp <- function(cn) ifelse(is.na(cn), "AMP", paste0("AMP×", cn))
    fmt_inv <- function(sz) ifelse(is.na(sz), "INV", paste0("INV(", format(sz, big.mark=","), " bp)"))
    
    mut_backfill <- vapply(seq_len(n_ra), function(i) {
      t <- mut_type[i]
      if (is.na(t)) return(NA_character_)
      switch(t,
             "DEL" = fmt_bp_del(del_size[i]),
             "INS" = fmt_bp_ins(ins_size[i], ins_seq[i]),
             "SUB" = fmt_sub(sub_size[i], sub_new_seq[i]),
             "AMP" = fmt_amp(amp_copies[i]),
             "INV" = fmt_inv(inv_size[i]),
             "SNP" = scalar_or(snp_label[i], "SNP"),
             t
      )
    }, character(1))
    
    mutation <- ifelse(!is.na(snp_label), snp_label,
                       ifelse(!is.na(mut_backfill), mut_backfill,
                              ifelse(!is.na(ra_pred) & nzchar(ra_pred), ra_pred, "RA")
                       ))
    
    have_aa <- !is.na(aa_ref_seq) & !is.na(aa_new_seq) &
      !is.na(cod_ref_seq) & !is.na(cod_new_seq) & !is.na(cod_number)
    annotation <- ifelse(
      have_aa,
      paste0(aa_ref_seq, cod_number, aa_new_seq, " (", cod_ref_seq, "→", cod_new_seq, ")"),
      ifelse(!is.na(gene_pos) & nzchar(gene_pos), gene_pos, snp_type)
    )
    
    gene <- ifelse(!is.na(gene_name) & nzchar(gene_name), gene_name, locus_tag)
    gene <- arrowize(gene)
    description <- arrowize(gene_prod)
    
    pos_str <- fmt_pos(ra$position)
    ins_pos <- ra$ra_insert_position
    pos_str <- ifelse(!is.na(ins_pos) & ins_pos > 0L, paste0(pos_str, ":", ins_pos), pos_str)
    
    keep_ra <- rep(TRUE, n_ra)
    if (is.finite(min_freq) && min_freq > 0) {
      keep_ra <- !is.na(ra$ev_frequency) & (ra$ev_frequency > min_freq)
    }
    
    ra_out <- data.frame(
      evidence    = ra$type[keep_ra],
      `seq id`    = ra$seq_id[keep_ra],
      position    = pos_str[keep_ra],
      mutation    = mutation[keep_ra],
      freq        = fmt_freq(ra$ev_frequency[keep_ra]),
      annotation  = annotation[keep_ra],
      gene        = gene[keep_ra],
      description = description[keep_ra],
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  
  # ======================================
  # Phase B — MC/JC structural predictions
  # (read directly from gd$events; aggregate tags from evidence parents)
  # ======================================
  struct_out <- data.frame()
  if (isTRUE(include_structural)) {
    ev_list <- gd$events
    
    # Build id -> idx map for all events
    all_ids <- vapply(ev_list, function(e) scalar_or(e$id, NA_integer_), integer(1))
    idx_by_id <- setNames(seq_along(all_ids), as.character(all_ids))
    
    # Helper: collect gene/annotation/product from parent evidence (MC/JC)
    collect_structural_tags <- function(parent_ids) {
      ann  <- character(0); gene <- character(0); desc <- character(0)
      strand <- character(0)  # capture gene_strand if present
      
      for (pid in parent_ids) {
        idx <- idx_by_id[[as.character(pid)]]
        if (is.null(idx) || is.na(idx)) next
        evp <- ev_list[[idx]]
        tags <- if (is.null(evp$tags)) list() else evp$tags
        
        gp  <- scalar_or(tags[["gene_position"]], NA_character_)
        gn  <- scalar_or(tags[["gene_name"]],    NA_character_)
        lt  <- scalar_or(tags[["locus_tag"]],    NA_character_)
        pr  <- scalar_or(tags[["gene_product"]], NA_character_)
        gs  <- scalar_or(tags[["gene_strand"]],  NA_character_)
        
        if (!is.na(gp) && nzchar(gp))  ann  <- c(ann,  gp)
        if (!is.na(gn) && nzchar(gn))  gene <- c(gene, gn) else if (!is.na(lt) && nzchar(lt)) gene <- c(gene, lt)
        if (!is.na(pr) && nzchar(pr))  desc <- c(desc, pr)
        if (!is.na(gs) && nzchar(gs))  strand <- c(strand, gs)
      }
      
      # Append arrows by strand if available (coarse, but helpful)
      # If any strand says "-", add " ←"; if "+", add " →".
      finalize_gene <- function(x, strands) {
        if (!length(x)) return(NA_character_)
        g <- paste(unique(x), collapse = " | ")
        if (length(strands)) {
          if (any(strands %in% c("-", "-1", -1))) g <- paste0(g, " \u2190")
          else if (any(strands %in% c("+", "1", 1))) g <- paste0(g, " \u2192")
        }
        g
      }
      
      list(
        annotation  = if (length(ann))  paste(unique(ann),  collapse = " | ") else NA_character_,
        gene        = finalize_gene(gene, strand),
        description = if (length(desc)) paste(unique(desc), collapse = ", ")  else NA_character_
      )
    }
    
    # Select mutation events that have MC and/or JC parents, but no RA
    mut_indices <- which(vapply(ev_list, function(e) identical(e$kind, "mutation"), logical(1)))
    keep_idx <- integer(0)
    ev_str   <- character(0)
    
    for (j in mut_indices) {
      m <- ev_list[[j]]
      parents <- scalar_or(m$parent_ids, integer(0))
      if (!length(parents)) next
      
      p_types <- vapply(parents, function(pid) {
        jj <- idx_by_id[[as.character(pid)]]
        if (is.null(jj) || is.na(jj)) return(NA_character_)
        scalar_or(ev_list[[jj]]$type, NA_character_)
      }, character(1))
      p_types <- unique(na.omit(p_types))
      
      has_RA <- "RA" %in% p_types
      has_MC <- "MC" %in% p_types
      has_JC <- "JC" %in% p_types
      
      if (!has_RA && (has_MC || has_JC)) {
        keep_idx <- c(keep_idx, j)
        ev_str   <- c(ev_str, paste(sort(intersect(p_types, c("MC","JC"))), collapse = " "))
      }
    }
    
    if (length(keep_idx)) {
      n_s <- length(keep_idx)
      seq_id     <- character(n_s)
      position   <- character(n_s)
      mutation   <- character(n_s)
      annotation <- character(n_s)
      gene       <- character(n_s)
      description <- character(n_s)
      
      fmt_bp_del <- function(s) ifelse(is.na(s), "DEL", ifelse(s == 1L, "Δ1 bp", paste0("Δ", format(s, big.mark=","), " bp")))
      fmt_bp_ins <- function(s, seq) {
        if (!is.na(s)) {
          if (s == 1L) "+1 bp" else paste0("+", format(s, big.mark=","), " bp")
        } else if (!is.na(seq) && nzchar(seq) && nchar(seq) == 1L) {
          paste0("+", seq)
        } else if (!is.na(seq) && nzchar(seq) && nchar(seq) > 1L) {
          paste0("+", nchar(seq), " bp")
        } else "INS"
      }
      fmt_sub <- function(sz, newseq) {
        if (is.na(sz)) "SUB"
        else if (sz == 1L && !is.na(newseq) && nzchar(newseq)) paste0("1 bp→", newseq)
        else paste0("SUB(", sz, " bp)")
      }
      fmt_amp <- function(cn) ifelse(is.na(cn), "AMP", paste0("AMP×", cn))
      fmt_inv <- function(sz) ifelse(is.na(sz), "INV", paste0("INV(", format(sz, big.mark=","), " bp)"))
      
      for (ii in seq_along(keep_idx)) {
        j <- keep_idx[ii]
        m <- ev_list[[j]]
        
        # Locus
        seq_id[ii]   <- scalar_or(m$seq_id, NA_character_)
        position[ii] <- fmt_pos(scalar_or(m$position, NA_integer_))
        
        # Label by mutation type
        t <- scalar_or(m$type, NA_character_)
        if (t == "DEL") {
          mutation[ii] <- fmt_bp_del(scalar_or(m$del_size, NA_integer_))
        } else if (t == "INS") {
          mutation[ii] <- fmt_bp_ins(scalar_or(m$ins_size, NA_integer_),
                                     scalar_or(m$ins_seq, NA_character_))
        } else if (t == "SUB") {
          mutation[ii] <- fmt_sub(scalar_or(m$sub_size, NA_integer_),
                                  scalar_or(m$sub_new_seq, NA_character_))
        } else if (t == "AMP") {
          mutation[ii] <- fmt_amp(scalar_or(m$amp_new_copy_number, NA_integer_))
        } else if (t == "INV") {
          mutation[ii] <- fmt_inv(scalar_or(m$inv_size, NA_integer_))
        } else if (t == "SNP") {
          mutation[ii] <- "SNP"
        } else {
          mutation[ii] <- scalar_or(t, "MUT")
        }
        
        # Aggregate annotation/gene/product from parent evidence
        parents <- scalar_or(m$parent_ids, integer(0))
        
        mtags <- if (is.null(m$tags)) list() else m$tags
        
        mut_ann  <- scalar_or(mtags[["gene_position"]], NA_character_)
        mut_gene <- scalar_or(mtags[["gene_name"]],    NA_character_)
        mut_ltag <- scalar_or(mtags[["locus_tag"]],    NA_character_)
        mut_desc <- scalar_or(mtags[["gene_product"]], NA_character_)
        mut_str  <- scalar_or(mtags[["gene_strand"]],  NA_character_)
        
        # If gene_name is missing but locus_tag exists, use locus_tag
        if (is.na(mut_gene) || !nzchar(mut_gene)) mut_gene <- mut_ltag
        
        # Collect MC/JC tags (secondary sources)
        agp <- collect_structural_tags(parents)
        
        # Choose annotation, gene, desc in priority order:
        annotation[ii] <- arrowize(
          if (!is.na(mut_ann)  && nzchar(mut_ann))  mut_ann  else agp$annotation
        )
        
        # For directional arrow: breseq uses ← or → based on gene_strand
        g_final <- if (!is.na(mut_gene) && nzchar(mut_gene)) mut_gene else agp$gene
        if (!is.na(mut_str) && mut_str == "<")  g_final <- paste0(g_final, " \u2190")
        if (!is.na(mut_str) && mut_str == ">")  g_final <- paste0(g_final, " \u2192")
        
        gene[ii] <- arrowize(g_final)
        
        description[ii] <- if (!is.na(mut_desc) && nzchar(mut_desc))
          mut_desc
        else
          agp$description
      }
      
      # breseq displays 100% for these fixed structural predictions (matches your HTML)
      
      freq <- rep("100.0%", n_s)
      
      struct_out <- data.frame(
        evidence    = ev_str,
        `seq id`    = seq_id,
        position    = position,
        mutation    = mutation,
        freq        = freq,
        annotation  = annotation,
        gene        = gene,
        description = description,   # <-- fixed
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      
    }
  }
  
  # ==========================
  # Merge, order, and return
  # ==========================
  out <- rbind(ra_out, struct_out)
  if (nrow(out) == 0L) return(out)
  
  # Stable sort: seq id, then numeric position
  numpos <- suppressWarnings(as.integer(gsub("[:,].*$", "", gsub(",", "", out$position))))
  o <- order(out$`seq id`, numpos, out$mutation)
  out[o, , drop = FALSE]
}








#' Columns: evidence, seq id, position, mutation, freq, annotation, gene, description
#'
#' @param gd genome_entity_gd parsed by parse_gd_annotated()
#' @param min_freq numeric between 0 and 1, applied to RA evidence (e.g., 0.20 for >20%)
#' @param include_structural logical; include MC/JC-backed structural predictions
#' @param join how to present multi-valued tags: "slash" (default), "pipe", or "newline"
#' @return data.frame with breseq-like predicted mutations
#' @keywords internal
predicted_mutations_orig <- function(gd, min_freq = 0, include_structural = TRUE,
                                join = c("slash","pipe","newline")) {
  
  join <- match.arg(join)

  # Unborked scalar/tag helpers
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
  
  # Phase A — RA (polymorphism) predictions
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

  # Phase B — MC/JC structural predictions
  # (read directly from gd$events; aggregate tags from evidence parents)
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
  
  # Merge, order, and return
  out <- rbind(ra_out, struct_out)
  if (nrow(out) == 0L) return(out)
  
  # Stable sort: seq id, then numeric position
  numpos <- suppressWarnings(as.integer(gsub("[:,].*$", "", gsub(",", "", out$position))))
  o <- order(out$`seq id`, numpos, out$mutation)
  out[o, , drop = FALSE]
}




#' Columns: evidence, type, seq_id, position, mutation, freq, annotation, gene, description
#'
#' @param gd genome_entity_gd parsed by parse_gd_annotated()
#' @param min_freq numeric between 0 and 1, applied to RA evidence (e.g., 0.20 for >20%)
#' @param include_structural logical; include MC/JC-backed structural predictions
#' @param join how to present multi-valued tags: "slash" (default), "pipe", or "newline"
#' @return data.frame with breseq-like predicted mutations (now with 3-letter `type`)
#' @keywords internal
predicted_mutations_int <- function(gd, min_freq = 0, include_structural = TRUE,
                                join = c("slash","pipe","newline")) {
  join <- match.arg(join)
  
  # Unborked scalar/tag helpers
  scalar_or <- function(a, b) {
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
    if (join == "slash")        gsub("\\|", " / ", s, fixed = TRUE)
    else if (join == "newline") gsub("\\|", "\n",  s, fixed = TRUE)
    else s  # "pipe"
  }
  fmt_pos  <- function(x) formatC(x, big.mark = ",", format = "d")
  fmt_freq <- function(x) ifelse(is.na(x), NA, sprintf("%.1f%%", 100 * as.numeric(x)))
  
  # Phase A — RA (polymorphism) predictions
  ev <- gd_events_table(gd, kinds = "evidence", expand_tags = FALSE)
  ra <- subset(ev, type == "RA")
  n_ra <- nrow(ra)
  ra_out <- data.frame()
  
  if (n_ra > 0L) {
    # Pull tag fields
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
    
    # Build id -> idx map
    ev_ids   <- vapply(gd$events, function(e) scalar_or(e$id, NA_integer_), integer(1))
    idx_by_id<- setNames(seq_along(ev_ids), as.character(ev_ids))
    
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
             "SNP" = if (!is.na(snp_label[i])) snp_label[i] else "SNP",
             t
      )
    }, character(1))
    
    # Final labels and TYPE for RA rows
    mutation <- ifelse(v_not_na(snp_label), snp_label,
                       ifelse(v_not_na(mut_backfill), mut_backfill,
                              ifelse(v_nzchar(ra_pred), ra_pred, "RA")
                       ))
    
    
    have_aa <- !is.na(aa_ref_seq) & !is.na(aa_new_seq) &
      !is.na(cod_ref_seq) & !is.na(cod_new_seq) & !is.na(cod_number)
    
    annotation <- ifelse(
      have_aa,
      paste0(aa_ref_seq, cod_number, aa_new_seq, " (", cod_ref_seq, "→", cod_new_seq, ")"),
      ifelse(v_nzchar(gene_pos), gene_pos, snp_type)
    )
    
    gene <- ifelse(v_nzchar(gene_name), gene_name, locus_tag)
    gene <- arrowize(gene)
    description <- arrowize(gene_prod)
    
    pos_str <- fmt_pos(ra$position)
    ins_pos <- ra$ra_insert_position
    pos_str <- ifelse(v_pos_gt0(ins_pos), paste0(pos_str, ":", ins_pos), pos_str)
    
    keep_ra <- rep(TRUE, n_ra)
    if (is.finite(min_freq) && min_freq > 0) {
      keep_ra <- v_not_na(ra$ev_frequency) & (ra$ev_frequency > min_freq)
    }
    
    # --- NEW: compute 3-letter TYPE for RA rows (evidence-agnostic) ---
    # Prefer explicit SNP when snp_label is present; else use linked mut_type.
    type_ra <- ifelse(!is.na(snp_label), "SNP", mut_type)
    
    ra_out <- data.frame(
      evidence    = ra$type[keep_ra],        # "RA"
      type        = type_ra[keep_ra],        # <-- added
      seq_id      = ra$seq_id[keep_ra],
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
  
  # Phase B — MC/JC structural predictions (unchanged logic; now emit `type`)
  struct_out <- data.frame()
  if (isTRUE(include_structural)) {
    ev_list <- gd$events
    
    all_ids <- vapply(ev_list, function(e) scalar_or(e$id, NA_integer_), integer(1))
    idx_by_id <- setNames(seq_along(all_ids), as.character(all_ids))
    
    collect_structural_tags <- function(parent_ids) {
      ann <- character(0); gene <- character(0); desc <- character(0); strand <- character(0)
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
      seq_id      <- character(n_s)
      position    <- character(n_s)
      mutation    <- character(n_s)
      annotation  <- character(n_s)
      gene        <- character(n_s)
      description <- character(n_s)
      type_vec    <- character(n_s)  # <-- new
      
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
        
        seq_id[ii]   <- scalar_or(m$seq_id, NA_character_)
        position[ii] <- fmt_pos(scalar_or(m$position, NA_integer_))
        
        t <- scalar_or(m$type, NA_character_)
        type_vec[ii] <- t  # <-- capture the three-letter type
        
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
        
        parents <- scalar_or(m$parent_ids, integer(0))
        mtags <- if (is.null(m$tags)) list() else m$tags
        mut_ann  <- scalar_or(mtags[["gene_position"]], NA_character_)
        mut_gene <- scalar_or(mtags[["gene_name"]],    NA_character_)
        mut_ltag <- scalar_or(mtags[["locus_tag"]],    NA_character_)
        mut_desc <- scalar_or(mtags[["gene_product"]], NA_character_)
        mut_str  <- scalar_or(mtags[["gene_strand"]],  NA_character_)
        if (is.na(mut_gene) || !nzchar(mut_gene)) mut_gene <- mut_ltag
        
        agp <- collect_structural_tags(parents)
        annotation[ii] <- arrowize(if (!is.na(mut_ann) && nzchar(mut_ann))  mut_ann else agp$annotation)
        
        g_final <- if (!is.na(mut_gene) && nzchar(mut_gene)) mut_gene else agp$gene
        if (!is.na(mut_str) && mut_str == "<") g_final <- paste0(g_final, " \u2190")
        if (!is.na(mut_str) && mut_str == ">") g_final <- paste0(g_final, " \u2192")
        gene[ii] <- arrowize(g_final)
        
        description[ii] <- if (!is.na(mut_desc) && nzchar(mut_desc)) mut_desc else agp$description
      }
      
      freq <- rep("100.0%", n_s)
      
      struct_out <- data.frame(
        evidence    = ev_str,
        type        = type_vec,      # <-- added
        seq_id      = seq_id,
        position    = position,
        mutation    = mutation,
        freq        = freq,
        annotation  = annotation,
        gene        = gene,
        description = description,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Merge, order, and return
  out <- rbind(ra_out, struct_out)
  if (nrow(out) == 0L) return(out)
  
  # Stable sort: seq_id, then numeric position
  numpos <- suppressWarnings(as.integer(gsub("[:,].*$", "", gsub(",", "", out$position))))
  o <- order(out$seq_id, numpos, out$mutation)
  out[o, , drop = FALSE]
}





#' Public wrapper: predicted mutations (opinionated defaults; args pass-through)
#'
#' Calls predicted_mutations_int(), announces row count, and prints a preview:
#' - If dplyr is installed: prints tibble (contract unaffected).
#' - If dplyr is not installed: prints base data.frame.
#' Any preview/logging errors are swallowed; the table is still returned.
#'
#' @param gd genome_entity_gd
#' @param ... reserved for future options passed to *_int (kept for forward-compat)
#' @param min_freq numeric scalar, keep rows with freq >= min_freq (default 0)
#' @param include_structural logical, include JC/MC/etc. (default TRUE)
#' @param join one of c("slash","pipe","newline") for multi-item fields (default "slash")
#' @return data.frame/tibble with predicted mutations (superset schema)
#' @export
predicted_mutations <- function(gd, ...,
                                min_freq = 0,
                                include_structural = TRUE,
                                join = c("slash", "pipe", "newline")) {
  gd_assert(gd)
  join <- match.arg(join)
  
  # Build once via internal, regardless of downstream printing path
  tbl <- predicted_mutations_int(
    gd,
    min_freq = min_freq,
    include_structural = include_structural,
    join = join,
    ...
  )
  
  # Announce + print a friendly preview, but never fail the call if printing goes sideways
  tryCatch(
    {
      n <- NROW(tbl)
      # row count announcement
      cli::cli_inform(c(
        "i" = sprintf("predicted_mutations: %d row%s.", n, if (n == 1L) "" else "s")
      ))
      
      # choose printing strategy
      has_dplyr <- requireNamespace("dplyr", quietly = TRUE)
      n_show <- min(n, 25L)
      
      if (has_dplyr) {
        # tibble path (do not modify tbl; only coerce for printing)
        cli::cli_inform(c(
          "i" = sprintf("Showing top %d as a tibble.", n_show)
        ))
        print(utils::head(dplyr::as_tibble(tbl), n_show))
        tbl <- dplyr::as_tibble(tbl)
      } else {
        cli::cli_inform(c(
          "!" = "{.pkg dplyr} not detected; showing a base data.frame preview.",
          "i" = "Install with: {.code install.packages('dplyr')} to prefer tibble display."
        ))
        print(utils::head(tbl, n_show))
      }
    },
    error = function(e) {
      # Stay serene; log and proceed to return the table
      cli::cli_warn(c(
        "!" = "Non-fatal issue while printing the preview.",
        ">" = conditionMessage(e)
      ))
    }
  )
  
  # Always return the fully built table (but quietly)
  invisible(tbl)
}

#' Fallback enrichment for predicted_mutations() output using hoisted features
#'
#' This function does *not* change existing, non-missing fields. It only
#' fills gaps using the reference features in `gd$features`.
#'
#' @param gd          genome_entity_gd
#' @param tbl         data.frame from predicted_mutations(gd)
#' @param want_distance logical; compute intergenic distance when gene is NA
#' @return data.frame with enriched columns (type/gene/annotation/distance if missing)
#' @keywords internal
pm_fallback_enrich <- function(gd, tbl, want_distance = TRUE) {
  stopifnot(inherits(gd, "genome_entity_gd"))
  if (!is.data.frame(tbl) || !nrow(tbl)) return(tbl)
  
  # Defensive column presence
  need <- c("evidence","seq_id","position","mutation","freq","annotation","gene","description")
  miss <- setdiff(need, names(tbl))
  if (length(miss)) {
    stop("pm_fallback_enrich(): tbl is missing required columns: ",
         paste(miss, collapse = ", "))
  }
  
  # 0) Extract numeric position (left of ":" if present; strip commas)
  pos_num <- suppressWarnings(
    as.integer(gsub("[:,].*$", "", gsub(",", "", tbl$position)))
  )
  seq_id  <- as.character(tbl$seq_id)
  
  # 1) Fallback 'type' if absent (or NA) by reading the human 'mutation' label
  if (!"type" %in% names(tbl)) {
    tbl$type <- NA_character_
  }
  fallback_type <- function(mut_string) {
    if (is.na(mut_string) || !nzchar(mut_string)) return(NA_character_)
    s <- trimws(mut_string)
    if (grepl("^Δ\\s*\\d+\\s*bp\\s*$", s)) return("DEL")
    if (grepl("^\\+\\s*\\d+\\s*bp\\s*$", s)) return("INS")
    if (grepl("→", s, fixed = TRUE)) {
      parts <- strsplit(s, "→", fixed = TRUE)[[1]]
      lhs <- gsub("[^ACGTN]", "", parts[1], perl = TRUE)
      rhs <- gsub("[^ACGTN\\*]", "", parts[2], perl = TRUE)
      if (nchar(lhs) == 1L && nchar(rhs) == 1L) return("SNP")
      return("SUB")
    }
    # Last resort: leave NA (e.g., AMP/INV without explicit tokens)
    NA_character_
  }
  # Only fill where type is missing/NA
  for (i in seq_len(nrow(tbl))) {
    if (is.na(tbl$type[i]) || !nzchar(tbl$type[i])) {
      tbl$type[i] <- fallback_type(tbl$mutation[i])
    }
  }
  
  # 2) Feature lookup (overlaps and nearest)
  # We will use gd$features (GFF3-like) where type=="CDS" or "gene"
  feats <- gd$features
  if (!is.data.frame(feats) || !nrow(feats)) return(tbl)  # nothing to do
  
  # Normalize feature columns we will use
  f_seq   <- as.character(feats$seqname)
  f_start <- as.integer(feats$start)
  f_end   <- as.integer(feats$end)
  f_strand<- as.character(feats$strand)
  f_type  <- as.character(feats$type)
  f_name  <- as.character(feats$Name)
  f_id    <- as.character(feats$ID)
  
  # We'll prefer CDS names; fall back to gene name or ID
  is_gene_like <- f_type %in% c("gene","CDS")
  # Pre-filter to speed up
  idx_keep <- which(is_gene_like & gd_not_na(f_start) & gd_not_na(f_end))
  if (!length(idx_keep)) return(tbl)
  
  f_seq   <- f_seq[idx_keep]
  f_start <- f_start[idx_keep]
  f_end   <- f_end[idx_keep]
  f_strand<- f_strand[idx_keep]
  f_type  <- f_type[idx_keep]
  f_name  <- f_name[idx_keep]
  f_id    <- f_id[idx_keep]
  
  # Helper: choose display symbol
  feat_symbol <- function(i) {
    nm <- f_name[i]
    if (!is.na(nm) && nzchar(nm)) return(nm)
    id <- f_id[i]
    if (!is.na(id) && nzchar(id)) return(id)
    return(NA_character_)
  }
  
  # For each row with missing gene and/or annotation, fill from overlap;
  # if intergenic and want_distance=TRUE, fill distance as "intergenic (+d/-d)"
  if (!"distance" %in% names(tbl)) {
    tbl$distance <- NA_character_
  }
  
  for (i in seq_len(nrow(tbl))) {
    # Skip rows that already have a gene + annotation
    need_gene <- !(gd_nzchar(tbl$gene[i]))
    need_ann  <- !(gd_nzchar(tbl$annotation[i]))
    
    if (!need_gene && !need_ann && (!want_distance || gd_nzchar(tbl$distance[i]))) next
    
    sid <- seq_id[i]
    pos <- pos_num[i]
    if (is.na(sid) || is.na(pos)) next
    
    # Candidate features on same contig
    hit_idx <- which(f_seq == sid & f_start <= pos & f_end >= pos)
    
    if (length(hit_idx)) {
      # Prefer CDS if available at this locus
      cds_idx <- hit_idx[f_type[hit_idx] == "CDS"]
      pick <- if (length(cds_idx)) cds_idx[1L] else hit_idx[1L]
      
      if (need_gene) {
        sym <- feat_symbol(pick)
        # Append arrow if strand available
        if (!is.na(f_strand[pick]) && f_strand[pick] %in% c("+","-")) {
          arrow <- if (f_strand[pick] == "+") " \u2192" else " \u2190"
          sym <- paste0(sym, arrow)
        }
        tbl$gene[i] <- sym
      }
      
      if (need_ann) {
        # Create a coarse position string "coding (x/y nt)" is not trivial here,
        # so fall back to "within <symbol>" to be explicit.
        sym <- feat_symbol(pick)
        tbl$annotation[i] <- if (!is.na(sym) && nzchar(sym)) {
          paste0("within ", sym)
        } else {
          "within feature"
        }
      }
      
      # distance for overlapping case is 0 (optional)
      if (want_distance) {
        tbl$distance[i] <- "0"
      }
      
    } else if (isTRUE(want_distance)) {
      # Intergenic: find nearest upstream and downstream gene-like features
      same_seq <- which(f_seq == sid)
      if (!length(same_seq)) next
      
      # distances to feature intervals
      # If pos < start: distance = start - pos (upstream feature)
      # If pos > end  : distance = pos - end  (downstream feature)
      d <- rep(NA_integer_, length(same_seq))
      for (k in seq_along(same_seq)) {
        j <- same_seq[k]
        if (pos < f_start[j]) d[k] <- f_start[j] - pos
        else if (pos > f_end[j]) d[k] <- pos - f_end[j]
        else d[k] <- 0L
      }
      
      # nearest feature by distance
      finite_idx <- which(!is.na(d))
      if (length(finite_idx)) {
        kmin <- finite_idx[ which.min(d[finite_idx]) ]
        j    <- same_seq[kmin]
        dist <- d[kmin]
        
        # Sign: negative if upstream (pos < start), positive if downstream (pos > end)
        sign_ch <- if (pos < f_start[j]) "-" else if (pos > f_end[j]) "+" else "0"
        if (need_ann) {
          sym <- feat_symbol(j)
          tbl$annotation[i] <- if (!is.na(sym) && nzchar(sym)) {
            paste0("intergenic (", sign_ch, dist, ") near ", sym)
          } else {
            paste0("intergenic (", sign_ch, dist, ")")
          }
        }
        if (gd_not_na(dist)) {
          tbl$distance[i] <- paste0(sign_ch, dist)
        }
        
        # If gene is completely missing, provide the nearest symbol (without implying overlap)
        if (need_gene) {
          sym <- feat_symbol(j)
          if (!is.na(sym) && nzchar(sym)) {
            tbl$gene[i] <- sym
          }
        }
      }
    }
  }
  
  # 3) Final hygiene: ensure columns in preferred order
  # evidence, type, seq_id, position, mutation, freq, annotation, gene, description, (optional) distance
  keep <- c("evidence","type","seq_id","position","mutation","freq","annotation","gene","description",
            intersect("distance", names(tbl)))
  extra <- setdiff(names(tbl), keep)
  tbl <- tbl[c(keep, extra)]
  
  tbl
}


#' Return the standardized predicted mutations tibble
#' Contract: (type, evidence, seq_id, position, end, ref, alt, freq,
#'            gene, locus_tag, strand, annotation, product, distance, tags, row_id)
#' @keywords internal
pm_tbl <- function(gd) {
  gd_assert(gd)
  
  # Require dplyr for tibble return type
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    cli::cli_abort(
      c(
        "x" = "The {.pkg dplyr} package is required to return a tibble.",
        "i" = "Install with: {.code install.packages('dplyr')}"
      )
    )
  }
  
  tb <- predicted_mutations_int(gd)
  
  # Enforce your pared-down column contract
  required_cols <- c("evidence","type","seq_id","position","mutation","freq","annotation","gene")
  missing <- setdiff(required_cols, names(tb))
  if (length(missing)) {
    cli::cli_abort(
      c(
        "x" = "pm_tbl: required column(s) missing.",
        "!" = paste("Missing:", paste(missing, collapse = ", ")),
        "i" = "Confirm that {.fn predicted_mutations} produces these fields upstream."
      )
    )
  }
  
  # Return a tibble (explicitly, via dplyr)
  dplyr::as_tibble(tb[, required_cols, drop = FALSE])
}


#' Quick human scan (top-n). Does not alter data, just prints.
#' Use cli for a little context; return the tibble (invisibly) for piping.
#' @keywords internal
pm_view <- function(gd, n = 25) {
  gd_assert(gd)
  n <- as.integer(`%||%`(n, 25L))
  
  has_dplyr <- requireNamespace("dplyr", quietly = TRUE)
  
  if (has_dplyr) {
    # Contracted tibble path
    tb <- pm_tbl(gd)  # will cli_abort if something else is wrong
    cli::cli_inform(c(
      "i" = sprintf("predicted_mutations (contracted tibble): %d rows; showing top %d.", nrow(tb), n)
    ))
    print(utils::head(tb, n))
    return(invisible(tb))
  }
  
  # Fallback: no dplyr installed — show something useful, don’t crash.
  cli::cli_inform(c(
    "!" = "{.pkg dplyr} not detected; showing a base preview from {.fn predicted_mutations}.",
    "i" = "Install with: {.code install.packages('dplyr')} to enable tibble output and focused columns."
  ))
  df <- predicted_mutations(gd)
  # Be defensive: if the object is huge, head() keeps it polite.
  print(utils::head(df, n))
  invisible(df)
}

#' Focus by exact gene symbol (case-sensitive, exact match).
#' Returns a tibble (0+ rows). Logs counts.
#' @keywords internal
pm_focus_gene <- function(gd, gene) {
  gd_assert(gd)
  tb <- pm_tbl(gd)
  out <- tb[!is.na(tb$gene) & tb$gene == gene, , drop = FALSE]
  cli::cli_inform(c("i" = paste0("pm_focus_gene: '", gene, "' → ", nrow(out), " row(s).")))
  out
}

#' Focus by genomic range (inclusive) on a specific seq_id, using the *contracted* table.
#' Behavior: position-only filtering (no 'end' dependence), exact seq_id match.
#' @keywords internal
pm_focus_roi <- function(gd, seq_id, start, end) {
  gd_assert(gd)
  stopifnot(!missing(seq_id), !missing(start), !missing(end))
  
  # Contracted 8-col view from your wrapper path
  tb <- pm_tbl(gd)
  
  # Normalize and validate inputs
  start <- suppressWarnings(as.integer(start))
  end   <- suppressWarnings(as.integer(end))
  if (is.na(start) || is.na(end)) {
    cli::cli_abort(c(
      "x" = "pm_focus_roi: start/end must be coercible to integer.",
      "i" = sprintf("Got start=%s, end=%s.", as.character(substitute(start)), as.character(substitute(end)))
    ))
  }
  if (start > end) { tmp <- start; start <- end; end <- tmp }
  
  # Coerce table fields defensively (position may be character/factor)
  pos <- suppressWarnings(as.integer(as.character(tb$position)))
  sid <- as.character(tb$seq_id)
  target_sid <- as.character(seq_id)
  
  # Build mask; treat NA position as non-hit
  hit <- (!is.na(pos)) & (sid == target_sid) & (pos >= start) & (pos <= end)
  if (length(hit) == 0L) hit <- rep(FALSE, nrow(tb))  # 0-row safety
  
  out <- tb[hit, , drop = FALSE]
  
  cli::cli_inform(c(
    "i" = paste0("pm_focus_roi: seq_id=", target_sid, " [", start, ", ", end, "] → ",
                 nrow(out), " row(s).")
  ))
  
  # Preserve return type preference
  if (requireNamespace("dplyr", quietly = TRUE)) {
    return(dplyr::as_tibble(out))
  }
  out
}



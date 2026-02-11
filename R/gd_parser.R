#' Parse a breseq annotated.gd file with type-aware dispatch
#'
#' Parses a breseq annotated.gd file (produced by gdtools ANNOTATE) and creates
#' a genome_entity_gd object. The parser uses type-aware dispatch to handle
#' different mutation types (SNP, DEL, INS, SUB, MOB, AMP, CON, INV), evidence
#' types (RA, MC, JC, UN), and validation types (TSEQ, PFLP, etc.). Events are
#' automatically binned by kind (mutation/evidence/validation) for efficient filtering.
#'
#' @param gd_path Path to the annotated.gd file (not raw output.gd)
#' @param entity A genome_entity object (reference genome; fields will be hoisted)
#' @param strict Logical; if TRUE, stop on inconsistencies; if FALSE, warn and continue
#' @param fasta_path Optional path to FASTA file for provenance checksums
#' @param gff3_path Optional path to GFF3 file for provenance checksums
#' @param gbk_path Optional path to GenBank file for provenance checksums
#' @return A genome_entity_gd object with parsed events binned by kind in
#'   \code{provenance$by_kind} (mutation_idx, evidence_idx, validation_idx)
#' @export
#' @examples
#' \dontrun{
#'   ref <- read_genome("reference.gbk")
#'   gd <- parse_gd_annotated("annotated.gd", ref, strict = TRUE)
#'   length(gd$events)
#'   gd$provenance$by_kind$mutation_idx  # indices of mutation events
#' }
parse_gd_annotated <- function(gd_path, entity, strict = TRUE,
                               fasta_path = NULL, gff3_path = NULL, 
                               gbk_path = NULL) {
  if (!file.exists(gd_path)) stop("File does not exist: ", gd_path)
  stopifnot(inherits(entity, "genome_entity"))
  
  # ---- local scalar helpers 
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  first_non_na <- function(x) {
    if (is.null(x) || length(x) == 0L) return(NA)
    i <- which(!is.na(x))
    if (length(i)) x[i[1]] else x[1]
  }
  scalar_coalesce <- function(a, b) {
    if (is.null(a))                  return(b)
    if (length(a) == 0L)             return(b)
    if (length(a) == 1L && is.na(a)) return(b)
    a
  }
  
  # ---- local helpers used by the parser -
  .normalize_seq_id <- function(s) {
    if (is.null(s)) return(NA_character_)
    s <- as.character(s); s <- trimws(s)
    if (identical(s, ".") || !nzchar(s)) NA_character_ else s
  }
  .parse_pair <- function(x, sep = "/", as = c("int","num","char")) {
    as <- match.arg(as)
    if (is.null(x) || !length(x)) return(c(NA, NA))
    y <- strsplit(x[1], sep, fixed = TRUE)[[1]]
    if (length(y) < 2) y <- c(y, NA_character_)
    if (as == "int") suppressWarnings(as.integer(y))
    else if (as == "num") suppressWarnings(as.numeric(y))
    else y
  }
  .as_int <- function(x) suppressWarnings(as.integer(if (length(x)) x else NA))
  .as_num <- function(x) suppressWarnings(as.numeric(if (length(x)) x else NA))
  .parse_tags <- function(tags_raw) {
    out <- list()
    if (!length(tags_raw)) return(out)
    for (t in tags_raw) {
      kv <- strsplit(t, "=", fixed = TRUE)[[1]]
      k  <- kv[1]
      v  <- if (length(kv) > 1) kv[2] else NA_character_
      if (k %in% names(out)) out[[k]] <- c(out[[k]], v) else out[[k]] <- v
    }
    out
  }
  # Skip bounds if seq_id is NA; otherwise validate contig and position
  .bounds_check <- function(i, contig, pos, strict, contig_lengths) {
    if (is.null(contig) || is.na(contig)) return(invisible(TRUE))  # skip
    if (!contig %in% names(contig_lengths)) {
      msg <- sprintf("GD contig '%s' not found in locked reference; no liftover attempted.", contig)
      if (strict) stop(msg) else return(invisible(NULL))
    }
    ln <- contig_lengths[[contig]]
    if (!is.null(pos) && !is.na(pos) && (pos < 1L || pos > ln)) {
      msg <- sprintf("Event %d position %s out of bounds for contig %s [1..%d].", i, pos, contig, ln)
      if (strict) stop(msg)
    }
    invisible(TRUE)
  }
  
  # ---------- read and pre-filter 
  lines <- readLines(gd_path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  
  if (!is_annotated_gd(lines)) {
    stop(
      "This appears to be a raw breseq output.gd (annotation tags not detected).\n",
      "Supply an annotated.gd produced by gdtools ANNOTATE together with the reference."
    )
  }
  
  header_idx   <- grep("^#", lines)
  body_idx     <- setdiff(seq_along(lines), header_idx)
  header_lines <- lines[header_idx]
  body_lines   <- lines[body_idx]
  
  # Header key/values
  header <- list()
  for (h in header_lines) {
    h2 <- sub("^#+=*", "", h)
    h2 <- trimws(h2)
    sp <- strsplit(h2, "[[:space:]]+", perl = TRUE)[[1]]
    key <- sp[1]; val <- paste(sp[-1], collapse = " ")
    if (!nzchar(key)) next
    if (key %in% names(header)) header[[key]] <- c(header[[key]], val) else header[[key]] <- val
  }
  
  # Provenance + reference
  ref_manifest   <- reference_manifest_from_genome_entity(entity, fasta_path, gff3_path, gbk_path)
  contig_lengths <- setNames(as.numeric(entity$metadata$length_bp), entity$metadata$seqname)
  
  MUT_TYPES <- c("SNP","SUB","DEL","INS","MOB","AMP","CON","INV")
  EVI_TYPES <- c("RA","MC","JC","UN")
  VAL_TYPES <- c("TSEQ","PFLP","RFLP","PFGE","PHYL","CURA","NOTE","FPOS")
  
  # ---------- type-specific parsers 
  .parse_mut <- function(row, type) {
    # mutation: type, id, parent-ids, seq_id, position, ...
    idx <- 4L
    f   <- function(k) if (length(row) >= k) row[k] else NA_character_
    
    parent_ids_raw <- if (length(row) >= 3L) row[3L] else NA_character_
    parent_ids     <- if (is.na(parent_ids_raw) || parent_ids_raw %in% c(".", "")) integer(0) else {
      suppressWarnings(as.integer(strsplit(parent_ids_raw, ",", fixed = TRUE)[[1]]))
    }
    
    if (type == "SNP") {
      seq_id   <- .normalize_seq_id(f(idx + 0L))
      position <- .as_int(f(idx + 1L))
      new_seq  <- f(idx + 2L)
      return(list(
        seq_id = seq_id, position = position, contig = seq_id,
        snp_ref_base = NA_character_, snp_alt_base = new_seq,
        col6 = new_seq, col7 = NA, col8 = NA,
        parent_ids_raw = parent_ids_raw, parent_ids = parent_ids,
        .fixed_end = idx + 2L
      ))
    }
    
    if (type == "SUB") {
      seq_id   <- .normalize_seq_id(f(idx + 0L))
      position <- .as_int(f(idx + 1L))
      size     <- .as_int(f(idx + 2L))
      new_seq  <- f(idx + 3L)
      return(list(
        seq_id = seq_id, position = position, contig = seq_id,
        sub_size = size, sub_new_seq = new_seq,
        col6 = f(idx + 2L), col7 = f(idx + 3L), col8 = NA,
        parent_ids_raw = parent_ids_raw, parent_ids = parent_ids,
        .fixed_end = idx + 3L
      ))
    }
    
    if (type == "DEL") {
      seq_id   <- .normalize_seq_id(f(idx + 0L))
      position <- .as_int(f(idx + 1L))
      size     <- .as_int(f(idx + 2L))
      return(list(
        seq_id = seq_id, position = position, contig = seq_id,
        del_size = size,
        col6 = f(idx + 2L), col7 = NA, col8 = NA,
        parent_ids_raw = parent_ids_raw, parent_ids = parent_ids,
        .fixed_end = idx + 2L
      ))
    }
    
    if (type == "INS") {
      seq_id   <- .normalize_seq_id(f(idx + 0L))
      position <- .as_int(f(idx + 1L))
      new_seq  <- f(idx + 2L)
      return(list(
        seq_id = seq_id, position = position, contig = seq_id,
        ins_seq = new_seq,
        ins_size = if (!is.na(new_seq)) nchar(new_seq) else NA_integer_,
        col6 = new_seq, col7 = NA, col8 = NA,
        parent_ids_raw = parent_ids_raw, parent_ids = parent_ids,
        .fixed_end = idx + 2L
      ))
    }
    
    if (type == "MOB") {
      seq_id   <- .normalize_seq_id(f(idx + 0L))
      position <- .as_int(f(idx + 1L))
      rname    <- f(idx + 2L)
      strand   <- .as_int(f(idx + 3L))
      dup      <- .as_int(f(idx + 4L))
      return(list(
        seq_id = seq_id, position = position, contig = seq_id,
        mob_repeat_name = rname, mob_strand = strand, mob_duplication_size = dup,
        col6 = rname, col7 = f(idx + 3L), col8 = f(idx + 4L),
        parent_ids_raw = parent_ids_raw, parent_ids = parent_ids,
        .fixed_end = idx + 4L
      ))
    }
    
    if (type == "AMP") {
      seq_id     <- .normalize_seq_id(f(idx + 0L))
      position   <- .as_int(f(idx + 1L))
      size       <- .as_int(f(idx + 2L))
      new_copies <- .as_int(f(idx + 3L))
      return(list(
        seq_id = seq_id, position = position, contig = seq_id,
        amp_size = size, amp_new_copy_number = new_copies,
        col6 = f(idx + 2L), col7 = f(idx + 3L), col8 = NA,
        parent_ids_raw = parent_ids_raw, parent_ids = parent_ids,
        .fixed_end = idx + 3L
      ))
    }
    
    if (type == "CON") {
      seq_id   <- .normalize_seq_id(f(idx + 0L))
      position <- .as_int(f(idx + 1L))
      size     <- .as_int(f(idx + 2L))
      region   <- f(idx + 3L)
      return(list(
        seq_id = seq_id, position = position, contig = seq_id,
        con_size = size, con_region = region,
        col6 = f(idx + 2L), col7 = f(idx + 3L), col8 = NA,
        parent_ids_raw = parent_ids_raw, parent_ids = parent_ids,
        .fixed_end = idx + 3L
      ))
    }
    
    if (type == "INV") {
      seq_id   <- .normalize_seq_id(f(idx + 0L))
      position <- .as_int(f(idx + 1L))
      size     <- .as_int(f(idx + 2L))
      return(list(
        seq_id = seq_id, position = position, contig = seq_id,
        inv_size = size,
        col6 = f(idx + 2L), col7 = NA, col8 = NA,
        parent_ids_raw = parent_ids_raw, parent_ids = parent_ids,
        .fixed_end = idx + 2L
      ))
    }
    
    stop("Unhandled mutation type: ", type)
  }
  
  .parse_evidence <- function(row, type) {
    # evidence: type, evidence-id, then seq_id/... (type-specific)
    idx <- 3L
    f   <- function(k) if (length(row) >= k) row[k] else NA_character_
    
    
    if (type == "RA") {
      # RA: seq_id, position, insert_position, ref_base, new_base
      toks <- if (length(row) >= idx) row[idx:length(row)] else character(0)
      
      .is_base_char <- function(x) {
        if (is.null(x) || length(x) == 0L) return(FALSE)
        x <- as.character(x[1])
        # Accept DNA bases and '.' as used by breseq for placeholders
        grepl("^[ACGTN\\.]$", x, ignore.case = FALSE)
      }
      .is_nonneg <- function(x) !is.na(x) && x >= 0L
      .has_sid   <- function(s) { s <- .normalize_seq_id(s); !is.na(s) && nzchar(s) }
      
      attempt <- function(o) {
        need <- o + 5L
        if (length(toks) < need) return(NULL)
        sid <- .normalize_seq_id(toks[o + 1L])
        pos <- .as_int(toks[o + 2L])
        ins <- .as_int(toks[o + 3L])
        ref <- toks[o + 4L]
        new <- toks[o + 5L]
        
        ok  <- .is_nonneg(ins) && .is_base_char(ref) && .is_base_char(new)
        # score: prefer real seq_id, then prefer ref not "0" (guards against shifted digits)
        score <- (if (.has_sid(sid)) 2L else 0L) + (ref != "0")
        list(ok = ok, score = score,
             sid = sid, pos = pos, ins = ins, ref = ref, new = new, used = need)
      }
      
      cand <- list(attempt(0L), attempt(1L), attempt(2L))
      cand <- Filter(Negate(is.null), cand)
      ok   <- Filter(function(z) isTRUE(z$ok), cand)
      
      pick <- NULL
      if (length(ok)) {
        # choose the highest score; ties go to the earliest offset
        sc   <- vapply(ok, function(z) z$score, integer(1))
        pick <- ok[[ which.max(sc) ]]
      } else if (length(cand)) {
        # fall back to offset 0 so tags still get parsed
        pick <- cand[[1L]]
      } else {
        return(list(.fixed_end = idx - 1L))  # no fixed fields; let the caller parse tags
      }
      
      return(list(
        seq_id = pick$sid, position = pick$pos, contig = pick$sid,
        ra_insert_position = pick$ins,
        ra_ref_base        = pick$ref,
        ra_new_base        = pick$new,
        .fixed_end         = (idx - 1L) + pick$used
      ))
    }
    
    if (type == "MC") {
      seq_id      <- .normalize_seq_id(f(idx + 0L))
      start       <- .as_int(f(idx + 1L))
      end         <- .as_int(f(idx + 2L))
      start_range <- .as_int(f(idx + 3L))
      end_range   <- .as_int(f(idx + 4L))
      return(list(
        seq_id = seq_id, contig = seq_id,
        start = start, end = end, mc_start_range = start_range, mc_end_range = end_range,
        .fixed_end = idx + 4L
      ))
    }
    
    if (type == "JC") {
      # JC has 7 fixed scalars; try offsets 0..2 to avoid stray placeholder columns.
      toks <- if (length(row) >= idx) row[idx:length(row)] else character(0)
      .is_strand <- function(x) !is.na(x) && x %in% c(-1L, 1L)
      .valid_pos <- function(x) !is.na(x) && x >= 0L
      
      attempt <- function(o) {
        need <- o + 7L
        if (length(toks) < need) return(NULL)
        s1_id <- .normalize_seq_id(toks[o + 1L])
        s1_po <- .as_int(toks[o + 2L])
        s1_st <- .as_int(toks[o + 3L])
        s2_id <- .normalize_seq_id(toks[o + 4L])
        s2_po <- .as_int(toks[o + 5L])
        s2_st <- .as_int(toks[o + 6L])
        ovlp  <- .as_int(toks[o + 7L])
        
        ok <- .is_strand(s1_st) && .is_strand(s2_st) &&
          .valid_pos(s1_po) && .valid_pos(s2_po) &&
          .valid_pos(ovlp)
        
        list(ok = ok,
             s1_id = s1_id, s1_po = s1_po, s1_st = s1_st,
             s2_id = s2_id, s2_po = s2_po, s2_st = s2_st,
             ovlp  = ovlp,
             used  = need)
      }
      
      cand <- list(attempt(0L), attempt(1L), attempt(2L))
      pick <- NULL
      for (c in cand) { if (!is.null(c) && isTRUE(c$ok)) { pick <- c; break } }
      if (is.null(pick)) {
        pick <- attempt(0L)
        if (is.null(pick)) return(list(.fixed_end = idx - 1L))
      }
      
      return(list(
        side_1_seq_id = pick$s1_id,
        side_1_position = pick$s1_po,
        side_1_strand = pick$s1_st,
        side_2_seq_id = pick$s2_id,
        side_2_position = pick$s2_po,
        side_2_strand = pick$s2_st,
        jc_overlap = pick$ovlp,
        .fixed_end = (idx - 1L) + pick$used
      ))
    }
    
    if (type == "UN") {
      seq_id <- .normalize_seq_id(f(idx + 0L))
      start  <- .as_int(f(idx + 1L))
      end    <- .as_int(f(idx + 2L))
      return(list(
        seq_id = seq_id, contig = seq_id, start = start, end = end,
        .fixed_end = idx + 2L
      ))
    }
    
    stop("Unhandled evidence type: ", type)
  }
  
  .rectify_for_constructor <- function(ev) {
    # Universal columns the constructor/validator require
    ev$type      <- ev$type      %||% NA_character_
    ev$id        <- ev$id        %||% NA_integer_
    ev$seq_id    <- ev$seq_id    %||% NA_character_
    ev$position  <- ev$position  %||% NA_integer_
    ev$contig    <- ev$contig    %||% ev$seq_id
    ev$tags      <- ev$tags      %||% list()
    ev$raw_line  <- ev$raw_line  %||% NA_character_
    ev$hash      <- ev$hash      %||% NA_character_
    
    # Legacy col6/7/8
    ev$col6 <- if (!is.null(ev$col6)) ev$col6 else NA
    ev$col7 <- if (!is.null(ev$col7)) ev$col7 else NA
    ev$col8 <- if (!is.null(ev$col8)) ev$col8 else NA
    
    # Legacy rank (not in GD spec for mutations)
    ev$rank <- ev$rank %||% NA_integer_
    
    # Evidence flags
    if (!is.null(ev$kind) && ev$kind == "evidence") {
      ev$is_evidence    <- TRUE
      ev$evidence_class <- ev$type
    } else {
      ev$is_evidence    <- FALSE
      ev$evidence_class <- NA_character_
    }
    
    # Descriptors (ensure presence)
    defaults_chr <- c("snp_alt_base","snp_ref_base","ins_seq","sub_new_seq","mob_repeat_name","con_region")
    defaults_int <- c("del_size","ins_size","sub_size","mob_strand","mob_duplication_size",
                      "amp_size","amp_new_copy_number","con_size","inv_size")
    for (nm in defaults_chr) if (is.null(ev[[nm]])) ev[[nm]] <- NA_character_
    for (nm in defaults_int) if (is.null(ev[[nm]])) ev[[nm]] <- NA_integer_
    
    # Evidence numerics (ensure presence)
    evid_num <- c("ev_frequency","ev_quality","ev_ref_cov_1","ev_ref_cov_2",
                  "ev_new_cov_1","ev_new_cov_2","ev_tot_cov_1","ev_tot_cov_2",
                  "ev_alignment_overlap","ev_cov_minus","ev_cov_plus",
                  "ev_pos_start","ev_pos_end")
    for (nm in evid_num) if (is.null(ev[[nm]])) ev[[nm]] <- NA_real_
    
    # pos_start/pos_end from tags if present
    if (is.na(ev$ev_pos_start) && length(ev$tags)) {
      ev$ev_pos_start <- suppressWarnings(as.integer(ev$tags[["position_start"]][1]))
    }
    if (is.na(ev$ev_pos_end) && length(ev$tags)) {
      ev$ev_pos_end <- suppressWarnings(as.integer(ev$tags[["position_end"]][1]))
    }
    
    # Evidence rows that lack a single canonical (seq_id, position)
    if (!is.null(ev$kind) && ev$kind == "evidence") {
      if (ev$type == "MC" || ev$type == "UN") {
        ev$position <- scalar_coalesce(ev$position, ev$start)
        ev$contig   <- scalar_coalesce(ev$contig,   ev$seq_id)
      } else if (ev$type == "JC") {
        # Choose a representative locus deterministically
        sid <- first_non_na(c(ev$seq_id, ev$side_1_seq_id, ev$side_2_seq_id))
        pos <- first_non_na(c(ev$position, ev$side_1_position, ev$side_2_position, ev$start, ev$end))
        ev$seq_id   <- scalar_coalesce(ev$seq_id,   sid)
        ev$position <- scalar_coalesce(ev$position, pos)
        ev$contig   <- scalar_coalesce(ev$contig,   ev$seq_id)
      }
    }
    
    ev
  }
  
  # ---------- parse body main loop
  events   <- vector("list", length(body_lines))
  kind_vec <- character(length(body_lines))
  
  for (i in seq_along(body_lines)) {
    raw <- body_lines[i]
    row <- strsplit(raw, "\t", fixed = TRUE)[[1]]
    if (length(row) < 2L) {
      msg <- sprintf("Malformed GD line %d: too few fields.", i)
      if (strict) stop(msg) else next
    }
    
    type <- row[1]
    id   <- .as_int(row[2])
    
    if (type %in% MUT_TYPES) {
      ev   <- .parse_mut(row, type)
      kind <- "mutation"
      if (!is.null(ev$seq_id) && !is.na(ev$seq_id) && !is.null(ev$position) && !is.na(ev$position)) {
        .bounds_check(i, ev$seq_id, ev$position, strict, contig_lengths)
      }
      
    } else if (type %in% EVI_TYPES) {
      ev   <- .parse_evidence(row, type)
      kind <- "evidence"
      
      # evidence-specific bounds
      if (type %in% c("RA","UN")) {
        sid <- ev$seq_id
        p   <- first_non_na(c(ev$position, ev$start))
        if (!is.null(sid) && !is.na(sid) && !is.null(p) && !is.na(p)) {
          .bounds_check(i, sid, p, strict, contig_lengths)
        }
      } else if (type == "MC") {
        sid <- ev$seq_id
        if (!is.null(sid) && !is.na(sid) && !is.na(ev$start)) .bounds_check(i, sid, ev$start, strict, contig_lengths)
        if (!is.null(sid) && !is.na(sid) && !is.na(ev$end))   .bounds_check(i, sid, ev$end,   strict, contig_lengths)
      } else if (type == "JC") {
        if (!is.na(ev$side_1_seq_id) && !is.na(ev$side_1_position))
          .bounds_check(i, ev$side_1_seq_id, ev$side_1_position, strict, contig_lengths)
        if (!is.na(ev$side_2_seq_id) && !is.na(ev$side_2_position))
          .bounds_check(i, ev$side_2_seq_id, ev$side_2_position, strict, contig_lengths)
      }
      
    } else if (type %in% VAL_TYPES) {
      ev <- list(seq_id = .normalize_seq_id(if (length(row) >= 3L) row[3L] else NA_character_),
                 .fixed_end = length(row))
      kind <- "validation"
      
    } else {
      if (strict) stop(sprintf("Unknown GD type '%s' at line %d.", type, i)) else next
    }
    
    # Slice trailing tags using sentinel (falls back to row length)
    fixed_end <- scalar_coalesce(ev$.fixed_end, length(row))
    tags_raw  <- if (length(row) > fixed_end) row[(fixed_end + 1L):length(row)] else character(0)
    taglist   <- .parse_tags(tags_raw)
    
    # Envelope
    ev$type     <- type
    ev$id       <- id
    ev$tags     <- taglist
    ev$raw_line <- raw
    ev$kind     <- kind
    
    # Evidence numerics (if applicable)
    if (kind == "evidence") {
      ev$ev_frequency <- .as_num(taglist[["frequency"]])
      ev$ev_quality   <- .as_num(taglist[["quality"]])
      ref_pair <- .parse_pair(taglist[["ref_cov"]], sep = "/", as = "int")
      new_pair <- .parse_pair(taglist[["new_cov"]], sep = "/", as = "int")
      tot_pair <- .parse_pair(taglist[["tot_cov"]], sep = "/", as = "int")
      ev$ev_ref_cov_1 <- ref_pair[1]; ev$ev_ref_cov_2 <- ref_pair[2]
      ev$ev_new_cov_1 <- new_pair[1]; ev$ev_new_cov_2 <- new_pair[2]
      ev$ev_tot_cov_1 <- tot_pair[1]; ev$ev_tot_cov_2 <- tot_pair[2]
      ev$ev_alignment_overlap <- .as_int(taglist[["alignment_overlap"]])
      ev$ev_cov_minus         <- .as_int(taglist[["coverage_minus"]])
      ev$ev_cov_plus          <- .as_int(taglist[["coverage_plus"]])
    }
    
    # Hash from raw fixed slice + tags (pre-rectification)
    row_fixed <- row[seq_len(fixed_end)]
    ev$hash <- canonical_event_hash(
      fixed_fields = row_fixed,
      taglist      = taglist
    )
    
    # Drop sentinel and rectify for constructor's rectangle
    ev$.fixed_end <- NULL
    ev <- .rectify_for_constructor(ev)
    
    # Store
    events[[i]] <- ev
    kind_vec[i] <- kind
  }
  
  # Package object
  gd_checksum <- gd_digest(gd_path, file = TRUE)
  
  obj <- new_genome_entity_gd(
    header    = header,
    events    = events,
    file      = list(path = normalizePath(gd_path), checksum = gd_checksum),
    entity    = entity,
    reference = ref_manifest,
    strict    = strict
  )
  
  # Optional: keep kind indices for convenience
  obj$provenance$by_kind <- list(
    mutation_idx   = which(kind_vec == "mutation"),
    evidence_idx   = which(kind_vec == "evidence"),
    validation_idx = which(kind_vec == "validation")
  )
  
  validate_genome_entity_gd(obj, strict = strict)
}

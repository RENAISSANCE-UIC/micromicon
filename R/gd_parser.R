#' Parse a breseq annotated.gd bound to a Mode A genome_entity
#'
#' @param gd_path Path to annotated.gd (annotated only; raw output.gd is rejected)
#' @param entity A locked genome_entity created by micromicon::read_genome()
#' @param strict If TRUE, halt on any mismatch; if FALSE, warn loudly and continue
#' @param fasta_path Optional path to the FASTA used for Mode A, for file checksum
#' @param gff3_path  Optional path to the GFF3 used for Mode A, for file checksum
#' @param gbk_path   Optional path to the GBK used for Mode A, for file checksum
#' @return An object of class genome_entity_gd
#' @export
parse_gd_annotated <- function(gd_path, entity, strict = TRUE,
                               fasta_path = NULL, gff3_path = NULL, 
                               gbk_path = NULL) {
  if (!file.exists(gd_path)) stop("File does not exist: ", gd_path)
  stopifnot(inherits(entity, "genome_entity"))
  
  lines <- readLines(gd_path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  
  # Fail fast: annotated only (presence of annotation tags)
  if (!is_annotated_gd(lines)) {
    stop(
      "This appears to be a raw breseq output.gd (annotation tags not detected).\n",
      "Please run breseq and supply the annotated.gd plus either the original input GBK ",
      "or the output FASTA+GFF3."
    )
  }
  
  header_idx   <- grep("^#", lines)
  body_idx     <- setdiff(seq_along(lines), header_idx)
  header_lines <- lines[header_idx]
  body_lines   <- lines[body_idx]
  
  # Header map (keep all #=KEY value pairs)
  header <- list()
  for (h in header_lines) {
    h2 <- sub("^#+=*", "", h)
    h2 <- trimws(h2)
    sp <- strsplit(h2, "[[:space:]]+", perl = TRUE)[[1]]
    key <- sp[1]; val <- paste(sp[-1], collapse = " ")
    if (!nzchar(key)) next
    if (key %in% names(header)) header[[key]] <- c(header[[key]], val) else header[[key]] <- val
  }
  
  # Provenance from Mode A
  ref_manifest   <- reference_manifest_from_genome_entity(entity, fasta_path, gff3_path, gbk_path)
  contig_lengths <- setNames(as.numeric(entity$metadata$length_bp), entity$metadata$seqname)
  
  events <- vector("list", length(body_lines))
  for (i in seq_along(body_lines)) {
    raw <- body_lines[i]
    row <- strsplit(raw, "\t", fixed = TRUE)[[1]]
    n   <- length(row)
    
    # At minimum we need type,id/rank,seq_id,position (5 cols)
    if (n < 5L) {
      msg <- sprintf("Malformed GD line %d: expected ≥5 tab-separated fields, got %d.", i, n)
      if (strict) stop(msg) else cli::cli_warn(msg)
    }
    
    # Rectangular slice: always first 8 fixed columns (pad with NA), rest -> tags
    fixed <- rep(NA_character_, 8L)
    upto  <- min(8L, n)
    fixed[1:upto] <- row[1:upto]
    tags_raw  <- if (n > 8L) row[(8L + 1L):n] else character(0)
    
    type     <- fixed[1]
    id       <- suppressWarnings(as.integer(fixed[2]))
    rank     <- suppressWarnings(as.integer(fixed[3]))
    seq_id   <- fixed[4]                                   # contig / reference fragment id
    position <- suppressWarnings(as.integer(fixed[5]))
    
    # Parse tag name=value into a list (preserve multiplicity)
    taglist <- list()
    if (length(tags_raw)) {
      for (t in tags_raw) {
        kv  <- strsplit(t, "=", fixed = TRUE)[[1]]
        k   <- kv[1]
        val <- if (length(kv) > 1) kv[2] else NA_character_
        if (k %in% names(taglist)) taglist[[k]] <- c(taglist[[k]], val) else taglist[[k]] <- val
      }
    }
    
    # Enforce contig identity and bounds vs Mode A (no liftover/remap)
    contig <- seq_id
    if (!nzchar(contig)) {
      msg <- sprintf("Event %d lacks contig (seq_id in column 4). Refusing to infer.", i)
      if (strict) stop(msg) else cli::cli_warn(msg)
    } else if (!contig %in% names(contig_lengths)) {
      msg <- sprintf("GD contig '%s' not found in locked reference; no liftover will be attempted.", contig)
      if (strict) stop(msg) else cli::cli_warn(msg)
    }
    if (!is.na(position) && nzchar(contig) && contig %in% names(contig_lengths)) {
      ln <- contig_lengths[[contig]]
      if (position < 1L || position > ln) {
        msg <- sprintf("Event %d position %d out of bounds for contig %s [1..%d].", i, position, contig, ln)
        if (strict) stop(msg) else cli::cli_warn(msg)
      }
    }
    
    # Raw payloads (rectangular provenance)
    col6 <- fixed[6]
    col7 <- fixed[7]
    col8 <- fixed[8]
    
    # ------------------ Descriptive mutation fields (NA by default) ------------------
    snp_alt_base          <- NA_character_
    snp_ref_base          <- NA_character_
    
    del_size              <- NA_integer_
    
    ins_seq               <- NA_character_
    ins_size              <- NA_integer_
    
    sub_size              <- NA_integer_
    sub_new_seq           <- NA_character_
    
    mob_repeat_name       <- NA_character_
    mob_strand            <- NA_integer_     # 1 or -1
    mob_duplication_size  <- NA_integer_
    
    # Populate per type (no mutation of provenance)
    if (identical(type, "SNP")) {
      snp_alt_base <- col6                                  # GD column 6 is new_seq (ALT)
      snp_ref_base <- (taglist[["ref_seq"]] %||% NA_character_)[1]  # REF from tag
    } else if (identical(type, "DEL")) {
      suppressWarnings({
        maybe_len <- as.integer(col6)                       # GD column 6 is size
        if (!is.na(maybe_len)) del_size <- maybe_len
      })
    } else if (identical(type, "INS")) {
      ins_seq  <- col6                                      # GD column 6 is inserted sequence
      ins_size <- if (!is.na(ins_seq)) nchar(ins_seq) else NA_integer_
    } else if (identical(type, "SUB")) {
      suppressWarnings({
        maybe_len <- as.integer(col6)                       # GD column 6 is size
        if (!is.na(maybe_len)) sub_size <- maybe_len
      })
      sub_new_seq <- col7                                   # GD column 7 is new sequence
    } else if (identical(type, "MOB")) {
      mob_repeat_name     <- col6                           # repeat_name
      suppressWarnings(mob_strand <- as.integer(col7))      # 1 / -1
      suppressWarnings(mob_duplication_size <- as.integer(col8))
    }
    
    # ------------------ Evidence fields (NA by default) ------------------------------
    is_evidence    <- nchar(type) == 2L                     # RA/JC/MC/UN are 2 letters
    evidence_class <- if (is_evidence) type else NA_character_
    
    ev_frequency   <- .as_num(taglist[["frequency"]])       # common in RA/JC
    ev_quality     <- .as_num(taglist[["quality"]])         # RA
    
    ref_cov_pair   <- .parse_pair(taglist[["ref_cov"]], sep = "/", as = "int")
    new_cov_pair   <- .parse_pair(taglist[["new_cov"]], sep = "/", as = "int")
    tot_cov_pair   <- .parse_pair(taglist[["tot_cov"]], sep = "/", as = "int")
    
    ev_ref_cov_1   <- ref_cov_pair[1]
    ev_ref_cov_2   <- ref_cov_pair[2]
    ev_new_cov_1   <- new_cov_pair[1]
    ev_new_cov_2   <- new_cov_pair[2]
    ev_tot_cov_1   <- tot_cov_pair[1]
    ev_tot_cov_2   <- tot_cov_pair[2]
    
    ev_alignment_overlap <- .as_int(taglist[["alignment_overlap"]])  # JC
    ev_cov_minus         <- .as_int(taglist[["coverage_minus"]])     # JC
    ev_cov_plus          <- .as_int(taglist[["coverage_plus"]])      # JC
    
    # Some annotated.gd include explicit position_start/position_end tags (DEL/INS etc.)
    ev_pos_start <- .as_int(taglist[["position_start"]])
    ev_pos_end   <- .as_int(taglist[["position_end"]])
    
    # ------------------ Assemble event ----------------------------------------------
    ev <- list(
      # rectangular core
      type     = type,
      id       = id,
      rank     = rank,
      seq_id   = seq_id,
      position = position,
      col6     = col6,
      col7     = col7,
      col8     = col8,
      contig   = contig,
      
      # full fidelity
      tags     = taglist,
      raw_line = raw,
      
      # descriptive mutation fields (NA when not applicable)
      snp_alt_base         = snp_alt_base,
      snp_ref_base         = snp_ref_base,
      del_size             = del_size,
      ins_seq              = ins_seq,
      ins_size             = ins_size,
      sub_size             = sub_size,
      sub_new_seq          = sub_new_seq,
      mob_repeat_name      = mob_repeat_name,
      mob_strand           = mob_strand,
      mob_duplication_size = mob_duplication_size,
      
      # evidence flags + descriptors (NA when not applicable)
      is_evidence          = is_evidence,
      evidence_class       = evidence_class,
      ev_frequency         = ev_frequency,
      ev_quality           = ev_quality,
      ev_ref_cov_1         = ev_ref_cov_1,
      ev_ref_cov_2         = ev_ref_cov_2,
      ev_new_cov_1         = ev_new_cov_1,
      ev_new_cov_2         = ev_new_cov_2,
      ev_tot_cov_1         = ev_tot_cov_1,
      ev_tot_cov_2         = ev_tot_cov_2,
      ev_alignment_overlap = ev_alignment_overlap,
      ev_cov_minus         = ev_cov_minus,
      ev_cov_plus          = ev_cov_plus,
      ev_pos_start         = ev_pos_start,
      ev_pos_end           = ev_pos_end
    )
    
    # Identity hash built from rectangular fixed + tags (ignore descriptive fields)
    ev$hash <- canonical_event_hash(
      fixed_fields = c(ev$type, ev$id, ev$rank, ev$seq_id, ev$position, ev$col6, ev$col7, ev$col8),
      taglist = ev$tags
    )
    
    events[[i]] <- ev
  }
  
  gd_checksum <- gd_digest(gd_path, file = TRUE)
  
  obj <- new_genome_entity_gd(
    header     = header,
    events     = events,
    file       = list(path = normalizePath(gd_path), checksum = gd_checksum),
    entity     = entity,
    reference  = ref_manifest,
    strict     = strict
  )
  
  validate_genome_entity_gd(obj, strict = strict)
}




# GD HELPERS

# Note: %||% operator is defined in R/operators.R

#' Compute MD5 checksum for file or object
#' @keywords internal
gd_fingerprint <- function(x, file = FALSE, algo = "md5") {
  if (!identical(tolower(algo), "md5")) {
    cli::cli_warn(c(
      "Only 'md5' is supported without external dependencies",
      "i" = "Ignoring algo = '{algo}'"
    ))
  }
  
  if (isTRUE(file)) {
    # tools::md5sum returns a named vector; unname to keep old shape
    return(unname(tools::md5sum(x)))
  }
  
  # In-memory object: serialize to a stable RDS on disk, md5 that file
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  # version = 2 for stability across R >= 3.x; matches common practice for reproducibility
  saveRDS(x, tf, version = 2, compress = FALSE)
  unname(tools::md5sum(tf))
}

#' @keywords internal
gd_digest <- function(x, file = FALSE, algo = "xxhash64") {
  gd_fingerprint(x, file = file, algo = "md5")
}

#' Check if GD file is annotated
#' @keywords internal
is_annotated_gd <- function(lines) {
  if (!length(lines)) return(FALSE)
  header <- grep("^#", lines, value = TRUE)
  body   <- grep("^[^#]", lines, value = TRUE)
  if (!length(header) || !length(body)) return(FALSE)
  anno_keys <- c("gene","locus_tag","product","protein_id","aa_position",
                 "codon_pos","gene_position","gene_name")
  rx <- paste0("(^|\\t)(", paste(anno_keys, collapse = "|"), ")=")
  any(grepl(rx, body, perl = TRUE))
}

#' Compute canonical hash for GD event
#' @keywords internal
canonical_event_hash <- function(fixed_fields, taglist) {
  fixed_str <- paste(fixed_fields, collapse = "\t")
  if (length(taglist)) {
    keys <- sort(names(taglist))
    tag_pairs <- unlist(lapply(keys, function(k) {
      vals <- taglist[[k]]
      if (length(vals) == 0L) return(sprintf("%s=", k))
      paste(sprintf("%s=%s", k, vals), collapse = ";")
    }), use.names = FALSE)
    basis <- paste(fixed_str, paste(tag_pairs, collapse = "\t"), sep = "\t")
  } else {
    basis <- fixed_str
  }
  gd_digest(basis)
}

#' Create reference manifest from genome entity
#' @keywords internal
reference_manifest_from_genome_entity <- function(entity, fasta_path = NULL,
                                                  gff3_path = NULL, gbk_path = NULL) {
  stopifnot(inherits(entity, "genome_entity"))
  contigs <- data.frame(
    name   = entity$metadata$seqname,
    length = as.numeric(entity$metadata$length_bp),
    stringsAsFactors = FALSE
  )
  checksums <- list()
  if (!is.null(fasta_path) && file.exists(fasta_path)) checksums$fasta <- gd_digest(fasta_path, file = TRUE)
  if (!is.null(gff3_path)  && file.exists(gff3_path))  checksums$gff3  <- gd_digest(gff3_path,  file = TRUE)
  if (!is.null(gbk_path)   && file.exists(gbk_path))   checksums$gbk   <- gd_digest(gbk_path,   file = TRUE)
  
  if (length(checksums) == 0L) {
    # Synthesize deterministic digests when file paths aren't available
    seqs <- entity$sequences$dna_raw
    nm <- names(seqs); ord <- order(nm)
    seq_payload <- paste(paste0(">", nm[ord]), seqs[ord], collapse = "\n")
    checksums$fasta_synth <- gd_digest(seq_payload)
    if (is.data.frame(entity$features)) {
      cols <- intersect(c("seqname","start","end","strand","type","ID","Name","Alias","locus_tag","product"),
                        names(entity$features))
      feat_payload <- paste(utils::capture.output(print(entity$features[, cols, drop = FALSE])), collapse = "\n")
      checksums$gff3_synth <- gd_digest(feat_payload)
    }
  }
  
  list(contigs = contigs, checksums = checksums)
}

#' Parse pair of values
#' @keywords internal
.parse_pair <- function(x, sep = "/", as = c("int","num")) {
  if (is.null(x) || !length(x) || is.na(x[1])) return(c(NA, NA))
  as <- match.arg(as)
  sp <- strsplit(x[1], sep, fixed = TRUE)[[1]]
  if (length(sp) != 2) return(c(NA, NA))
  if (as == "int") {
    suppressWarnings(c(as.integer(sp[1]), as.integer(sp[2])))
  } else {
    suppressWarnings(c(as.numeric(sp[1]), as.numeric(sp[2])))
  }
}

#' Coerce to numeric
#' @keywords internal
.as_num <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  suppressWarnings(as.numeric(x[1]))
}

#' Coerce to integer
#' @keywords internal
.as_int <- function(x) {
  if (is.null(x) || !length(x)) return(NA_integer_)
  suppressWarnings(as.integer(x[1]))
}


#' Check if column has non-NA values
#' @keywords internal
gd_col_has <- function(row, nm) !is.null(row[[nm]]) && !all(is.na(row[[nm]]))

#' Get first value from row column
#' @keywords internal
gd_get1    <- function(row, nm) as.character(row[[nm]][1])

#' Check if row has qualifiers
#' @keywords internal
gd_qual_has <- function(row) gd_col_has(row, "qualifiers") && length(row$qualifiers[[1]]) > 0

#' Get qualifier value from row
#' @keywords internal
gd_get_qual <- function(row, key, default = NA_character_) {
  if (!gd_qual_has(row)) return(default)
  q <- row$qualifiers[[1]]
  if (is.null(names(q)) || !(key %in% names(q))) return(default)
  val <- q[[key]]
  if (length(val) == 0) return(default)
  as.character(val[[1]])
}

#' Check if position is within CDS row
#' @keywords internal
gd_pos_in_row <- function(cds_row, pos) {
  # Fast path: single contiguous CDS
  if (!gd_col_has(cds_row, "location_type") || identical(as.character(cds_row$location_type[1]), "single")) {
    s <- gd_parse_int(cds_row$start[1]); e <- gd_parse_int(cds_row$end[1])
    return(!is.na(s) && !is.na(e) && pos >= s && pos <= e)
  }
  
  # Multi-segment via 'ranges' list-col (PGAP)
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
  
  # Fallback: parse 'location_string' like "complement(123..456,789..999)"
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

#' Filter CDS features covering position
#' @keywords internal
gd_filter_cds_covering <- function(feat, seq_id, pos) {
  cand <- feat[feat$seqname == seq_id & feat$type == "CDS", , drop = FALSE]
  keep <- vapply(seq_len(nrow(cand)), function(i) gd_pos_in_row(cand[i, , drop = FALSE], pos), logical(1))
  cand[keep, , drop = FALSE]
}

#' Check if value is not NA
#' @keywords internal
gd_not_na <- function(x) {
  !is.na(x)
}

#' Check if character value is non-NA and non-empty
#' @keywords internal
gd_nzchar <- function(x) {
  !is.na(x) & nzchar(x)
}

#' Parse integer from string, stripping non-numeric characters
#' @keywords internal
gd_parse_int <- function(x) {
  if (is.numeric(x)) return(as.integer(x))
  x <- as.character(x)
  x <- gsub("[^0-9-]", "", x, perl = TRUE)
  suppressWarnings(as.integer(x))
}



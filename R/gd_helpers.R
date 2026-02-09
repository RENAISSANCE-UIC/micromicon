# GD HELPERS

`%||%` <- function(a, b) if (!is.null(a)) a else b

gd_digest <- function(x, file = FALSE, algo = "xxhash64") {
  if (file) digest::digest(file = x, algo = algo) else digest::digest(x, algo = algo)
}

# Conservative heuristic: annotated.gd should have gene-level tags somewhere
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

.as_num <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  suppressWarnings(as.numeric(x[1]))
}

.as_int <- function(x) {
  if (is.null(x) || !length(x)) return(NA_integer_)
  suppressWarnings(as.integer(x[1]))
}

# Build a reference manifest from a genome_entity (Mode A)
# If you can pass the file paths you used to create Mode A, we store true file checksums.
# Otherwise we synthesize checksums deterministically from object content.
# DO WE NEED THIS??? ----
reference_manifest_from_genome_entity <- function(entity, 
                                                  fasta_path = NULL, 
                                                  gff3_path = NULL, 
                                                  gbk_path = NULL) {
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
  
  # If no file paths available, synthesize
  if (length(checksums) == 0L) {
    # fasta-like digest from dna_raw
    seqs <- entity$sequences$dna_raw
    # ensure deterministic order
    nm <- names(seqs)
    ord <- order(nm)
    seq_payload <- paste(paste0(">", nm[ord]), seqs[ord], collapse = "\n")
    checksums$fasta_synth <- gd_digest(seq_payload)
    # gff3-like digest from features
    if (is.data.frame(entity$features)) {
      # stable key subset
      cols <- intersect(c("seqname","start","end","strand","type","ID","Name","Alias","locus_tag","product"), names(entity$features))
      feat_payload <- paste(utils::capture.output(print(entity$features[ , cols, drop = FALSE])), collapse = "\n")
      checksums$gff3_synth <- gd_digest(feat_payload)
    }
  }
  
  list(contigs = contigs, checksums = checksums)
}

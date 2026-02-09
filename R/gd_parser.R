
#' Compute a fast, stable checksum for a string or file
gd_digest <- function(x, file = FALSE, algo = "xxhash64") {
  if (file) {
    return(digest::digest(file = x, algo = algo))
  } else {
    return(digest::digest(x, algo = algo))
  }
}

# Minimal, opinionated check that a GD is "annotated" (not raw)
# Heuristic: at least one body line must contain typical annotation tags,
# and header must exist. Adjust the tag set to your taste.
.is_annotated_gd <- function(lines) {
  if (!length(lines)) return(FALSE)
  header <- grep("^#", lines, value = TRUE)
  body   <- grep("^[^#]", lines, value = TRUE)
  if (!length(header) || !length(body)) return(FALSE)
  
  # Search for any tag that only appears when reference annotation is available
  anno_keys <- c("gene", "locus_tag", "product", "protein_id", "aa_position",
                 "codon_pos", "gene_position", "gene_product", "gene_name")
  # Does any body line contain key= pattern for those keys?
  rx <- paste0("(^|\\t)(", paste(anno_keys, collapse = "|"), ")=")
  any(grepl(rx, body, perl = TRUE))
}

# Canonicalize event fields and tags into a stable fingerprint
.canonical_fingerprint <- function(fixed_fields, taglist) {
  # fixed: a character vector of the first up-to-7 columns in order
  fixed_str <- paste(fixed_fields, collapse = "\t")
  
  if (length(taglist)) {
    # Sort keys for determinism, and within each key keep the *value order*
    # to preserve multiplicity semantics across repeats of the same tag.
    keys <- sort(names(taglist))
    tag_pairs <- unlist(lapply(keys, function(k) {
      vals <- taglist[[k]]
      if (length(vals) == 0L) return(sprintf("%s=", k))
      paste(sprintf("%s=%s", k, vals), collapse = ";")
    }), use.names = FALSE)
    tags_str <- paste(tag_pairs, collapse = "\t")
    basis <- paste(fixed_str, tags_str, sep = "\t")
  } else {
    basis <- fixed_str
  }
  gd_digest(basis)
}

#' Parse a breseq annotated.gd into a genome_entity_gd
#'
#' @param path Path to annotated.gd (annotated ONLY; raw output.gd is rejected)
#' @param reference Optional Mode A reference metadata, a list with:
#'   - contigs: data.frame(name, length)
#'   - checksums: list(fasta = <md5/xxh>, gff3 = <md5/xxh> | gbk = <md5/xxh>)
#' @param strict If TRUE, fail on any mismatch; if FALSE, warn and continue.
#' @return An object of class "genome_entity_gd"
#' @export
parse_gd_annotated <- function(path, reference = NULL, strict = TRUE) {
  if (!file.exists(path)) {
    stop("File does not exist: ", path)
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]  # drop empties; retain order
  
  # --- INPUT DISCIPLINE: annotated only
  if (!.is_annotated_gd(lines)) {
    stop(paste0(
      "This appears to be a raw breseq output.gd (no annotation tags detected).\n",
      "Run breseq with annotation and supply `annotated.gd` plus either the ",
      "original input GBK or the output FASTA+GFF3.\n",
      "Refusal conditions: no gene/locus/product-style tags found in body."
    ))
  }
  
  header_idx <- grep("^#", lines)
  body_idx   <- setdiff(seq_along(lines), header_idx)
  header_lines <- lines[header_idx]
  body_lines   <- lines[body_idx]
  
  # --- HEADER: strip '#=' prefix (support '#', '#=' and '#==')
  header <- list()
  for (h in header_lines) {
    h2  <- sub("^#+=*", "", h)
    h2  <- trimws(h2)
    # Split on first run of whitespace to get key and remainder
    sp  <- strsplit(h2, "[[:space:]]+", perl = TRUE)[[1]]
    key <- sp[1]
    val <- paste(sp[-1], collapse = " ")
    if (is.null(key) || !nzchar(key)) next
    
    if (key %in% names(header)) {
      header[[key]] <- c(header[[key]], val)
    } else {
      header[[key]] <- val
    }
  }
  
  # --- BODY → events
  events <- vector("list", length(body_lines))
  
  for (i in seq_along(body_lines)) {
    raw <- body_lines[i]
    row <- strsplit(raw, "\t", fixed = TRUE)[[1]]
    n   <- length(row)
    
    if (n < 6L) {
      msg <- sprintf("Malformed GD line %d: expected ≥6 tab-separated fields, got %d.", i, n)
      if (strict) stop(msg) else warning(msg)
    }
    
    fixed <- row[1:min(7, n)]
    tags  <- if (n > 7) row[8:n] else character(0)
    
    # Parse tags into list; preserve multiplicity
    taglist <- list()
    if (length(tags)) {
      for (t in tags) {
        kv  <- strsplit(t, "=", fixed = TRUE)[[1]]
        k   <- kv[1]
        val <- if (length(kv) > 1) kv[2] else NA_character_
        if (k %in% names(taglist)) {
          taglist[[k]] <- c(taglist[[k]], val)
        } else {
          taglist[[k]] <- val
        }
      }
    }
    
    # Assemble event; coerce fixed numeric fields when present
    ev <- list(
      type          = fixed[1],
      id            = suppressWarnings(as.integer(fixed[2])),
      rank          = suppressWarnings(as.integer(fixed[3])),
      genome        = suppressWarnings(as.integer(fixed[4])),
      position      = suppressWarnings(as.integer(fixed[5])),
      ref_seq_fixed = fixed[6],
      alt_or_len    = if (length(fixed) >= 7) fixed[7] else NA_character_,
      tags          = taglist,
      raw_line      = raw
    )
    ev$hash <- .canonical_fingerprint(
      fixed_fields = c(ev$type, ev$id, ev$rank, ev$genome, ev$position, ev$ref_seq_fixed, ev$alt_or_len),
      taglist = ev$tags
    )
    
    events[[i]] <- ev
  }
  
  # --- Build object
  gd_checksum <- gd_digest(path, file = TRUE)
  obj <- new_genome_entity_gd(
    header   = header,
    events   = events,
    file     = list(path = normalizePath(path), checksum = gd_checksum),
    reference = reference,
    strict   = strict
  )
  
  # Validate against reference if provided
  obj <- validate_genome_entity_gd(obj, strict = strict)
  
  obj
}



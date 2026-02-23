
#' Assert length-one logical and return a scalar
#' @keywords internal
gd_as_bool1 <- function(x, what = "value") {
  if (length(x) != 1L || (!isTRUE(is.na(x)) && !is.logical(x))) {
    cli::cli_abort("{.arg {what}} must be a length-one logical")
  }
  as.logical(x)
}


#' Bounds check for 1-indexed string positions
#' @keywords internal
gd_check_pos1 <- function(pos, n, what = "position") {
  if (!is.numeric(pos) || length(pos) != 1L || is.na(pos) || pos < 1L || pos > n) {
    cli::cli_abort("{.arg {what}} out of range [1, {n}] (got {pos})")
  }
  invisible(pos)
}

# ---- DNA utilities ----

#' Reverse-complement a DNA string
#' Falls back to a project-level reverse_complement() if present.
#' @keywords internal
gd_revcomp <- function(dna) {
  if (exists("reverse_complement", mode = "function")) {
    return(reverse_complement(dna))
  }
  if (!is.character(dna) || length(dna) != 1L) {
    cli::cli_abort("{.arg dna} must be a length-one character string")
  }
  s <- strsplit(dna, "", fixed = TRUE)[[1]]
  s <- chartr("ACGTNacgtn", "TGCANtgcan", s)
  paste0(rev(s), collapse = "")
}

#' Complement a base (vectorized; accepts A/C/G/T/N)
#' @keywords internal
gd_complement_base <- function(b) {
  if (!is.character(b)) cli::cli_abort("{.arg b} must be character")
  out <- chartr("ACGTNacgtn", "TGCANtgcan", b)
  bad <- is.na(out) | !nzchar(out)
  if (any(bad)) {
    cli::cli_abort("Invalid base(s) for complement: {paste(unique(b[bad]), collapse = ', ')}")
  }
  toupper(out)
}

#' Replace a single character at a 1-based position
#' @keywords internal
gd_replace_substring1 <- function(s, pos, newchar) {
  if (!is.character(s) || length(s) != 1L) cli::cli_abort("{.arg s} must be length-one character")
  n <- nchar(s, type = "chars", allowNA = FALSE)
  gd_check_pos1(pos, n, "pos")
  pre  <- if (pos > 1L) substr(s, 1L, pos - 1L) else ""
  post <- if (pos < n) substr(s, pos + 1L, n) else ""
  paste0(pre, newchar, post)
}

#' Validate DNA (ACGTN only), returns uppercase no-whitespace
#' @keywords internal
gd_validate_dna <- function(s, allow_n = TRUE) {
  if (!is.character(s) || length(s) != 1L) cli::cli_abort("{.arg s} must be length-one character")
  s2 <- toupper(gsub("\\s+", "", s, perl = TRUE))
  pat <- if (allow_n) "^[ACGTN]*$" else "^[ACGT]*$"
  if (!grepl(pat, s2, perl = TRUE)) cli::cli_abort("Non-DNA characters detected in sequence")
  s2
}


# ---- AA / codon utilities ----

#' Minimal DNA -> AA translation (standard code, '*' for stop, 'X' for ambiguity)
#' Assumes length multiple of 3.
#' @keywords internal
gd_translate_cds <- function(dna) {
  dna <- gd_validate_dna(dna, allow_n = TRUE)
  n <- nchar(dna)
  if (n %% 3L != 0L) cli::cli_abort("CDS length must be divisible by 3 (got {n})")
  codons <- substring(dna, seq(1L, n - 2L, by = 3L), seq(3L, n, by = 3L))
  tbl <- c(
    TTT="F", TTC="F", TTA="L", TTG="L", TCT="S", TCC="S", TCA="S", TCG="S",
    TAT="Y", TAC="Y", TAA="*", TAG="*", TGT="C", TGC="C", TGA="*", TGG="W",
    CTT="L", CTC="L", CTA="L", CTG="L", CCT="P", CCC="P", CCA="P", CCG="P",
    CAT="H", CAC="H", CAA="Q", CAG="Q", CGT="R", CGC="R", CGA="R", CGG="R",
    ATT="I", ATC="I", ATA="I", ATG="M", ACT="T", ACC="T", ACA="T", ACG="T",
    AAT="N", AAC="N", AAA="K", AAG="K", AGT="S", AGC="S", AGA="R", AGG="R",
    GTT="V", GTC="V", GTA="V", GTG="V", GCT="A", GCC="A", GCA="A", GCG="A",
    GAT="D", GAC="D", GAA="E", GAG="E", GGT="G", GGC="G", GGA="G", GGG="G"
  )
  is_amb <- grepl("N", codons, fixed = TRUE)
  aa <- tbl[codons]
  aa[is.na(aa)] <- "X"
  aa[is_amb]    <- "X"
  paste0(unname(aa), collapse = "")
}

#' Codon window and AA index from CDS position (1-based)
#' @keywords internal
gd_codon_window <- function(cds_pos) {
  stopifnot(length(cds_pos) == 1L, is.finite(cds_pos), cds_pos >= 1L)
  codon_start <- ((cds_pos - 1L) %/% 3L) * 3L + 1L
  aa_index    <- ((cds_pos - 1L) %/% 3L) + 1L
  list(codon_start = codon_start, aa_index = aa_index)
}

# ---- light formatting ----

#' Format integer with thousands separator
#' @keywords internal
gd_fmt_int <- function(x) {
  ifelse(is.na(x), NA_character_, formatC(as.integer(x), big.mark = ",", format = "d"))
}

#' Format frequency (0 to 1) as percentage with 1 decimal
#' @keywords internal
gd_fmt_pct <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.1f%%", 100 * as.numeric(x)))
}

#' Normalize arrows / slashes spacing in gene/product fields
#' @keywords internal
gd_arrowize <- function(x) {
  if (is.null(x)) return(x)
  x <- gsub("\\s+", " ", x, perl = TRUE)
  x <- gsub("/", " / ", x, perl = TRUE)
  trimws(x)
}

# ---- class guard 

#' Assert that x is a genome_entity_gd
#' @keywords internal
gd_assert <- function(x, arg = "x") {
  if (!inherits(x, "genome_entity_gd")) {
    cli::cli_abort("{arg} must be a 'genome_entity_gd' (got class: {paste(class(x), collapse = '/')})")
  }
  invisible(x)
}


# ---- Vector-safe predicates for presenters/assemblers ----

#' TRUE for every element that is not NA
#' @keywords internal
v_not_na <- function(x) {
  !is.na(x)
}

#' TRUE for every element that is non-NA and non-empty string
#' @keywords internal
v_nzchar <- function(x) {
  x <- as.character(x)
  !is.na(x) & nzchar(x)
}

#' TRUE for every element that is numeric and greater than zero
#' @keywords internal
v_pos_gt0 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  !is.na(x) & (x > 0)
}


# ---- New tooling ----


# Compute the codon/AA effect of a single-base substitution at a locus for a given gene
# Evidence-agnostic: will look up RA to get ref/alt; for future extension, allow passing ref/alt.
#' @keywords internal
gd_compute_snp_effect <- function(gd, seq_id, position, gene) {
  gd_assert(gd, "gd")
  # Normalize inputs
  seq_id   <- as.character(seq_id)
  position <- as.integer(gsub("[^0-9]", "", position))
  
  # 1) Find the RA evidence at this locus to obtain alleles
  ev_tbl <- gd_events_table(gd, kinds = "evidence", expand_tags = FALSE)
  ra <- ev_tbl[ev_tbl$type == "RA" & ev_tbl$seq_id == seq_id & ev_tbl$position == position, , drop = FALSE]
  if (nrow(ra) != 1L) {
    return(list(
      ok = FALSE, reason = sprintf("expected 1 RA at %s:%d; found %d", seq_id, position, nrow(ra))
    ))
  }
  ref <- ra$ra_ref_base
  alt <- ra$ra_new_base
  if (is.na(ref) || is.na(alt) || nchar(ref) != 1L || nchar(alt) != 1L) {
    return(list(ok = FALSE, reason = "non-single-base alleles or missing ref/alt"))
  }
  
  # 2) Resolve feature/gene and map genomic -> CDS
  feat <- gd_resolve_feature(gd, gene)  # prefers CDS; falls back to gene
  strand <- feat$strand[1L]
  
  cds_pos <- map_genomic_to_cds(gd, gene, position)
  cds_pos <- if (is.list(cds_pos)) scalar_or(cds_pos$cds_pos, NA_integer_) else cds_pos
  if (is.na(cds_pos)) {
    return(list(ok = FALSE, reason = "position_not_in_CDS"))
  }
  
  # 3) Retrieve reference CDS and AA strings
  dna_ref <- get_gene_dna(gd, gene)
  aa_ref  <- get_gene_aa(gd, gene)
  if (is.na(dna_ref) || nchar(dna_ref) < 3L) {
    return(list(ok = FALSE, reason = "invalid_reference_CDS"))
  }
  
  # 4) Strand-aware allele to splice into CDS
  # For minus-strand genes, the CDS is oriented 5'->3' in gene direction, so complement the ALT.
  alt_cds <- if (identical(strand, "+")) alt else gd_complement_base(alt)
  
  # 5) Edit the CDS at that CDS position
  dna_mut <- gd_replace_substring1(dna_ref, cds_pos, alt_cds)
  
  # 6) Identify codon window and translate
  w <- gd_codon_window(cds_pos)
  codon_ref <- toupper(substr(dna_ref, w$codon_start, w$codon_start + 2L))
  codon_new <- toupper(substr(dna_mut, w$codon_start, w$codon_start + 2L))
  
  aa_mut <- gd_translate_cds(dna_mut)
  aa_ref_i <- substr(aa_ref, w$aa_index, w$aa_index)
  aa_new_i <- substr(aa_mut, w$aa_index, w$aa_index)
  
  consequence <- if (identical(aa_ref_i, aa_new_i)) "synonymous" else if (identical(aa_new_i, "*")) "nonsense" else "missense"
  
  list(
    ok         = TRUE,
    gene       = gene,
    strand     = strand,
    genomic    = list(seq_id = seq_id, position = position, ref = ref, alt = alt),
    cds_pos    = as.integer(cds_pos),
    aa_index   = as.integer(w$aa_index),
    codon_ref  = codon_ref,
    codon_new  = codon_new,
    aa_ref     = aa_ref_i,
    aa_new     = aa_new_i,
    consequence = consequence
  )
}

# printer for gd_compute_snp_effect()
#' @keywords internal
gd_snp_effect_row <- function(gd, seq_id, position, gene) {
  r <- gd_compute_snp_effect(gd, seq_id, position, gene)
  if (!isTRUE(r$ok)) {
    return(data.frame(
      seq_id = as.character(seq_id),
      position_num = as.integer(gsub("[^0-9]", "", position)),
      gene   = as.character(gene),
      verified = FALSE,
      consequence = NA_character_,
      aa_change  = NA_character_,
      codon_change = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    seq_id = r$genomic$seq_id,
    position_num = r$genomic$position,
    gene   = r$gene,
    verified = TRUE,
    consequence = r$consequence,
    aa_change   = paste0(r$aa_ref, r$aa_index, r$aa_new),
    codon_change= paste0(r$codon_ref, "→", r$codon_new),
    stringsAsFactors = FALSE
  )
}

# Prefer CDS exact symbol; fall back to gene exact symbol.
#' @keywords internal
gd_resolve_feature <- function(gd, symbol) {
  gd_assert(gd, "gd")
  # Fast path: O(1) lookup via cds_hash (prefers CDS, built at load time)
  h <- gd$indices$cds_hash
  if (!is.null(h) && exists(symbol, envir = h, inherits = FALSE)) {
    row_idx <- get(symbol, envir = h, inherits = FALSE)
    return(gd$features[row_idx, , drop = FALSE])
  }
  # Slow path: regex search (fallback for objects without cds_hash)
  cds <- try(search_features(gd, type = "CDS", pattern = paste0("^", symbol, "$")), silent = TRUE)
  if (!inherits(cds, "try-error") && !is.null(cds) && nrow(cds) >= 1L) {
    return(cds[1L, , drop = FALSE])
  }
  gen <- try(search_features(gd, type = "gene", pattern = paste0("^", symbol, "$")), silent = TRUE)
  if (!inherits(gen, "try-error") && !is.null(gen) && nrow(gen) >= 1L) {
    return(gen[1L, , drop = FALSE])
  }
  stop("gd_resolve_feature: could not resolve feature for symbol: ", symbol, call. = FALSE)
}

#' @keywords internal
gd_classify_snp <- function(gd, seq_id, position, gene) {
  v <- gd_verify_snp(gd, seq_id, position, gene)
  if (!isTRUE(v$ok)) {
    return(c(v, consequence = NA_character_))
  }
  cons <- if (identical(v$aa_ref, v$aa_new)) "synonymous" else if (identical(v$aa_new, "*")) "nonsense" else "missense"
  c(v, consequence = cons)
}


#' Verify a single SNP by editing the CDS and re-translating
#'
#' @param gd genome_entity_gd
#' @param seq_id character or integer contig id (as used in GD)
#' @param position integer 1-based genomic position
#' @param gene feature symbol to verify against (prefer CDS symbol)
#' @return tibble-like list/data.frame with fields:
#'   ok, reason, seq_id, position, gene, strand, cds_pos, aa_index,
#'   codon_ref, codon_new, aa_ref, aa_new, consequence
#' @keywords internal
gd_verify_snp <- function(gd, seq_id, position, gene) {
  gd_assert(gd, "gd")
  
  # Evidence table – we still key off RA to get ref/alt bases
  ev_tbl <- gd_events_table(gd, kinds = "evidence", expand_tags = FALSE)
  
  # Normalize inputs
  seq_id_chr <- as.character(seq_id)
  pos_int    <- suppressWarnings(as.integer(gsub(",", "", position)))
  
  # Find the single RA evidence at (seq_id, position)
  ra <- subset(ev_tbl, type == "RA" & seq_id == seq_id_chr & position == pos_int)
  if (nrow(ra) != 1L) {
    return(tibble::tibble(
      ok = FALSE, reason = sprintf("expected 1 RA at %s:%s; found %d", seq_id_chr, pos_int, nrow(ra)),
      seq_id = seq_id_chr, position = pos_int, gene = gene,
      strand = NA_character_,
      cds_pos = NA_integer_, aa_index = NA_integer_,
      codon_ref = NA_character_, codon_new = NA_character_,
      aa_ref = NA_character_, aa_new = NA_character_,
      consequence = NA_character_
    ))
  }
  
  ref <- ra$ra_ref_base
  alt <- ra$ra_new_base
  
  # Resolve feature, mapping, and sequences
  feat   <- gd_resolve_feature(gd, gene)
  strand <- feat$strand[1L]
  
  cds_pos <- map_genomic_to_cds(gd, gene, pos_int)
  cds_pos <- if (is.list(cds_pos)) scalar_or(cds_pos$cds_pos, NA_integer_) else cds_pos
  if (is.na(cds_pos)) {
    return(tibble::tibble(
      ok = FALSE, reason = "noncoding_or_intergenic",
      seq_id = seq_id_chr, position = pos_int, gene = gene,
      strand = strand,
      cds_pos = NA_integer_, aa_index = NA_integer_,
      codon_ref = NA_character_, codon_new = NA_character_,
      aa_ref = NA_character_, aa_new = NA_character_,
      consequence = NA_character_
    ))
  }
  
  dna_ref <- get_gene_dna(gd, gene)
  aa_ref  <- get_gene_aa(gd, gene)
  
  # Strand discipline: complement allele for minus strand
  alt_cds <- if (identical(strand, "+")) alt else gd_revcomp(alt)
  
  # Edit the CDS at cds_pos (1-based)
  dna_mut <- gd_replace_substring1(dna_ref, cds_pos, alt_cds)
  
  # Compute codon window and translate
  w <- gd_codon_window(cds_pos)
  codon_ref <- toupper(substr(dna_ref, w$codon_start, w$codon_start + 2L))
  codon_new <- toupper(substr(dna_mut, w$codon_start, w$codon_start + 2L))
  
  aa_mut <- gd_translate_cds(dna_mut)
  aa_r   <- substr(aa_ref, w$aa_index, w$aa_index)
  aa_n   <- substr(aa_mut, w$aa_index, w$aa_index)
  
  consequence <- if (identical(aa_r, aa_n)) {
    "synonymous"
  } else if (identical(aa_n, "*")) {
    "nonsense"
  } else {
    "missense"
  }
  
  tibble::tibble(
    ok = TRUE, reason = NA_character_,
    seq_id = seq_id_chr, position = pos_int, gene = gene, strand = strand,
    cds_pos = as.integer(cds_pos), aa_index = as.integer(w$aa_index),
    codon_ref = codon_ref, codon_new = codon_new,
    aa_ref = aa_r, aa_new = aa_n,
    consequence = consequence
  )
}


# Return all CDS overlapping pos; if none, return nearest upstream/downstream CDS with distances
# GENOME coordinate geometry
#' @keywords internal
gd_features_at <- function(gd, seq_id, pos, feature_type = "CDS", 
                           max_candidates = 10L) {
  gd_assert(gd, "gd")
  # Pull CDS catalogue once (assumes you have a cached feature table accessor)
  # If you already have a function to list features, use that; otherwise:
  feats <- search_features(gd, type = feature_type, pattern = ".*")  # all CDS
  if (!nrow(feats)) return(list(overlap = feats[0,], upstream = NULL, downstream = NULL))

  # Restrict to contig
  f <- feats[feats$seqname == as.character(seq_id), , drop = FALSE]
  if (!nrow(f)) return(list(overlap = f[0,], upstream = NULL, downstream = NULL))
  
  # Normalize coordinates: start <= end; get overlaps
  s <- pmin(f$start, f$end); e <- pmax(f$start, f$end)
  overlaps <- (pos >= s) & (pos <= e)
  
  if (any(overlaps)) {
    return(list(
      overlap   = f[overlaps, , drop = FALSE],
      upstream  = NULL,
      downstream= NULL
    ))
  }
  
  # No overlap: compute signed distances to interval edges
  # distance is min(|pos - s|, |pos - e|) but signed by direction w.r.t. feature
  # We'll also retain which side (upstream/downstream) in genomic order.
  d_left  <- pos - e   # negative if feature end is to the right of pos
  d_right <- s - pos   # positive if feature start is to the right of pos
  
  # Upstream candidate: features entirely before pos => e < pos (so d_left < 0); pick max e (closest)
  upstream_idx   <- which(e < pos)
  downstream_idx <- which(s > pos)
  
  upstream <- if (length(upstream_idx)) {
    j <- upstream_idx[which.max(e[upstream_idx])]
    uf <- f[j, , drop = FALSE]
    uf$distance <- pos - e[j]  # positive distance to the left edge (in bp)
    uf
  } else NULL
  
  downstream <- if (length(downstream_idx)) {
    j <- downstream_idx[which.min(s[downstream_idx])]
    df <- f[j, , drop = FALSE]
    df$distance <- s[j] - pos  # positive distance to the right edge (in bp)
    df
  } else NULL
  
  list(overlap = f[0, , drop = FALSE], upstream = upstream, downstream = downstream)
}


# Classify a genomic position as coding/intergenic and emit a breseq-like label
#' @keywords internal
gd_locate <- function(gd, seq_id, pos) {
  gd_assert(gd, "gd")
  
  hit <- gd_features_at(gd, seq_id, pos, feature_type = "CDS")
  
  # If overlapping one or more CDS: report coding with per-CDS local coords
  if (nrow(hit$overlap)) {
    # For each overlapping CDS, compute CDS-relative position (1-based)
    ann <- character(nrow(hit$overlap))
    genes <- character(nrow(hit$overlap))

    for (i in seq_len(nrow(hit$overlap))) {
      # Try multiple column names for gene identifier (GFF sources vary)
      gene <- if ("gene" %in% names(hit$overlap)) hit$overlap$gene[i]
              else if ("Name" %in% names(hit$overlap)) hit$overlap$Name[i]
              else if ("locus_tag" %in% names(hit$overlap)) hit$overlap$locus_tag[i]
              else if ("ID" %in% names(hit$overlap)) hit$overlap$ID[i]
              else NA_character_
      strd <- hit$overlap$strand[i]
      gstart <- hit$overlap$start[i]; gend <- hit$overlap$end[i]
      # Normalize to CDS order (strand-aware)
      if (identical(strd, "+") || isTRUE(strd == 1L)) {
        cds_pos <- pos - min(gstart, gend) + 1L
      } else {
        cds_pos <- max(gstart, gend) - pos + 1L
      }
      cds_len <- abs(gend - gstart) + 1L
      ann[i] <- sprintf("coding (%d/%d nt)", cds_pos, cds_len)
      genes[i] <- gene
    }

    return(list(
      region = "coding",
      label  = paste(ann, collapse = " | "),
      genes  = paste(genes, collapse = " | ")
    ))
  }
  
  # Intergenic: compute signed distances to nearest upstream & downstream genes
  up <- hit$upstream
  dn <- hit$downstream
  
  # If no CDS on either side (weird contig ends), degrade gracefully
  if (is.null(up) && is.null(dn)) {
    return(list(region = "intergenic", label = "intergenic (NA/NA)", genes = NA_character_))
  }
  
  # Strand-aware sign convention:
  # If nearest gene is '+' strand and mutation is downstream of its end -> '+<dist>'
  # If nearest gene is '-' strand and mutation is downstream of its start w.r.t. transcription -> '+<dist>'
  # For a succinct approximation, we’ll sign by *genomic* left/right, and append arrows to gene names.
  # This keeps display consistent with what you’re already emitting (gene arrows via strand).
  
  # Helper to extract gene name from feature row (tries multiple column names)
  .get_gene_name <- function(feat_row) {
    if ("gene" %in% names(feat_row)) feat_row$gene[1]
    else if ("Name" %in% names(feat_row)) feat_row$Name[1]
    else if ("locus_tag" %in% names(feat_row)) feat_row$locus_tag[1]
    else if ("ID" %in% names(feat_row)) feat_row$ID[1]
    else NA_character_
  }

  up_label <- NA_character_; dn_label <- NA_character_
  up_gene  <- NA_character_; dn_gene  <- NA_character_

  if (!is.null(up)) {
    # pos is to the right of 'up' gene; distance stored as positive
    up_label <- paste0("+", up$distance[1])
    up_gene  <- .get_gene_name(up)
    if (!is.null(up$strand[1]) && up$strand[1] %in% c("-", -1L)) up_gene <- paste0(up_gene, " \u2190") else up_gene <- paste0(up_gene, " \u2192")
  }
  if (!is.null(dn)) {
    # pos is to the left of 'dn' gene; distance stored as positive
    dn_label <- paste0("+", dn$distance[1])
    dn_gene  <- .get_gene_name(dn)
    if (!is.null(dn$strand[1]) && dn$strand[1] %in% c("-", -1L)) dn_gene <- paste0(dn_gene, " \u2190") else dn_gene <- paste0(dn_gene, " \u2192")
  }
  
  # Compose breseq-like “intergenic (+a/+b)” and gene pair “Gleft / Gright”
  label <- sprintf("intergenic (%s/%s)",
                   ifelse(is.na(up_label), "NA", up_label),
                   ifelse(is.na(dn_label), "NA", dn_label))
  genes <- paste(na.omit(c(up_gene, dn_gene)), collapse = " / ")
  if (!nzchar(genes)) genes <- NA_character_
  
  list(region = "intergenic", label = label, genes = genes)
}


# Robust integer parse for positions that may contain commas
#' @keywords internal
.pm_parse_pos1 <- function(x) {
  if (is.numeric(x)) return(as.integer(x))
  x <- as.character(x)
  x <- gsub("[^0-9-]", "", x, perl = TRUE)
  suppressWarnings(as.integer(x))
}


# Defensive RA allele accessor: handles column name variants you saw in practice
#' @keywords internal
gd_ra_alleles <- function(ra_row) {
  cols <- names(ra_row)
  pick <- function(primary, alt) {
    if (primary %in% cols) ra_row[[primary]] else if (alt %in% cols) ra_row[[alt]] else NA_character_
  }
  ref <- pick("ra_ref_base",  "snp_ref_base")
  alt <- pick("ra_new_base",  "snp_alt_base")
  list(ref = as.character(ref)[1], alt = as.character(alt)[1])
}

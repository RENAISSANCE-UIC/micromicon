#' Get DNA Sequence for a Gene
#'
#' @description
#' Extract the DNA sequence for a specific gene feature in 5'-to-3' canonical orientation.
#' The sequence is automatically reverse-complemented if the gene is on the minus strand.
#'
#' @param x A genome_entity object
#' @param gene Gene identifier: gene name, locus_tag, feature index (integer), or feature row (data.frame)
#' @param ... Additional arguments
#'
#' @return Character vector (DNA sequence) or DNAString if Biostrings is available
#' @export
get_gene_dna <- function(x, gene, ...) {
  UseMethod("get_gene_dna")
}

#' @export
get_gene_dna.genome_entity <- function(x, gene, format = c("character", "DNAString"), ...) {
  format <- match.arg(format)
  validate_genome_entity(x)

  # Resolve gene to feature row
  feat <- .resolve_gene(x, gene)

  # Extract DNA sequence
  full_seq <- x$sequences$dna_raw[[feat$seqname]]

  if (is.null(full_seq)) {
    cli::cli_abort("Sequence '{feat$seqname}' not found in genome_entity")
  }

  # Validate coordinates
  if (feat$start < 1 || feat$end > nchar(full_seq)) {
    cli::cli_abort("Gene coordinates out of bounds: {feat$start}-{feat$end}")
  }

  # Extract sequence, using precomputed reverse complement for minus-strand genes
  # when available (avoids per-call reverse_complement() on potentially long genes).
  if (!is.na(feat$strand) && feat$strand == "-") {
    dna_rev <- x$sequences$dna_rev[[feat$seqname]]
    if (!is.null(dna_rev)) {
      L <- nchar(full_seq)
      dna <- substr(dna_rev, L - feat$end + 1L, L - feat$start + 1L)
    } else {
      dna <- reverse_complement(substr(full_seq, feat$start, feat$end))
    }
  } else {
    dna <- substr(full_seq, feat$start, feat$end)
  }

  if (format == "DNAString") {
    bio_gateway <- create_bioconductor_gateway()
    if (!bio_gateway$is_available()) {
      cli::cli_abort(c(
        "Biostrings required for DNAString format",
        "i" = "Install with: BiocManager::install('Biostrings')"
      ))
    }
    dna_set <- bio_gateway$to_dnastringset(c(gene = dna))
    return(dna_set[[1]])
  }

  dna
}

#' @export
get_gene_dna.default <- function(x, gene, ...) {
  cli::cli_abort("get_gene_dna() not implemented for class {.cls {class(x)[1]}}")
}


#' Get Amino Acid Sequence for a Gene
#'
#' @description
#' Extract and translate the amino acid sequence for a CDS feature.
#' Uses the genetic code specified in the feature's transl_table attribute,
#' or defaults to genetic code 11 (bacterial/archaeal).
#'
#' In bacterial genomes, alternative start codons (GTG, TTG, CTG, ATT, ATC, ATA)
#' are translated as Methionine (M) when they occur as the first codon, even though
#' they normally code for other amino acids (e.g., GTG codes for V internally, but M at start).
#'
#' @param x A genome_entity object
#' @param gene Gene identifier: gene name, locus_tag, feature index (integer), or feature row (data.frame)
#' @param genetic_code Genetic code to use (default = "11" for bacteria). Overrides feature transl_table.
#' @param fix_start_codon Logical, whether to fix alternative start codons to M (default = TRUE)
#' @param ... Additional arguments
#'
#' @return Character vector (amino acid sequence)
#' @export
get_gene_aa <- function(x, gene, ...) {
  UseMethod("get_gene_aa")
}

#' @rdname get_gene_aa
#' @export
get_gene_aa.genome_entity <- function(x, gene, genetic_code = NULL, fix_start_codon = TRUE, ...) {
  validate_genome_entity(x)

  # Resolve gene to feature row
  feat <- .resolve_gene(x, gene)

  # Check if it's a CDS
  if (is.na(feat$type) || feat$type != "CDS") {
    cli::cli_warn("Feature is not a CDS (type = '{feat$type}'). Translation may be incorrect.")
  }

  # Fast path: use pre-stored GenBank /translation= when available and no
  # caller overrides are in play (genetic_code = NULL, fix_start_codon = TRUE).
  # GenBank translations are authoritative: already M-normalised, correct
  # genetic code, and exclude the stop codon.
  if (is.null(genetic_code) && isTRUE(fix_start_codon) &&
      "translation" %in% names(feat) &&
      !is.na(feat$translation) && nzchar(feat$translation)) {
    return(as.character(feat$translation))
  }

  # Slow path: extract DNA and translate
  dna <- get_gene_dna(x, gene, format = "character")

  # Determine genetic code
  if (is.null(genetic_code)) {
    if ("transl_table" %in% names(feat) && !is.na(feat$transl_table)) {
      genetic_code <- as.character(feat$transl_table)
    } else {
      genetic_code <- "11"  # Default: bacterial/archaeal
    }
  }

  aa <- translate_dna(dna, frame = 1, genetic_code = genetic_code, .internal = TRUE)

  # Fix alternative start codons (GTG, TTG, CTG -> M in bacteria)
  if (fix_start_codon && nchar(dna) >= 3) {
    start_codon <- toupper(substr(dna, 1, 3))
    alt_starts <- c("GTG", "TTG", "CTG", "ATT", "ATC", "ATA")
    if (start_codon %in% alt_starts && nchar(aa) > 0) {
      substr(aa, 1, 1) <- "M"
    }
  }

  aa
}

#' @export
get_gene_aa.default <- function(x, gene, ...) {
  cli::cli_abort("get_gene_aa() not implemented for class {.cls {class(x)[1]}}")
}


#' Get DNA Sequence for a Region of Interest (ROI)
#'
#' @description
#' Extract DNA sequence for an arbitrary genomic region.
#' This is a convenience wrapper around extract_by_coords.
#'
#' @param x A genome_entity object
#' @param contig Contig/chromosome name
#' @param start Start position (1-based, inclusive)
#' @param end End position (1-based, inclusive)
#' @param strand Strand ("+" or "-", default = "+")
#' @param ... Additional arguments
#'
#' @return Character vector (DNA sequence)
#' @export
get_roi_dna <- function(x, contig, start, end, strand = "+", ...) {
  UseMethod("get_roi_dna")
}

#' @export
get_roi_dna.genome_entity <- function(x, contig, start, end, strand = "+", ...) {
  validate_genome_entity(x)

  # Use existing extract_by_coords generic
  result <- extract_by_coords(x, contig = contig, start = start, end = end,
                              strand = strand, translate = FALSE)

  # Return just the sequence (strip names)
  as.character(result)
}

#' @export
get_roi_dna.default <- function(x, contig, start, end, strand = "+", ...) {
  cli::cli_abort("get_roi_dna() not implemented for class {.cls {class(x)[1]}}")
}


#' Get DNA Sequences for Multiple Regions of Interest (vectorised)
#'
#' @description
#' Extract DNA sequences for a batch of genomic regions in a single call.
#' Rows sharing the same contig are processed together with a single
#' \code{substring()} call, making this much faster than calling
#' \code{get_roi_dna()} in a loop.  Minus-strand regions use the precomputed
#' \code{sequences$dna_rev} when available (built by \code{new_genome_entity}).
#'
#' @param x A genome_entity object
#' @param contig Character vector of contig/chromosome names (length n)
#' @param start Integer vector of start positions, 1-based inclusive (length n)
#' @param end Integer vector of end positions, 1-based inclusive (length n)
#' @param strand Character vector of strand values ("+" or "-"); recycled to
#'   length n (default "+")
#' @param ... Additional arguments (currently unused)
#'
#' @return Character vector of length n
#' @keywords internal
get_roi_dna_vec <- function(x, contig, start, end, strand = "+", ...) {
  UseMethod("get_roi_dna_vec")
}

#' @keywords internal
get_roi_dna_vec.genome_entity <- function(x, contig, start, end, strand = "+", ...) {
  validate_genome_entity(x)

  n <- length(contig)
  result <- character(n)
  strand <- rep_len(strand, n)

  plus_idx  <- which(strand != "-")
  minus_idx <- which(strand == "-")

  # Plus strand: group by contig, one substring() call per contig
  if (length(plus_idx) > 0) {
    for (sn in unique(contig[plus_idx])) {
      rows <- plus_idx[contig[plus_idx] == sn]
      seq  <- x$sequences$dna_raw[[sn]]
      result[rows] <- if (is.null(seq)) NA_character_ else
        substring(seq, start[rows], end[rows])
    }
  }

  # Minus strand: use precomputed dna_rev for vectorised slicing when available
  if (length(minus_idx) > 0) {
    for (sn in unique(contig[minus_idx])) {
      rows    <- minus_idx[contig[minus_idx] == sn]
      seq     <- x$sequences$dna_raw[[sn]]
      if (is.null(seq)) {
        result[rows] <- NA_character_
      } else {
        dna_rev <- x$sequences$dna_rev[[sn]]
        if (!is.null(dna_rev)) {
          L <- nchar(seq)
          result[rows] <- substring(dna_rev, L - end[rows] + 1L, L - start[rows] + 1L)
        } else {
          for (i in rows) {
            result[i] <- reverse_complement(substr(seq, start[i], end[i]))
          }
        }
      }
    }
  }

  result
}

#' @keywords internal
get_roi_dna_vec.default <- function(x, contig, start, end, strand = "+", ...) {
  cli::cli_abort("get_roi_dna_vec() not implemented for class {.cls {class(x)[1]}}")
}


#' Translate DNA to Protein
#'
#' @description
#' Translate DNA sequence to amino acid sequence using specified genetic code and reading frame.
#' This function performs standard genetic code translation without special handling
#' of alternative start codons.
#'
#' **Note:** If you are translating a complete CDS from a bacterial genome, use
#' \code{\link{get_gene_aa}} instead, which correctly handles alternative start codons
#' (GTG, TTG, CTG, etc.) that translate to Methionine when used as start codons.
#'
#' @param dna Character string of DNA sequence
#' @param frame Reading frame (1, 2, or 3, default = 1)
#' @param genetic_code Genetic code to use (default = "11" for bacterial/archaeal).
#'   Can be NCBI genetic code number ("1", "11", etc.) or a named codon table vector.
#'
#' @return Protein sequence (single-letter amino acids)
#' @param .internal Internal use only - suppress informational message
#' @seealso \code{\link{get_gene_aa}} for CDS translation with alternative start codon handling
#' @keywords internal
translate_dna <- function(dna, frame = 1, genetic_code = "11", .internal = FALSE) {
  # One-time informational message per session (only for user-facing calls)
  if (!.internal && !isTRUE(getOption("micromicon.translate_dna.informed"))) {
    cli::cli_inform(c(
      "i" = "{.fn translate_dna} uses standard genetic code translation",
      "i" = "Alternative start codons (GTG, TTG, etc.) are NOT converted to M",
      ">" = "For complete CDS translation, use {.fn get_gene_aa} instead"
    ))
    options(micromicon.translate_dna.informed = TRUE)
  }

  # Validate frame
  if (!frame %in% 1:3) {
    cli::cli_abort("frame must be 1, 2, or 3")
  }

  # Adjust for frame
  if (frame > 1) {
    dna <- substr(dna, frame, nchar(dna))
  }

  # Get codon table
  if (is.character(genetic_code) && length(genetic_code) == 1) {
    # NCBI genetic code number
    codon_table <- .get_genetic_code_table(genetic_code)
  } else if (is.character(genetic_code) && length(genetic_code) > 1) {
    # Custom codon table
    codon_table <- genetic_code
  } else {
    cli::cli_abort("genetic_code must be a string (NCBI code) or named character vector")
  }

  # Convert to uppercase
  dna <- toupper(dna)

  # Extract codons (groups of 3)
  n_codons <- floor(nchar(dna) / 3)
  if (n_codons == 0) return("")

  starts <- seq(1L, n_codons * 3L, by = 3L)
  codons <- substring(dna, starts, starts + 2L)

  # Translate
  aa <- codon_table[codons]

  # Handle unknown codons
  aa[is.na(aa)] <- "X"

  paste(aa, collapse = "")
}


#' Map Genomic Position to CDS Position
#'
#' @description
#' Convert a genomic position to a CDS position (1-based) for a gene.
#' Accounts for exons and strand orientation.
#'
#' @param x A genome_entity object
#' @param gene Gene identifier: gene name, locus_tag, feature index (integer), or feature row (data.frame)
#' @param genomic_pos Genomic position (1-based)
#' @param ... Additional arguments
#'
#' @return Integer CDS position (1-based), or NA if position is not in CDS
#' @export
map_genomic_to_cds <- function(x, gene, genomic_pos, ...) {
  UseMethod("map_genomic_to_cds")
}

#' @export
map_genomic_to_cds.genome_entity <- function(x, gene, genomic_pos, ...) {
  validate_genome_entity(x)

  # Resolve gene to feature row
  feat <- .resolve_gene(x, gene)

  # Validate genomic position is on same seqname
  # (In more complex cases with exons, you'd need to check exon coordinates)

  # For now, assume simple CDS without introns
  if (genomic_pos < feat$start || genomic_pos > feat$end) {
    return(NA_integer_)
  }

  # Calculate CDS position based on strand
  if (!is.na(feat$strand) && feat$strand == "-") {
    # Minus strand: count from end
    cds_pos <- feat$end - genomic_pos + 1L
  } else {
    # Plus strand: count from start
    cds_pos <- genomic_pos - feat$start + 1L
  }

  cds_pos
}

#' @export
map_genomic_to_cds.default <- function(x, gene, genomic_pos, ...) {
  cli::cli_abort("map_genomic_to_cds() not implemented for class {.cls {class(x)[1]}}")
}


#' Map CDS Position to Genomic Position
#'
#' @description
#' Convert a CDS position (1-based) to a genomic position for a gene.
#' Accounts for strand orientation (inverse of map_genomic_to_cds).
#'
#' @param x A genome_entity object
#' @param gene Gene identifier: gene name, locus_tag, feature index (integer), or feature row (data.frame)
#' @param cds_pos CDS position (1-based)
#' @param ... Additional arguments
#'
#' @return Integer genomic position (1-based), or NA if CDS position is out of range
#' @export
map_cds_to_genomic <- function(x, gene, cds_pos, ...) {
  UseMethod("map_cds_to_genomic")
}

#' @export
map_cds_to_genomic.genome_entity <- function(x, gene, cds_pos, ...) {
  validate_genome_entity(x)

  # Resolve gene to feature row
  feat <- .resolve_gene(x, gene)

  # Calculate gene length
  gene_len <- feat$end - feat$start + 1L

  # Validate CDS position
  if (cds_pos < 1L || cds_pos > gene_len) {
    return(NA_integer_)
  }

  # Calculate genomic position based on strand
  if (!is.na(feat$strand) && feat$strand == "-") {
    # Minus strand: count from end
    genomic_pos <- feat$end - cds_pos + 1L
  } else {
    # Plus strand: count from start
    genomic_pos <- feat$start + cds_pos - 1L
  }

  genomic_pos
}

#' @export
map_cds_to_genomic.default <- function(x, gene, cds_pos, ...) {
  cli::cli_abort("map_cds_to_genomic() not implemented for class {.cls {class(x)[1]}}")
}


#' Get Gene Information
#'
#' @description
#' Extract metadata for a gene feature including chromosome, coordinates, strand,
#' CDS ranges, and reading frame information.
#'
#' @param x A genome_entity object
#' @param gene Gene identifier: gene name, locus_tag, feature index (integer), or feature row (data.frame)
#' @param ... Additional arguments
#'
#' @return List with gene information:
#'   - chrom: Chromosome/contig name
#'   - start: Start position (1-based)
#'   - end: End position (1-based)
#'   - strand: Strand ("+" or "-")
#'   - type: Feature type (e.g., "CDS", "gene")
#'   - length: Gene length in bp
#'   - cds_ranges: List of CDS ranges (for now, just start-end)
#'   - frame: Reading frame offset (0, 1, or 2)
#'   - gene: Gene name (if available)
#'   - locus_tag: Locus tag (if available)
#'   - product: Product description (if available)
#' @export
gene_info <- function(x, gene, ...) {
  UseMethod("gene_info")
}

#' @export
gene_info.genome_entity <- function(x, gene, ...) {
  validate_genome_entity(x)

  # Resolve gene to feature row
  feat <- .resolve_gene(x, gene)

  # Extract frame if available (GFF3 phase field)
  frame <- if ("phase" %in% names(feat) && !is.na(feat$phase)) {
    as.integer(feat$phase)
  } else {
    0L
  }

  # Build info list
  info <- list(
    contig = feat$seqname,
    start = feat$start,
    end = feat$end,
    strand = if (!is.na(feat$strand)) feat$strand else "+",
    type = if (!is.na(feat$type)) feat$type else NA_character_,
    length = feat$end - feat$start + 1L,
    cds_ranges = list(c(start = feat$start, end = feat$end)),
    frame = frame,
    gene = if ("gene" %in% names(feat)) feat$gene else NA_character_,
    locus_tag = if ("locus_tag" %in% names(feat)) feat$locus_tag else NA_character_,
    product = if ("product" %in% names(feat)) feat$product else NA_character_
  )

  info
}

#' @export
gene_info.default <- function(x, gene, ...) {
  cli::cli_abort("gene_info() not implemented for class {.cls {class(x)[1]}}")
}


#' Validate Variant in Gene
#'
#' @description
#' Check if a genomic variant (position + reference base) is valid within a gene.
#' Validates that the position is within the gene and the reference base matches.
#'
#' @param x A genome_entity object
#' @param gene Gene identifier: gene name, locus_tag, feature index (integer), or feature row (data.frame)
#' @param genomic_pos Genomic position (1-based)
#' @param ref_base Expected reference base (single character)
#' @param ... Additional arguments
#'
#' @return Logical value (TRUE if valid, FALSE if not) with message attribute
#' @export
validate_variant_in_gene <- function(x, gene, genomic_pos, ref_base, ...) {
  UseMethod("validate_variant_in_gene")
}

#' @export
validate_variant_in_gene.genome_entity <- function(x, gene, genomic_pos, ref_base, ...) {
  validate_genome_entity(x)

  # Resolve gene to feature row
  feat <- .resolve_gene(x, gene)

  # Check if position is within gene
  if (genomic_pos < feat$start || genomic_pos > feat$end) {
    result <- FALSE
    attr(result, "message") <- sprintf(
      "Position %d is outside gene bounds [%d, %d]",
      genomic_pos, feat$start, feat$end
    )
    return(result)
  }

  # Get the reference base at this position
  full_seq <- x$sequences$dna_raw[[feat$seqname]]
  actual_base <- substr(full_seq, genomic_pos, genomic_pos)

  # Handle strand
  if (!is.na(feat$strand) && feat$strand == "-") {
    # For minus strand genes, need to reverse complement the reference
    actual_base <- reverse_complement(actual_base)
  }

  # Compare bases (case-insensitive)
  ref_base <- toupper(ref_base)
  actual_base <- toupper(actual_base)

  if (actual_base != ref_base) {
    result <- FALSE
    attr(result, "message") <- sprintf(
      "Reference base mismatch at position %d: expected '%s', found '%s'",
      genomic_pos, ref_base, actual_base
    )
    return(result)
  }

  # Valid
  result <- TRUE
  attr(result, "message") <- "Valid variant"
  result
}

#' @export
validate_variant_in_gene.default <- function(x, gene, genomic_pos, ref_base, ...) {
  cli::cli_abort("validate_variant_in_gene() not implemented for class {.cls {class(x)[1]}}")
}


# ---- Internal Helper Functions ----

#' Resolve Gene Identifier to Feature Row
#'
#' @description
#' Internal helper to resolve various gene identifiers to a single feature row.
#' Search order for character identifiers: locus_tag, gene, Name, Alias, ID
#'
#' @param entity A genome_entity object
#' @param gene Gene identifier: gene name, locus_tag, Name, Alias, ID, feature index (integer), or feature row (data.frame)
#'
#' @return Single-row data.frame (feature)
#' @keywords internal
.resolve_gene <- function(entity, gene) {
  # Case 1: Integer index
  if (is.numeric(gene) && length(gene) == 1) {
    idx <- as.integer(gene)
    if (idx < 1 || idx > nrow(entity$features)) {
      cli::cli_abort("Feature index {idx} out of range [1, {nrow(entity$features)}]")
    }
    return(entity$features[idx, , drop = FALSE])
  }

  # Case 2: Data frame (single row)
  if (is.data.frame(gene)) {
    if (nrow(gene) != 1) {
      cli::cli_abort("gene data.frame must have exactly 1 row")
    }
    return(gene)
  }

  # Case 3: Character (gene name or locus_tag)
  if (is.character(gene) && length(gene) == 1) {
    # Fast path: O(1) lookup via cds_hash (built at load time)
    h <- entity$indices$cds_hash
    if (!is.null(h) && exists(gene, envir = h, inherits = FALSE)) {
      row_idx <- get(gene, envir = h, inherits = FALSE)
      return(entity$features[row_idx, , drop = FALSE])
    }

    # Slow path: sequential scan (fallback for objects without cds_hash)
    # Try locus_tag first (more specific)
    if ("locus_tag" %in% names(entity$features)) {
      matches <- which(entity$features$locus_tag == gene)
      if (length(matches) == 1) {
        return(entity$features[matches[1], , drop = FALSE])
      } else if (length(matches) > 1) {
        result <- .resolve_multiple_matches(entity$features[matches, , drop = FALSE], gene, "locus_tag")
        return(result)
      }
    }

    # Try gene name
    if ("gene" %in% names(entity$features)) {
      matches <- which(entity$features$gene == gene)
      if (length(matches) == 1) {
        return(entity$features[matches[1], , drop = FALSE])
      } else if (length(matches) > 1) {
        result <- .resolve_multiple_matches(entity$features[matches, , drop = FALSE], gene, "gene")
        return(result)
      }
    }

    # Try Name (GFF3 display name)
    if ("Name" %in% names(entity$features)) {
      matches <- which(entity$features$Name == gene)
      if (length(matches) == 1) {
        return(entity$features[matches[1], , drop = FALSE])
      } else if (length(matches) > 1) {
        result <- .resolve_multiple_matches(entity$features[matches, , drop = FALSE], gene, "Name")
        return(result)
      }
    }

    # Try Alias (alternative names, common in breseq)
    if ("Alias" %in% names(entity$features)) {
      matches <- which(entity$features$Alias == gene)
      if (length(matches) == 1) {
        return(entity$features[matches[1], , drop = FALSE])
      } else if (length(matches) > 1) {
        result <- .resolve_multiple_matches(entity$features[matches, , drop = FALSE], gene, "Alias")
        return(result)
      }
    }

    # Try ID
    if ("ID" %in% names(entity$features)) {
      matches <- which(entity$features$ID == gene)
      if (length(matches) == 1) {
        return(entity$features[matches[1], , drop = FALSE])
      } else if (length(matches) > 1) {
        result <- .resolve_multiple_matches(entity$features[matches, , drop = FALSE], gene, "ID")
        return(result)
      }
    }

    cli::cli_abort("Gene '{gene}' not found in features")
  }

  cli::cli_abort("Invalid gene identifier type: {class(gene)[1]}")
}


#' Resolve Multiple Feature Matches
#'
#' @description
#' Internal helper to intelligently handle multiple feature matches.
#' When multiple features match (e.g., gene + CDS with same name at same location),
#' this function prefers CDS features and only warns for genuinely ambiguous cases.
#'
#' @param matched_features Data.frame of features that matched the query
#' @param gene_id Gene identifier that was searched
#' @param field_name Name of field that was searched ("gene", "locus_tag", or "ID")
#'
#' @return Single-row data.frame (the selected feature)
#' @keywords internal
.resolve_multiple_matches <- function(matched_features, gene_id, field_name) {
  # Check if all matches are at the same genomic location
  if ("start" %in% names(matched_features) && "end" %in% names(matched_features)) {
    same_location <- length(unique(matched_features$start)) == 1 &&
                     length(unique(matched_features$end)) == 1

    if (same_location && "type" %in% names(matched_features)) {
      # Multiple features at same location - likely gene + CDS + mRNA etc.
      # Prefer CDS for protein-coding genes
      cds_matches <- matched_features[matched_features$type == "CDS", , drop = FALSE]

      if (nrow(cds_matches) > 0) {
        # Found CDS at this location - use it without warning
        if (nrow(cds_matches) == 1) {
          return(cds_matches[1, , drop = FALSE])
        } else {
          # Multiple CDS at same location - this is genuinely ambiguous
          cli::cli_warn(paste0(
            "Multiple CDS features with ", field_name, " '{gene_id}' at same location. ",
            "Using first match. Consider using locus_tag or feature index for specificity."
          ))
          return(cds_matches[1, , drop = FALSE])
        }
      } else {
        # No CDS, but multiple features at same location (e.g., gene + rRNA)
        # Use first match but inform user
        types <- paste(unique(matched_features$type), collapse = ", ")
        cli::cli_inform(paste0(
          "Found ", nrow(matched_features), " features with ", field_name, " '{gene_id}' ",
          "at same location (types: ", types, "). Using first match (",
          matched_features$type[1], ")."
        ))
        return(matched_features[1, , drop = FALSE])
      }
    }
  }

  # Genuinely different features (different locations or no location info)
  # This is ambiguous - warn and use first
  cli::cli_warn(paste0(
    "Multiple distinct features with ", field_name, " '{gene_id}'. ",
    "Using first match. Consider using locus_tag or feature index for specificity."
  ))
  return(matched_features[1, , drop = FALSE])
}


#' Get Genetic Code Codon Table
#'
#' @description
#' Internal helper to get NCBI genetic code codon table.
#'
#' @param code_num NCBI genetic code number as string
#'
#' @return Named character vector of codon translations
#' @keywords internal
.get_genetic_code_table <- function(code_num = "11") {
  # For now, implement standard code (1) and bacterial code (11)
  # In future, could add more codes as needed

  if (code_num == "1") {
    # Standard genetic code
    return(c(
      TTT = "F", TTC = "F", TTA = "L", TTG = "L",
      TCT = "S", TCC = "S", TCA = "S", TCG = "S",
      TAT = "Y", TAC = "Y", TAA = "*", TAG = "*",
      TGT = "C", TGC = "C", TGA = "*", TGG = "W",
      CTT = "L", CTC = "L", CTA = "L", CTG = "L",
      CCT = "P", CCC = "P", CCA = "P", CCG = "P",
      CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
      CGT = "R", CGC = "R", CGA = "R", CGG = "R",
      ATT = "I", ATC = "I", ATA = "I", ATG = "M",
      ACT = "T", ACC = "T", ACA = "T", ACG = "T",
      AAT = "N", AAC = "N", AAA = "K", AAG = "K",
      AGT = "S", AGC = "S", AGA = "R", AGG = "R",
      GTT = "V", GTC = "V", GTA = "V", GTG = "V",
      GCT = "A", GCC = "A", GCA = "A", GCG = "A",
      GAT = "D", GAC = "D", GAA = "E", GAG = "E",
      GGT = "G", GGC = "G", GGA = "G", GGG = "G"
    ))
  } else if (code_num == "11") {
    # Bacterial, Archaeal and Plant Plastid Code
    # Same as standard code except for stop codons
    return(c(
      TTT = "F", TTC = "F", TTA = "L", TTG = "L",
      TCT = "S", TCC = "S", TCA = "S", TCG = "S",
      TAT = "Y", TAC = "Y", TAA = "*", TAG = "*",
      TGT = "C", TGC = "C", TGA = "*", TGG = "W",
      CTT = "L", CTC = "L", CTA = "L", CTG = "L",
      CCT = "P", CCC = "P", CCA = "P", CCG = "P",
      CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
      CGT = "R", CGC = "R", CGA = "R", CGG = "R",
      ATT = "I", ATC = "I", ATA = "I", ATG = "M",
      ACT = "T", ACC = "T", ACA = "T", ACG = "T",
      AAT = "N", AAC = "N", AAA = "K", AAG = "K",
      AGT = "S", AGC = "S", AGA = "R", AGG = "R",
      GTT = "V", GTC = "V", GTA = "V", GTG = "V",
      GCT = "A", GCC = "A", GCA = "A", GCG = "A",
      GAT = "D", GAC = "D", GAA = "E", GAG = "E",
      GGT = "G", GGC = "G", GGA = "G", GGG = "G"
    ))
  } else {
    cli::cli_abort("Genetic code {code_num} not implemented. Currently supports: 1 (standard), 11 (bacterial)")
  }
}


#' Reverse Complement DNA Sequence
#'
#' @description
#' Reverse complement a DNA sequence. Handles standard bases (A, T, G, C) and
#' IUPAC ambiguity codes (N, R, Y, S, W, K, M, B, V, D, H).
#'
#' @param seq Character string of DNA sequence
#' @return Character string of reverse complement sequence
#' @keywords internal
reverse_complement <- function(seq) {
  # Complement lookup
  complement_map <- c(
    A = "T", T = "A", G = "C", C = "G",
    a = "t", t = "a", g = "c", c = "g",
    N = "N", n = "n",
    R = "Y", Y = "R", S = "S", W = "W", K = "M", M = "K",
    B = "V", V = "B", D = "H", H = "D"
  )

  # Split into characters
  chars <- strsplit(seq, "")[[1]]

  # Complement
  comp_chars <- complement_map[chars]

  # Handle unknown characters
  comp_chars[is.na(comp_chars)] <- "N"

  # Reverse
  paste(rev(comp_chars), collapse = "")
}


# Note: %||% operator is defined in R/operators.R

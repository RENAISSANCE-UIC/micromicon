#' Enrich Mutation Table with Consequences
#'
#' @description
#' Adds molecular consequence annotations to a `predict_mutations()` table.
#' For SNPs in coding regions, computes reference/alternate sequences and
#' consequence (synonymous/missense/nonsense). For intergenic SNPs, extracts
#' a flanking DNA window.
#'
#' @param gd A GenomeData object
#' @param pm_tbl A data.frame from `predict_mutations()`, containing at minimum:
#'   `type`, `seq_id`, `position`, `mutation`, `gene`
#' @param flank Integer, number of bases to extract upstream/downstream for
#'   intergenic regions (default: 50)
#' @param quiet Logical, suppress informational messages (default: FALSE)
#'
#' @return A data.frame with additional columns:
#' \itemize{
#'   \item \code{dna_ref} - Reference DNA sequence (full CDS for coding; window for intergenic)
#'   \item \code{dna_alt} - Alternate DNA sequence with mutation applied
#'   \item \code{aa_ref} - Reference amino acid sequence (coding only)
#'   \item \code{aa_alt} - Alternate amino acid sequence (coding only)
#'   \item \code{codon_ref} - Reference codon at mutation site (coding only)
#'   \item \code{codon_alt} - Alternate codon at mutation site (coding only)
#'   \item \code{codon_new} - Alias for codon_alt (for consistency with gd_verify_snp)
#'   \item \code{consequence} - Mutation effect: "synonymous", "missense", "nonsense" (coding only)
#'   \item \code{region} - "coding" or "intergenic"
#'   \item \code{strand} - Gene strand (coding only)
#'   \item \code{qc_note} - Quality control notes (offset warnings, fallbacks, normalizations)
#' }
#'
#' @details
#' This function leverages existing micromicon utilities:
#' \itemize{
#'   \item \code{gd_verify_snp()} - Computes SNP effects for coding regions
#'   \item \code{gd_locate()} - Classifies regions as coding vs intergenic
#'   \item \code{get_gene_dna()}, \code{get_gene_aa()} - Extract gene sequences
#'   \item \code{get_roi_dna()} - Extract arbitrary genomic windows
#'   \item \code{translate_dna()} - Translate CDS to protein
#' }
#'
#' Currently supports SNP mutations only. Future versions will handle DEL, INS, SUB.
#'
#' @examples
#' \dontrun{
#' gd <- read_genomedata("data/sample.gd", ref_dir = "data/reference")
#' mutations <- predict_mutations(gd)
#' enriched <- pm_enrich_consequences(gd, mutations, flank = 50)
#' }
#'
#' @export
pm_enrich_consequences <- function(gd, pm_tbl, flank = 50L, quiet = FALSE) {
  gd_assert(gd, "gd")
  stopifnot(is.data.frame(pm_tbl))
  stopifnot(is.numeric(flank) && length(flank) == 1L && flank >= 0)

  # Validate required columns
  required_cols <- c("type", "seq_id", "position", "mutation", "gene")
  missing_cols <- setdiff(required_cols, names(pm_tbl))
  if (length(missing_cols) > 0) {
    stop(sprintf("pm_tbl missing required columns: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  if (!quiet) {
    cli::cli_inform(c(
      "i" = "pm_enrich_consequences(): enriching {nrow(pm_tbl)} mutation{?s}"
    ))
  }

  # Add new columns (handle empty table case)
  out <- pm_tbl
  n_rows <- nrow(out)
  out$dna_ref <- rep(NA_character_, n_rows)
  out$dna_alt <- rep(NA_character_, n_rows)
  out$aa_ref <- rep(NA_character_, n_rows)
  out$aa_alt <- rep(NA_character_, n_rows)
  out$codon_ref <- rep(NA_character_, n_rows)
  out$codon_alt <- rep(NA_character_, n_rows)
  out$codon_new <- rep(NA_character_, n_rows)
  out$consequence <- rep(NA_character_, n_rows)
  out$region <- rep(NA_character_, n_rows)
  out$strand <- rep(NA_character_, n_rows)
  out$qc_note <- rep(NA_character_, n_rows)

  # Process each row
  for (i in seq_len(nrow(out))) {
    row_type <- toupper(as.character(out$type[i]))

    if (row_type == "SNP") {
      enriched_row <- .pm_enrich_snp(gd, out[i, , drop = FALSE], flank, quiet)
      # Copy enriched fields back
      for (col in c("dna_ref", "dna_alt", "aa_ref", "aa_alt",
                    "codon_ref", "codon_alt", "codon_new", "consequence", "region", "strand", "qc_note")) {
        out[i, col] <- enriched_row[[col]]
      }
    }
    # Future: handle DEL, INS, SUB types here
  }

  if (!quiet) {
    n_coding <- sum(out$region == "coding", na.rm = TRUE)
    n_intergenic <- sum(out$region == "intergenic", na.rm = TRUE)
    cli::cli_inform(c(
      "v" = "Enriched {n_coding} coding and {n_intergenic} intergenic mutation{?s}"
    ))
  }

  out
}


#' Append QC Note to Row
#'
#' @description
#' Helper to accumulate QC notes in a row's qc_note column.
#'
#' @param row Single-row data.frame
#' @param note Character string note to append
#' @return Modified row with note appended
#' @keywords internal
.pm_append_qc_note <- function(row, note) {
  if (is.na(row$qc_note) || !nzchar(row$qc_note)) {
    row$qc_note <- note
  } else {
    row$qc_note <- paste(row$qc_note, note, sep = "; ")
  }
  row
}


#' Parse Annotation Geometry
#'
#' @description
#' Extracts region classification and coding geometry from annotation field.
#' Tries to parse patterns like "coding (45/1200 nt)" before falling back to gd_locate().
#'
#' @param annotation Character string from annotation column
#' @return List with region, coding_pos, coding_len (may be NA)
#' @keywords internal
.pm_parse_annotation_geometry <- function(annotation) {
  if (is.na(annotation) || !nzchar(annotation)) {
    return(list(region = NA_character_, coding_pos = NA_integer_, coding_len = NA_integer_))
  }

  # Pattern: "coding (45/1200 nt)" or "coding ( 45 / 1200 nt )"
  coding_match <- regexec("coding\\s*\\(\\s*(\\d+)\\s*/\\s*(\\d+)\\s*nt\\s*\\)", annotation)
  matches <- regmatches(annotation, coding_match)[[1]]

  if (length(matches) == 3) {
    return(list(
      region = "coding",
      coding_pos = as.integer(matches[2]),
      coding_len = as.integer(matches[3])
    ))
  }

  # Pattern: "intergenic (+123/-456)" or just "intergenic"
  if (grepl("intergenic", annotation, ignore.case = TRUE)) {
    return(list(region = "intergenic", coding_pos = NA_integer_, coding_len = NA_integer_))
  }

  # Unknown format
  list(region = NA_character_, coding_pos = NA_integer_, coding_len = NA_integer_)
}


#' Enrich Single SNP Row
#'
#' @description
#' Internal helper to enrich a single SNP mutation. Delegates to `gd_verify_snp()`
#' for coding regions or `get_roi_dna()` for intergenic regions.
#'
#' @param gd GenomeData object
#' @param row Single-row data.frame with mutation info
#' @param flank Integer, flanking window size for intergenic regions
#' @param quiet Logical, suppress warnings
#'
#' @return The input row with enriched columns filled in
#' @keywords internal
.pm_enrich_snp <- function(gd, row, flank, quiet = FALSE) {
  # Parse position (remove commas, handle P:k format)
  parsed <- .pm_parse_position(as.character(row$position), return_offset = TRUE)
  if (is.na(parsed$position)) {
    if (!quiet) {
      cli::cli_warn("Could not parse position: {row$position}")
    }
    row <- .pm_append_qc_note(row, "Position parse failed")
    return(row)
  }

  pos_parsed <- parsed$position

  # Warn if k != 0 for SNPs (indicates mutation within deletion/insertion context)
  if (parsed$offset != 0L) {
    qc_msg <- sprintf("SNP has non-zero offset (%d) - unusual", parsed$offset)
    if (!quiet) {
      cli::cli_warn(qc_msg)
    }
    row <- .pm_append_qc_note(row, qc_msg)
  }

  # Parse mutation string (A->C format)
  mut <- .pm_parse_snp_mutation(as.character(row$mutation))
  if (is.na(mut$ref) || is.na(mut$alt)) {
    row$consequence <- "unknown_format"
    row <- .pm_append_qc_note(row, "Unparseable mutation format")
    return(row)
  }

  # Classify region (try annotation field first, then gd_locate)
  seq_id_chr <- as.character(row$seq_id)

  # Try parsing annotation field for geometry
  ann_geo <- if ("annotation" %in% names(row) && !is.na(row$annotation)) {
    .pm_parse_annotation_geometry(as.character(row$annotation))
  } else {
    list(region = NA_character_, coding_pos = NA_integer_, coding_len = NA_integer_)
  }

  # Use annotation geometry if available, otherwise fall back to gd_locate()
  if (!is.na(ann_geo$region)) {
    loc_region <- ann_geo$region
  } else {
    loc <- tryCatch(
      gd_locate(gd, seq_id = seq_id_chr, pos = pos_parsed),
      error = function(e) {
        if (!quiet) cli::cli_warn("gd_locate failed for {seq_id_chr}:{pos_parsed}")
        list(region = NA_character_, genes = NA_character_)
      }
    )
    loc_region <- loc$region
    row <- .pm_append_qc_note(row, "Used gd_locate fallback")
  }

  if (is.na(loc_region) || loc_region == "intergenic") {
    # Intergenic: extract DNA window only
    start_pos <- max(1L, pos_parsed - flank)
    end_pos <- pos_parsed + flank

    dna_window <- tryCatch(
      get_roi_dna(gd, chrom = seq_id_chr, start = start_pos, end = end_pos, strand = "+"),
      error = function(e) NA_character_
    )

    if (!is.na(dna_window)) {
      row$dna_ref <- dna_window
      # Apply SNP at the center position (flank + 1 in 1-based)
      window_offset <- pos_parsed - start_pos + 1L
      row$dna_alt <- .pm_apply_snp_to_sequence(dna_window, window_offset, mut$alt)
    }

    row$region <- "intergenic"
    return(row)
  }

  # Coding region: use gd_verify_snp()
  gene <- as.character(row$gene)
  if (is.na(gene) || !nzchar(gene) || gene == "-") {
    # Fallback: search for overlapping CDS
    cds <- tryCatch(
      search_features(gd, type = "CDS", seqname = seq_id_chr,
                     start = pos_parsed, end = pos_parsed),
      error = function(e) data.frame()
    )

    if (nrow(cds) > 0) {
      gene <- if ("gene" %in% names(cds) && !is.na(cds$gene[1])) {
        as.character(cds$gene[1])
      } else if ("locus_tag" %in% names(cds) && !is.na(cds$locus_tag[1])) {
        as.character(cds$locus_tag[1])
      } else {
        NA_character_
      }
    }
  }

  if (!is.na(gene) && nzchar(gene)) {
    # Compute SNP effect using existing utility
    effect <- tryCatch(
      gd_verify_snp(gd, seq_id = seq_id_chr, position = pos_parsed, gene = gene),
      error = function(e) {
        if (!quiet) cli::cli_warn("gd_verify_snp failed for {gene}")
        data.frame(ok = FALSE)
      }
    )

    if (effect$ok) {
      # Get full reference CDS and AA sequences
      dna_ref_full <- tryCatch(
        get_gene_dna(gd, gene),
        error = function(e) NA_character_
      )

      aa_ref_full <- tryCatch(
        get_gene_aa(gd, gene),
        error = function(e) NA_character_
      )

      if (!is.na(dna_ref_full) && !is.na(aa_ref_full)) {
        row$dna_ref <- dna_ref_full

        # Create mutated CDS (apply SNP at cds_pos)
        cds_alt <- dna_ref_full
        substr(cds_alt, effect$cds_pos, effect$cds_pos) <-
          if (effect$strand == "-") {
            gd_complement_base(mut$alt)
          } else {
            mut$alt
          }
        row$dna_alt <- cds_alt

        # Translate mutated CDS
        aa_alt_raw <- tryCatch(
          translate_dna(cds_alt, frame = 1, genetic_code = "11", .internal = TRUE),
          error = function(e) NA_character_
        )

        # Normalize alternative start codons to match get_gene_aa() behavior
        # get_gene_aa() normalizes GTG, TTG, CTG, ATT, ATC, ATA → M at position 1
        aa_alt_full <- aa_alt_raw
        if (!is.na(aa_alt_raw) && nchar(aa_alt_raw) > 0 && nchar(cds_alt) >= 3) {
          first_codon <- toupper(substr(cds_alt, 1, 3))
          alt_start_codons <- c("GTG", "TTG", "CTG", "ATT", "ATC", "ATA")

          if (first_codon %in% alt_start_codons && substr(aa_alt_raw, 1, 1) != "M") {
            # Replace first AA with M
            aa_alt_full <- paste0("M", substr(aa_alt_raw, 2, nchar(aa_alt_raw)))
            qc_msg <- sprintf("Normalized alt start %s->M", first_codon)
            row <- .pm_append_qc_note(row, qc_msg)
          }
        }

        row$aa_ref <- aa_ref_full
        row$aa_alt <- aa_alt_full
        row$codon_ref <- effect$codon_ref
        row$codon_alt <- effect$codon_new
        row$codon_new <- effect$codon_new
        row$consequence <- effect$consequence
        row$region <- "coding"
        row$strand <- effect$strand
      }
    }
  }

  row
}


#' Parse Position String
#'
#' @description
#' Parses position strings from predict_mutations() output, handling:
#' - Comma-formatted numbers (e.g., "1,234,567")
#' - breseq P:k format where P is genomic position, k is offset
#'
#' @param position Character string position
#' @param return_offset Logical, if TRUE return list(position, offset) instead of just position
#' @return Integer position, or list(position, offset) if return_offset=TRUE. NA if parse fails.
#' @keywords internal
.pm_parse_position <- function(position, return_offset = FALSE) {
  if (is.na(position) || !nzchar(position)) {
    if (return_offset) {
      return(list(position = NA_integer_, offset = 0L))
    } else {
      return(NA_integer_)
    }
  }

  # Remove commas
  pos_clean <- gsub(",", "", position, fixed = TRUE)

  # Handle P:k format (breseq notation for within-deletion positions)
  # P is base position, k is offset indicator
  k_offset <- 0L
  if (grepl(":", pos_clean, fixed = TRUE)) {
    parts <- strsplit(pos_clean, ":", fixed = TRUE)[[1]]
    pos_clean <- parts[1]
    if (length(parts) > 1) {
      k_offset <- suppressWarnings(as.integer(parts[2]))
      if (is.na(k_offset)) k_offset <- 0L
    }
  }

  # Convert to integer
  pos_int <- suppressWarnings(as.integer(pos_clean))

  if (is.na(pos_int) || pos_int < 1) {
    if (return_offset) {
      return(list(position = NA_integer_, offset = k_offset))
    } else {
      return(NA_integer_)
    }
  }

  if (return_offset) {
    list(position = pos_int, offset = k_offset)
  } else {
    pos_int
  }
}


#' Parse SNP Mutation String
#'
#' @description
#' Parses SNP mutation strings like "A->C" or "A>C" to extract
#' reference and alternate bases. Supports HTML-escaped formats.
#'
#' @param mutation Character string mutation
#' @return List with `ref` and `alt` bases, or NAs if parse fails
#' @keywords internal
.pm_parse_snp_mutation <- function(mutation) {
  if (is.na(mutation) || !nzchar(mutation)) {
    return(list(ref = NA_character_, alt = NA_character_))
  }

  # Try different arrow formats (including HTML-escaped)
  sep_patterns <- c(
    "\u2192",    # → (unicode arrow)
    "->",        # ASCII arrow
    ">",         # Simple greater-than
    "-&gt;",     # HTML-escaped arrow
    "&gt;"       # HTML-escaped greater-than
  )

  for (sep in sep_patterns) {
    if (grepl(sep, mutation, fixed = TRUE)) {
      parts <- strsplit(mutation, sep, fixed = TRUE)[[1]]
      if (length(parts) == 2) {
        ref <- trimws(parts[1])
        alt <- trimws(parts[2])

        # Validate bases
        if (nchar(ref) == 1 && nchar(alt) == 1 &&
            ref %in% c("A", "C", "G", "T") &&
            alt %in% c("A", "C", "G", "T")) {
          return(list(ref = ref, alt = alt))
        }
      }
    }
  }

  # Parse failed
  list(ref = NA_character_, alt = NA_character_)
}


#' Apply SNP to Sequence
#'
#' @description
#' Simple helper to substitute a single base in a DNA string.
#'
#' @param seq Character string DNA sequence
#' @param pos Integer, 1-based position to substitute
#' @param alt_base Character, new base (A/C/G/T)
#' @return Modified sequence, or NA if position is out of bounds
#' @keywords internal
.pm_apply_snp_to_sequence <- function(seq, pos, alt_base) {
  if (is.na(seq) || is.na(pos) || is.na(alt_base)) return(NA_character_)
  if (pos < 1 || pos > nchar(seq)) return(NA_character_)

  paste0(
    substr(seq, 1, pos - 1),
    alt_base,
    substr(seq, pos + 1, nchar(seq))
  )
}

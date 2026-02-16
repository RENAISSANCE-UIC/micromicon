#' Enrich Mutation Table with Consequences
#'
#' @description
#' Adds molecular consequence annotations to a `predict_mutations()` table.
#' Supports SNPs, deletions (DEL), insertions (INS), and substitutions (SUB).
#'
#' For coding region mutations, computes reference/alternate sequences and
#' consequences. For intergenic mutations, extracts flanking DNA windows.
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
#' **Supported mutation types:**
#' \itemize{
#'   \item \strong{SNP}: synonymous, missense, nonsense
#'   \item \strong{DEL}: inframe_deletion, frameshift
#'   \item \strong{INS}: inframe_insertion, frameshift
#'   \item \strong{SUB}: complex (partial support)
#' }
#'
#' **Note**: Only coding region mutations get amino acid sequences. Intergenic
#' mutations only receive DNA sequences (no protein translation).
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
    } else if (row_type == "DEL") {
      enriched_row <- .pm_enrich_del(gd, out[i, , drop = FALSE], flank, quiet)
      for (col in c("dna_ref", "dna_alt", "aa_ref", "aa_alt",
                    "codon_ref", "codon_alt", "codon_new", "consequence", "region", "strand", "qc_note")) {
        out[i, col] <- enriched_row[[col]]
      }
    } else if (row_type == "INS") {
      enriched_row <- .pm_enrich_ins(gd, out[i, , drop = FALSE], flank, quiet)
      for (col in c("dna_ref", "dna_alt", "aa_ref", "aa_alt",
                    "codon_ref", "codon_alt", "codon_new", "consequence", "region", "strand", "qc_note")) {
        out[i, col] <- enriched_row[[col]]
      }
    } else if (row_type == "SUB") {
      enriched_row <- .pm_enrich_sub(gd, out[i, , drop = FALSE], flank, quiet)
      for (col in c("dna_ref", "dna_alt", "aa_ref", "aa_alt",
                    "codon_ref", "codon_alt", "codon_new", "consequence", "region", "strand", "qc_note")) {
        out[i, col] <- enriched_row[[col]]
      }
    }
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

  # Pattern: "coding (45/1200 nt)" or "coding (45/?? nt)" or "coding (45/ nt)"
  # The length part can be: digits, ??, or missing entirely
  # Use a more flexible pattern that captures any non-space characters or nothing
  coding_match <- regexec("coding\\s*\\(\\s*(\\d+)\\s*/\\s*([^\\s)]*?)\\s*nt\\s*\\)", annotation)
  matches <- regmatches(annotation, coding_match)[[1]]

  if (length(matches) >= 2) {
    # We have at least "coding" and position
    coding_pos <- as.integer(matches[2])
    # Try to parse length (may be NA if it's "??" or empty)
    coding_len <- if (length(matches) >= 3 && nzchar(matches[3]) && !grepl("\\?", matches[3])) {
      as.integer(matches[3])
    } else {
      NA_integer_
    }

    return(list(
      region = "coding",
      coding_pos = coding_pos,
      coding_len = coding_len
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
  gene <- .pm_resolve_gene(gd, as.character(row$gene), row$seq_id, row$position)

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


#' Clean and Resolve Gene Name
#'
#' @description
#' Cleans gene name from predict_mutations() output (removes arrows, whitespace)
#' and resolves to actual gene identifier. Falls back to position-based lookup.
#'
#' @param gd GenomeData object
#' @param gene_raw Raw gene string from mutation table (may have arrows)
#' @param seq_id Sequence ID
#' @param position Position (can be comma-formatted or P:k format)
#'
#' @return Cleaned gene name, or NA if not resolvable
#' @keywords internal
.pm_resolve_gene <- function(gd, gene_raw, seq_id, position) {
  # Clean gene name: remove arrows (←, →, <-, ->), trim whitespace
  gene_clean <- gene_raw
  gene_clean <- gsub("←|→|<-|->", "", gene_clean)
  gene_clean <- trimws(gene_clean)

  # If gene is NA or empty after cleaning, try position-based lookup
  if (is.na(gene_clean) || !nzchar(gene_clean)) {
    # Parse position
    pos_parsed <- .pm_parse_position(position)
    if (!is.na(pos_parsed)) {
      # Use gd_locate to find gene at this position
      loc <- tryCatch(
        gd_locate(gd, seq_id = as.character(seq_id), pos = pos_parsed),
        error = function(e) list(genes = NA_character_)
      )
      if (!is.na(loc$genes) && length(loc$genes) > 0) {
        gene_clean <- loc$genes[1]  # Take first gene if multiple
      }
    }
  }

  # Verify gene exists - try to find it
  if (!is.na(gene_clean) && nzchar(gene_clean)) {
    # Try finding by pattern (searches Name, ID, Alias, gene, locus_tag, product)
    features <- tryCatch(
      search_features(gd, pattern = gene_clean, type = "CDS"),
      error = function(e) data.frame()
    )

    if (nrow(features) > 0) {
      # Use the ID field for most reliable lookup
      if ("ID" %in% names(features) && !is.na(features$ID[1])) {
        return(as.character(features$ID[1]))
      }
      # Fall back to Name
      if ("Name" %in% names(features) && !is.na(features$Name[1])) {
        return(as.character(features$Name[1]))
      }
      # Fall back to locus_tag
      if ("locus_tag" %in% names(features) && !is.na(features$locus_tag[1])) {
        return(as.character(features$locus_tag[1]))
      }
      # Last resort: use cleaned name
      return(gene_clean)
    }
  }

  # Couldn't resolve
  NA_character_
}


#' Enrich Single DEL Row
#'
#' @description
#' Internal helper to enrich a single deletion mutation.
#'
#' @param gd GenomeData object
#' @param row Single-row data.frame with mutation info
#' @param flank Integer, flanking window size for intergenic regions
#' @param quiet Logical, suppress warnings
#'
#' @return The input row with enriched columns filled in
#' @keywords internal
.pm_enrich_del <- function(gd, row, flank, quiet = FALSE) {
  # Parse annotation to determine region
  ann_geo <- if ("annotation" %in% names(row) && !is.na(row$annotation)) {
    .pm_parse_annotation_geometry(as.character(row$annotation))
  } else {
    list(region = NA_character_, coding_pos = NA_integer_, coding_len = NA_integer_)
  }

  # Parse deletion annotation for range
  del_info <- .pm_parse_deletion_annotation(as.character(row$annotation))

  if (!is.na(ann_geo$region) && ann_geo$region == "coding" && !is.na(del_info$start) && !is.na(del_info$end)) {
    # Coding deletion
    # Resolve gene name (clean arrows, fall back to position lookup)
    gene <- .pm_resolve_gene(gd, as.character(row$gene), row$seq_id, row$position)

    if (is.na(gene)) {
      row <- .pm_append_qc_note(row, "Could not resolve gene name")
      return(row)
    }

    # Get reference sequences
    dna_ref_full <- tryCatch(get_gene_dna(gd, gene), error = function(e) NA_character_)
    aa_ref_full <- tryCatch(get_gene_aa(gd, gene), error = function(e) NA_character_)

    if (!is.na(dna_ref_full) && !is.na(aa_ref_full)) {
      row$dna_ref <- dna_ref_full
      row$aa_ref <- aa_ref_full

      # Apply deletion (remove bases from start to end)
      if (del_info$start > 1) {
        before_del <- substr(dna_ref_full, 1, del_info$start - 1)
      } else {
        before_del <- ""
      }

      if (del_info$end < nchar(dna_ref_full)) {
        after_del <- substr(dna_ref_full, del_info$end + 1, nchar(dna_ref_full))
      } else {
        after_del <- ""
      }

      dna_alt <- paste0(before_del, after_del)
      row$dna_alt <- dna_alt

      # Translate mutated sequence
      aa_alt <- tryCatch(
        translate_dna(dna_alt, frame = 1, genetic_code = "11", .internal = TRUE),
        error = function(e) NA_character_
      )

      # Normalize alt-start if needed
      if (!is.na(aa_alt) && nchar(aa_alt) > 0 && nchar(dna_alt) >= 3) {
        first_codon <- toupper(substr(dna_alt, 1, 3))
        alt_start_codons <- c("GTG", "TTG", "CTG", "ATT", "ATC", "ATA")
        if (first_codon %in% alt_start_codons && substr(aa_alt, 1, 1) != "M") {
          aa_alt <- paste0("M", substr(aa_alt, 2, nchar(aa_alt)))
          row <- .pm_append_qc_note(row, sprintf("Normalized alt start %s->M", first_codon))
        }
      }

      row$aa_alt <- aa_alt

      # Determine consequence
      del_size <- del_info$end - del_info$start + 1
      if (del_size %% 3 == 0) {
        row$consequence <- "inframe_deletion"
      } else {
        row$consequence <- "frameshift"
      }

      row$region <- "coding"
      # Get strand info
      strand_info <- tryCatch(gd_locate(gd, seq_id = as.character(row$seq_id), pos = del_info$start),
                             error = function(e) list(strand = NA_character_))
      row$strand <- strand_info$strand %||% NA_character_
    }
  } else {
    # Intergenic deletion - just provide DNA context, no AA
    pos_parsed <- .pm_parse_position(as.character(row$position))
    if (!is.na(pos_parsed)) {
      seq_id_chr <- as.character(row$seq_id)
      start_pos <- max(1L, pos_parsed - flank)
      end_pos <- pos_parsed + flank

      dna_window <- tryCatch(
        get_roi_dna(gd, chrom = seq_id_chr, start = start_pos, end = end_pos, strand = "+"),
        error = function(e) NA_character_
      )

      if (!is.na(dna_window)) {
        row$dna_ref <- dna_window
        # For intergenic, we can't easily apply the deletion without more info
        row$dna_alt <- NA_character_
        row <- .pm_append_qc_note(row, "Intergenic deletion - dna_alt not computed")
      }

      row$region <- "intergenic"
    }
  }

  row
}


#' Enrich Single INS Row
#'
#' @description
#' Internal helper to enrich a single insertion mutation.
#'
#' @param gd GenomeData object
#' @param row Single-row data.frame with mutation info
#' @param flank Integer, flanking window size for intergenic regions
#' @param quiet Logical, suppress warnings
#'
#' @return The input row with enriched columns filled in
#' @keywords internal
.pm_enrich_ins <- function(gd, row, flank, quiet = FALSE) {
  # Parse annotation
  ann_geo <- if ("annotation" %in% names(row) && !is.na(row$annotation)) {
    .pm_parse_annotation_geometry(as.character(row$annotation))
  } else {
    list(region = NA_character_, coding_pos = NA_integer_, coding_len = NA_integer_)
  }

  # Parse insertion
  ins_info <- .pm_parse_insertion_mutation(as.character(row$mutation))

  if (!is.na(ann_geo$region) && ann_geo$region == "coding" && !is.na(ins_info$position) && !is.na(ins_info$sequence)) {
    # Coding insertion
    gene <- .pm_resolve_gene(gd, as.character(row$gene), row$seq_id, row$position)

    if (is.na(gene)) {
      row <- .pm_append_qc_note(row, "Could not resolve gene name for insertion")
      return(row)
    }

    # Get reference sequences
    dna_ref_full <- tryCatch(get_gene_dna(gd, gene), error = function(e) NA_character_)
    aa_ref_full <- tryCatch(get_gene_aa(gd, gene), error = function(e) NA_character_)

    if (!is.na(dna_ref_full) && !is.na(aa_ref_full)) {
      row$dna_ref <- dna_ref_full
      row$aa_ref <- aa_ref_full

      # Apply insertion
      before_ins <- substr(dna_ref_full, 1, ins_info$position)
      after_ins <- substr(dna_ref_full, ins_info$position + 1, nchar(dna_ref_full))
      dna_alt <- paste0(before_ins, ins_info$sequence, after_ins)
      row$dna_alt <- dna_alt

      # Translate
      aa_alt <- tryCatch(
        translate_dna(dna_alt, frame = 1, genetic_code = "11", .internal = TRUE),
        error = function(e) NA_character_
      )

      # Normalize alt-start if needed
      if (!is.na(aa_alt) && nchar(aa_alt) > 0 && nchar(dna_alt) >= 3) {
        first_codon <- toupper(substr(dna_alt, 1, 3))
        alt_start_codons <- c("GTG", "TTG", "CTG", "ATT", "ATC", "ATA")
        if (first_codon %in% alt_start_codons && substr(aa_alt, 1, 1) != "M") {
          aa_alt <- paste0("M", substr(aa_alt, 2, nchar(aa_alt)))
          row <- .pm_append_qc_note(row, sprintf("Normalized alt start %s->M", first_codon))
        }
      }

      row$aa_alt <- aa_alt

      # Determine consequence
      if (nchar(ins_info$sequence) %% 3 == 0) {
        row$consequence <- "inframe_insertion"
      } else {
        row$consequence <- "frameshift"
      }

      row$region <- "coding"
      strand_info <- tryCatch(gd_locate(gd, seq_id = as.character(row$seq_id), pos = ins_info$position),
                             error = function(e) list(strand = NA_character_))
      row$strand <- strand_info$strand %||% NA_character_
    }
  } else {
    # Intergenic insertion
    pos_parsed <- .pm_parse_position(as.character(row$position))
    if (!is.na(pos_parsed)) {
      seq_id_chr <- as.character(row$seq_id)
      start_pos <- max(1L, pos_parsed - flank)
      end_pos <- pos_parsed + flank

      dna_window <- tryCatch(
        get_roi_dna(gd, chrom = seq_id_chr, start = start_pos, end = end_pos, strand = "+"),
        error = function(e) NA_character_
      )

      if (!is.na(dna_window)) {
        row$dna_ref <- dna_window
        row$dna_alt <- NA_character_
        row <- .pm_append_qc_note(row, "Intergenic insertion - dna_alt not computed")
      }

      row$region <- "intergenic"
    }
  }

  row
}


#' Enrich Single SUB Row
#'
#' @description
#' Internal helper to enrich a single substitution mutation.
#'
#' @param gd GenomeData object
#' @param row Single-row data.frame with mutation info
#' @param flank Integer, flanking window size for intergenic regions
#' @param quiet Logical, suppress warnings
#'
#' @return The input row with enriched columns filled in
#' @keywords internal
.pm_enrich_sub <- function(gd, row, flank, quiet = FALSE) {
  # Substitutions are complex - for now treat similarly to deletions
  # Future: implement full substitution logic
  row <- .pm_append_qc_note(row, "SUB mutations not fully implemented")

  # Parse annotation
  ann_geo <- if ("annotation" %in% names(row) && !is.na(row$annotation)) {
    .pm_parse_annotation_geometry(as.character(row$annotation))
  } else {
    list(region = NA_character_, coding_pos = NA_integer_, coding_len = NA_integer_)
  }

  if (!is.na(ann_geo$region)) {
    row$region <- ann_geo$region
  } else {
    row$region <- "intergenic"
  }

  row
}


#' Parse Deletion Annotation
#'
#' @description
#' Extracts deletion range from annotation like "coding (152-278 / 435 nt)".
#'
#' @param annotation Character string from annotation column
#' @return List with start, end, total_len
#' @keywords internal
.pm_parse_deletion_annotation <- function(annotation) {
  if (is.null(annotation) || length(annotation) == 0 || is.na(annotation) || !nzchar(annotation)) {
    return(list(start = NA_integer_, end = NA_integer_, total_len = NA_integer_))
  }

  # Pattern: "coding (152-278 / 435 nt)"
  del_match <- regexec("coding\\s*\\(\\s*(\\d+)\\s*-\\s*(\\d+)\\s*/\\s*([^\\s)]*?)\\s*nt\\s*\\)", annotation)
  matches <- regmatches(annotation, del_match)[[1]]

  if (length(matches) >= 3) {
    start <- as.integer(matches[2])
    end <- as.integer(matches[3])
    total_len <- if (length(matches) >= 4 && !grepl("\\?", matches[4])) {
      as.integer(matches[4])
    } else {
      NA_integer_
    }

    return(list(start = start, end = end, total_len = total_len))
  }

  # Couldn't parse
  list(start = NA_integer_, end = NA_integer_, total_len = NA_integer_)
}


#' Parse Insertion Mutation
#'
#' @description
#' Extracts insertion details from mutation field like "+ACGT" or "ins_ACGT".
#'
#' @param mutation Character string from mutation column
#' @return List with position (from annotation) and sequence
#' @keywords internal
.pm_parse_insertion_mutation <- function(mutation) {
  if (is.null(mutation) || length(mutation) == 0 || is.na(mutation) || !nzchar(mutation)) {
    return(list(position = NA_integer_, sequence = NA_character_))
  }

  # Pattern: "+ACGT" or "ins_ACGT"
  if (grepl("^\\+", mutation)) {
    seq <- sub("^\\+", "", mutation)
    return(list(position = NA_integer_, sequence = seq))
  }

  if (grepl("^ins_", mutation, ignore.case = TRUE)) {
    seq <- sub("^ins_", "", mutation, ignore.case = TRUE)
    return(list(position = NA_integer_, sequence = seq))
  }

  # Couldn't parse
  list(position = NA_integer_, sequence = NA_character_)
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

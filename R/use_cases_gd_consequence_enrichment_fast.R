#' Fast Mutation Consequence Enrichment
#'
#' @description
#' Optimized version of pm_enrich_consequences() with:
#' - Codon-level computation (don't translate full CDS)
#' - Gene-based batching (process all mutations per gene together)
#' - Minimal sequence operations
#'
#' @param gd A GenomeData object
#' @param pm_tbl A data.frame from `predict_variants()`
#' @param flank Integer, flanking bases for intergenic regions (default: 50)
#' @param quiet Logical, suppress messages (default: FALSE)
#'
#' @return Enriched data.frame with consequence columns
#' @keywords internal
pm_enrich_consequences_fast <- function(gd, pm_tbl, flank = 50L, quiet = FALSE) {
  gd_assert(gd, "gd")
  stopifnot(is.data.frame(pm_tbl))

  if (!quiet) {
    cli::cli_inform(c(
      "i" = "pm_enrich_consequences_fast(): enriching {nrow(pm_tbl)} mutation{?s}"
    ))
    start_time <- Sys.time()
  }

  # Initialize output columns
  out <- pm_tbl
  n_rows <- nrow(out)
  out$dna_ref <- rep(NA_character_, n_rows)
  out$dna_alt <- rep(NA_character_, n_rows)
  out$aa_ref <- rep(NA_character_, n_rows)
  out$aa_alt <- rep(NA_character_, n_rows)
  out$codon_ref <- rep(NA_character_, n_rows)
  out$codon_alt <- rep(NA_character_, n_rows)
  out$consequence <- if ("var_type" %in% names(pm_tbl)) pm_tbl$var_type else rep(NA_character_, n_rows)
  out$region <- rep(NA_character_, n_rows)
  out$strand <- rep(NA_character_, n_rows)
  out$qc_note <- rep(NA_character_, n_rows)

  # Group mutations by gene for batch processing
  out$..row_id <- seq_len(n_rows)

  # Separate intergenic mutations from coding ones (by annotation)
  out$..is_intergenic <- grepl("intergenic", out$annotation, ignore.case = TRUE)

  # Split into coding vs intergenic
  coding_rows <- out[!out$..is_intergenic, , drop = FALSE]
  intergenic_rows <- out[out$..is_intergenic, , drop = FALSE]

  # Split coding mutations by gene
  gene_groups <- if (nrow(coding_rows) > 0) {
    split(coding_rows, coding_rows$gene, drop = FALSE)
  } else {
    list()
  }

  # Global caches (shared across all genes)
  gene_dna_cache <- new.env(parent = emptyenv())
  gene_aa_cache <- new.env(parent = emptyenv())
  gene_meta_cache <- new.env(parent = emptyenv())  # Store strand, feature info

  # Process genes with progress bar
  gene_names <- names(gene_groups)

  if (!quiet && length(gene_names) > 0) {
    pb <- txtProgressBar(min = 0, max = length(gene_names), style = 3)
  }

  for (i in seq_along(gene_names)) {
    gene <- gene_names[i]
    group <- gene_groups[[gene]]

    if (is.na(gene) || !nzchar(gene) || nrow(group) == 0) {
      if (!quiet && length(gene_names) > 0) {
        setTxtProgressBar(pb, i)
      }
      next
    }

    # Process all mutations for this gene together
    enriched_group <- .pf_enrich_gene_batch(gd, group, flank, quiet,
                                            gene_dna_cache, gene_aa_cache, gene_meta_cache)

    # Write back to output
    for (col in c("dna_ref", "dna_alt", "aa_ref", "aa_alt",
                  "codon_ref", "codon_alt", "consequence",
                  "region", "strand", "qc_note")) {
      out[enriched_group$..row_id, col] <- enriched_group[[col]]
    }

    if (!quiet && length(gene_names) > 0) {
      setTxtProgressBar(pb, i)
    }
  }

  if (!quiet && length(gene_names) > 0) {
    close(pb)
    cat("\n")  # New line after progress bar
  }

  # Handle mutations that weren't enriched (intergenic or failed gene enrichment)
  # Check for region = NA (not enriched yet) or explicitly intergenic
  unenriched_idx <- which(is.na(out$region) | out$region == "intergenic")
  if (length(unenriched_idx) > 0) {
    for (i in unenriched_idx) {
      out[i, ] <- .pf_enrich_intergenic(gd, out[i, , drop = FALSE], flank, quiet)
    }
  }

  out$..row_id <- NULL
  out$..is_intergenic <- NULL

  if (!quiet) {
    elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
    n_coding <- sum(out$region == "coding", na.rm = TRUE)
    n_intergenic <- sum(out$region == "intergenic", na.rm = TRUE)
    cli::cli_inform(c(
      "v" = "Enriched {n_coding} coding and {n_intergenic} intergenic mutation{?s} in {elapsed}s"
    ))
  }

  out
}


#' Batch Enrich All Mutations for a Single Gene
#'
#' @description
#' Process all mutations in a gene together since they share:
#' - Same ref DNA sequence
#' - Same ref AA sequence
#' - Same strand
#' - Same genomic coordinates
#'
#' @keywords internal
.pf_enrich_gene_batch <- function(gd, group, flank, quiet,
                                   gene_dna_cache, gene_aa_cache, gene_meta_cache) {
  gene <- as.character(group$gene[1])
  gene_clean <- gsub("\\s*[\u2192\u2190]\\s*$", "", gene)
  gene_clean <- trimws(gene_clean)

  # Handle multi-gene annotations (e.g., "geneA|geneB" or "geneA / geneB")
  # Take the first gene in the list
  gene_clean <- gsub("\\|.*$", "", gene_clean)  # Remove everything after |
  gene_clean <- gsub("\\s*/\\s*.*$", "", gene_clean)  # Remove everything after /
  gene_clean <- trimws(gene_clean)

  # Get gene metadata once (cached)
  gene_meta <- .pf_get_gene_metadata(gd, gene_clean, gene_meta_cache)

  if (is.null(gene_meta) || is.na(gene_meta$strand)) {
    # Can't process this gene
    return(group)
  }

  # Get ref sequences once (cached)
  dna_ref <- .pm_cached_get_gene_dna(gd, gene_clean, gene_dna_cache)
  aa_ref <- .pm_cached_translate_gene(dna_ref, gene_clean, gene_aa_cache)

  if (is.na(dna_ref) || is.na(aa_ref)) {
    return(group)
  }

  # Fill in shared info for all mutations in this gene
  group$dna_ref <- dna_ref
  group$aa_ref <- aa_ref
  group$region <- "coding"
  group$strand <- gene_meta$strand

  # Process each mutation (now much faster - just codon-level work)
  for (i in seq_len(nrow(group))) {
    row_type <- toupper(as.character(group$type[i]))

    if (row_type == "SNP") {
      group[i, ] <- .pf_enrich_snp_fast(gd, group[i, , drop = FALSE],
                                        dna_ref, aa_ref, gene_meta, quiet)
    } else if (row_type == "DEL") {
      group[i, ] <- .pf_enrich_del_fast(gd, group[i, , drop = FALSE],
                                        dna_ref, aa_ref, gene_meta, quiet)
    } else if (row_type == "INS") {
      group[i, ] <- .pf_enrich_ins_fast(gd, group[i, , drop = FALSE],
                                        dna_ref, aa_ref, gene_meta, quiet)
    } else if (row_type %in% c("MOB", "AMP", "CON", "INV")) {
      group[i, ] <- .pf_enrich_structural_fast(group[i, , drop = FALSE])
    }
  }

  group
}


#' Get Gene Metadata (strand, feature info) - Cached
#'
#' @keywords internal
.pf_get_gene_metadata <- function(gd, gene, cache) {
  if (exists(gene, envir = cache, inherits = FALSE)) {
    return(cache[[gene]])
  }

  meta <- tryCatch({
    feat <- gd_resolve_feature(gd, gene)
    if (is.null(feat) || nrow(feat) == 0) {
      return(NULL)
    }

    # Only process CDS features (skip tRNA, rRNA, ncRNA, etc.)
    if (!is.na(feat$type[1]) && feat$type[1] != "CDS") {
      return(NULL)
    }

    list(
      strand = as.character(feat$strand[1]),
      seqname = as.character(feat$seqname[1]),
      start = feat$start[1],
      end = feat$end[1]
    )
  }, error = function(e) NULL)

  cache[[gene]] <- meta
  meta
}


#' Fast SNP Enrichment (Codon-Level Only)
#'
#' @description
#' HUGE optimization: Don't translate entire mutated CDS.
#' Just extract and translate the affected codon.
#'
#' @keywords internal
.pf_enrich_snp_fast <- function(gd, row, dna_ref, aa_ref, gene_meta, quiet) {
  # Parse position and mutation
  pos_parsed <- as.integer(gsub(",|:.*", "", row$position))

  mutation_str <- as.character(row$mutation)
  # Parse ref->alt (handle various arrow formats)
  mut_parts <- strsplit(gsub("\u2192|->|&gt;", "|", mutation_str), "\\|")[[1]]
  if (length(mut_parts) != 2) {
    return(row)
  }

  ref_base <- mut_parts[1]
  alt_base <- mut_parts[2]

  # Map genomic position to CDS position
  cds_pos_result <- tryCatch(
    map_genomic_to_cds(gd, as.character(row$gene), pos_parsed),
    error = function(e) NA_integer_
  )
  cds_pos <- if (is.list(cds_pos_result)) {
    scalar_or(cds_pos_result$cds_pos, NA_integer_)
  } else {
    cds_pos_result
  }

  if (is.na(cds_pos) || cds_pos < 1 || cds_pos > nchar(dna_ref)) {
    return(row)
  }

  # Compute codon boundaries
  aa_index <- ceiling(cds_pos / 3)
  codon_start <- (aa_index - 1) * 3 + 1
  codon_end <- min(codon_start + 2, nchar(dna_ref))

  # Extract reference codon
  codon_ref <- substr(dna_ref, codon_start, codon_end)
  row$codon_ref <- codon_ref

  # Apply SNP to codon (not entire CDS!)
  codon_offset <- cds_pos - codon_start + 1  # Position within codon (1-3)
  alt_base_cds <- if (gene_meta$strand == "-") {
    gd_complement_base(alt_base)
  } else {
    alt_base
  }

  codon_alt <- codon_ref
  substr(codon_alt, codon_offset, codon_offset) <- alt_base_cds
  row$codon_alt <- codon_alt

  # Translate just the two codons (not entire CDS!)
  aa_ref_char <- substr(aa_ref, aa_index, aa_index)
  aa_alt_char <- tryCatch({
    translate_dna(codon_alt, frame = 1, genetic_code = "11", .internal = TRUE)
  }, error = function(e) NA_character_)

  # For dna_alt, we still need full mutated CDS (but we can be lazy here)
  dna_alt <- dna_ref
  substr(dna_alt, cds_pos, cds_pos) <- alt_base_cds
  row$dna_alt <- dna_alt

  # Compute aa_alt (full protein) - but shortcut: copy ref and change one AA
  aa_alt <- aa_ref
  if (!is.na(aa_alt_char) && nchar(aa_alt_char) > 0) {
    substr(aa_alt, aa_index, aa_index) <- substr(aa_alt_char, 1, 1)
  }
  row$aa_alt <- aa_alt


  row
}


#' Fast DEL Enrichment
#'
#' @keywords internal
.pf_enrich_del_fast <- function(gd, row, dna_ref, aa_ref, gene_meta, quiet) {
  # Parse deletion info from annotation
  annotation <- as.character(row$annotation)
  del_info <- .pm_parse_deletion_annotation(annotation)

  if (is.na(del_info$start) || is.na(del_info$end)) {
    return(row)
  }

  # Apply deletion
  del_size <- del_info$end - del_info$start + 1

  if (del_info$start > 1) {
    before_del <- substr(dna_ref, 1, del_info$start - 1)
  } else {
    before_del <- ""
  }

  if (del_info$end < nchar(dna_ref)) {
    after_del <- substr(dna_ref, del_info$end + 1, nchar(dna_ref))
  } else {
    after_del <- ""
  }

  dna_alt <- paste0(before_del, after_del)
  row$dna_alt <- dna_alt


  # Translate altered sequence (still needed for aa_alt)
  aa_alt <- tryCatch(
    translate_dna(dna_alt, frame = 1, genetic_code = "11", .internal = TRUE),
    error = function(e) NA_character_
  )

  # Normalize alt start codon
  if (!is.na(aa_alt) && nchar(aa_alt) > 0 && nchar(dna_alt) >= 3) {
    first_codon <- toupper(substr(dna_alt, 1, 3))
    if (first_codon %in% c("GTG", "TTG", "CTG", "ATT", "ATC", "ATA") &&
        substr(aa_alt, 1, 1) != "M") {
      aa_alt <- paste0("M", substr(aa_alt, 2, nchar(aa_alt)))
    }
  }

  # Truncate at first stop codon (match original behavior)
  aa_alt <- .pm_truncate_at_stop(aa_alt)

  row$aa_alt <- aa_alt
  row
}


#' Fast INS Enrichment
#'
#' @keywords internal
.pf_enrich_ins_fast <- function(gd, row, dna_ref, aa_ref, gene_meta, quiet) {
  # Parse insertion info
  ins_info <- .pm_parse_insertion_mutation(as.character(row$mutation))

  if (is.na(ins_info$position) || is.na(ins_info$sequence)) {
    return(row)
  }

  # Apply insertion
  before_ins <- substr(dna_ref, 1, ins_info$position)
  after_ins <- substr(dna_ref, ins_info$position + 1, nchar(dna_ref))
  dna_alt <- paste0(before_ins, ins_info$sequence, after_ins)
  row$dna_alt <- dna_alt

  ins_size <- nchar(ins_info$sequence)

  # Translate altered sequence
  aa_alt <- tryCatch(
    translate_dna(dna_alt, frame = 1, genetic_code = "11", .internal = TRUE),
    error = function(e) NA_character_
  )

  # Normalize alt start
  if (!is.na(aa_alt) && nchar(aa_alt) > 0 && nchar(dna_alt) >= 3) {
    first_codon <- toupper(substr(dna_alt, 1, 3))
    if (first_codon %in% c("GTG", "TTG", "CTG", "ATT", "ATC", "ATA") &&
        substr(aa_alt, 1, 1) != "M") {
      aa_alt <- paste0("M", substr(aa_alt, 2, nchar(aa_alt)))
    }
  }

  # Truncate at first stop codon (match original behavior)
  aa_alt <- .pm_truncate_at_stop(aa_alt)

  row$aa_alt <- aa_alt
  row
}


#' Fast Structural Variant Enrichment
#'
#' @keywords internal
.pf_enrich_structural_fast <- function(row) {
  # Already have dna_ref, aa_ref from batch processing
  # Just set alt fields to NA
  row$dna_alt <- NA_character_
  row$aa_alt <- NA_character_
  row$codon_ref <- NA_character_
  row$codon_alt <- NA_character_

  variant_type <- toupper(as.character(row$type))
  row$qc_note <- sprintf("Structural variant (%s) - alternate sequences not computed", variant_type)

  row
}


#' Fast Intergenic Enrichment
#'
#' @keywords internal
.pf_enrich_intergenic <- function(gd, row, flank, quiet, dna_window = NULL) {
  # Extract DNA window (use proper position parser to handle colons)
  pos <- .pm_parse_position(as.character(row$position))
  seq_id <- as.character(row$seq_id)

  start_pos <- max(1L, pos - flank)
  end_pos <- pos + flank

  if (is.null(dna_window)) {
    dna_window <- tryCatch(
      get_roi_dna(gd, contig = seq_id, start = start_pos, end = end_pos, strand = "+"),
      error = function(e) NA_character_
    )
  }

  row$dna_ref <- dna_window
  row$region <- "intergenic"

  # Apply mutation to window based on type
  if (!is.na(dna_window) && !is.na(pos)) {
    window_offset <- pos - start_pos + 1L
    mut_type <- toupper(as.character(row$type))

    if (mut_type == "SNP") {
      # SNP: substitute single base
      mutation_str <- as.character(row$mutation)
      mut_parts <- strsplit(gsub("\u2192|->|&gt;", "|", mutation_str), "\\|")[[1]]
      if (length(mut_parts) == 2) {
        alt_base <- mut_parts[2]
        dna_alt <- dna_window
        substr(dna_alt, window_offset, window_offset) <- alt_base
        row$dna_alt <- dna_alt
      }

    } else if (mut_type == "DEL") {
      # DEL: remove bases from window
      mutation_str <- as.character(row$mutation)
      # Parse deletion size (e.g., "delta1 bp" or "delta5 bp")
      del_size <- 1L  # Default to 1
      if (grepl("\u0394(\\d+)", mutation_str)) {
        del_match <- regexpr("\u0394(\\d+)", mutation_str, perl = TRUE)
        del_text <- regmatches(mutation_str, del_match)
        if (length(del_text) > 0) {
          del_size <- as.integer(gsub("\u0394", "", del_text))
        }
      }

      # Apply deletion to window (with NA checks)
      if (!is.na(window_offset) && !is.na(del_size) &&
          window_offset >= 1 && window_offset + del_size - 1 <= nchar(dna_window)) {
        before_del <- if (window_offset > 1) substr(dna_window, 1, window_offset - 1) else ""
        after_del <- if (window_offset + del_size <= nchar(dna_window)) {
          substr(dna_window, window_offset + del_size, nchar(dna_window))
        } else {
          ""
        }
        row$dna_alt <- paste0(before_del, after_del)
      }

    } else if (mut_type == "INS") {
      # INS: insert bases into window
      mutation_str <- as.character(row$mutation)
      # Parse insertion sequence (e.g., "+ACGT" or "ins_ACGT")
      ins_seq <- NA_character_
      if (grepl("^\\+([ACGT]+)$", mutation_str)) {
        ins_seq <- gsub("^\\+", "", mutation_str)
      } else if (grepl("ins_([ACGT]+)", mutation_str)) {
        ins_match <- regexpr("ins_([ACGT]+)", mutation_str)
        ins_seq <- gsub("ins_", "", regmatches(mutation_str, ins_match))
      }

      # If we couldn't parse from mutation string, look it up from GD object
      if (is.na(ins_seq)) {
        events <- tryCatch(
          gd_events_table(gd),
          error = function(e) NULL
        )
        if (!is.null(events)) {
          matching_event <- events[events$type == "INS" &
                                  events$position == pos &
                                  events$seq_id == seq_id, ]
          if (nrow(matching_event) > 0 && !is.na(matching_event$ins_seq[1])) {
            ins_seq <- as.character(matching_event$ins_seq[1])
            # Add QC note that we looked it up
            row$qc_note <- paste0(
              if (!is.na(row$qc_note) && nzchar(row$qc_note)) paste0(row$qc_note, "; ") else "",
              "ins_seq from GD object"
            )
          }
        }
      }

      if (!is.na(ins_seq) && window_offset >= 1 && window_offset <= nchar(dna_window)) {
        before_ins <- substr(dna_window, 1, window_offset)
        after_ins <- if (window_offset < nchar(dna_window)) {
          substr(dna_window, window_offset + 1, nchar(dna_window))
        } else {
          ""
        }
        row$dna_alt <- paste0(before_ins, ins_seq, after_ins)
      }

    } else if (mut_type == "SUB") {
      # SUB: substitute sequence segment
      # For now, just copy dna_ref (SUB not fully implemented)
      row$dna_alt <- dna_window
      row$qc_note <- "SUB in intergenic - dna_alt = dna_ref (not fully implemented)"
    }
  }

  row
}

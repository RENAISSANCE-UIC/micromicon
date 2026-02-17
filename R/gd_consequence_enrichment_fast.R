#' Fast Mutation Consequence Enrichment
#'
#' @description
#' Optimized version of pm_enrich_consequences() with:
#' - Codon-level computation (don't translate full CDS)
#' - Gene-based batching (process all mutations per gene together)
#' - Minimal sequence operations
#'
#' @param gd A GenomeData object
#' @param pm_tbl A data.frame from `predict_mutations()`
#' @param flank Integer, flanking bases for intergenic regions (default: 50)
#' @param quiet Logical, suppress messages (default: FALSE)
#'
#' @return Enriched data.frame with consequence columns
#' @export
pm_enrich_consequences_fast <- function(gd, pm_tbl, flank = 50L, quiet = FALSE) {
  gd_assert(gd, "gd")
  stopifnot(is.data.frame(pm_tbl))

  if (!quiet) {
    cli::cli_inform(c(
      "i" = "pm_enrich_consequences_fast(): enriching {nrow(pm_tbl)} mutation{?s}"
    ))
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
  out$codon_new <- rep(NA_character_, n_rows)
  out$consequence <- rep(NA_character_, n_rows)
  out$region <- rep(NA_character_, n_rows)
  out$strand <- rep(NA_character_, n_rows)
  out$qc_note <- rep(NA_character_, n_rows)

  # Group mutations by gene for batch processing
  out$..row_id <- seq_len(n_rows)

  # Split by gene (NA genes handled separately)
  gene_groups <- split(out, out$gene, drop = FALSE)

  # Global caches (shared across all genes)
  gene_dna_cache <- new.env(parent = emptyenv())
  gene_aa_cache <- new.env(parent = emptyenv())
  gene_meta_cache <- new.env(parent = emptyenv())  # Store strand, feature info

  for (gene in names(gene_groups)) {
    group <- gene_groups[[gene]]

    if (is.na(gene) || !nzchar(gene) || nrow(group) == 0) {
      next
    }

    # Process all mutations for this gene together
    enriched_group <- .pf_enrich_gene_batch(gd, group, flank, quiet,
                                            gene_dna_cache, gene_aa_cache, gene_meta_cache)

    # Write back to output
    for (col in c("dna_ref", "dna_alt", "aa_ref", "aa_alt",
                  "codon_ref", "codon_alt", "codon_new", "consequence",
                  "region", "strand", "qc_note")) {
      out[enriched_group$..row_id, col] <- enriched_group[[col]]
    }
  }

  # Handle mutations without gene assignment (intergenic)
  intergenic_idx <- which(is.na(out$gene) | !nzchar(out$gene) | out$region == "intergenic")
  if (length(intergenic_idx) > 0) {
    for (i in intergenic_idx) {
      out[i, ] <- .pf_enrich_intergenic(gd, out[i, , drop = FALSE], flank, quiet)
    }
  }

  out$..row_id <- NULL

  if (!quiet) {
    n_coding <- sum(out$region == "coding", na.rm = TRUE)
    n_intergenic <- sum(out$region == "intergenic", na.rm = TRUE)
    cli::cli_inform(c(
      "v" = "Enriched {n_coding} coding and {n_intergenic} intergenic mutation{?s}"
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
  gene_clean <- gsub("\\s*[→←]\\s*$", "", gene)
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
  # Parse ref→alt (handle various arrow formats)
  mut_parts <- strsplit(gsub("→|->|&gt;", "|", mutation_str), "\\|")[[1]]
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
  row$codon_new <- codon_alt

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

  # Determine consequence
  if (!is.na(aa_alt_char)) {
    if (aa_alt_char == "*") {
      row$consequence <- "nonsense"
    } else if (aa_ref_char == aa_alt_char) {
      row$consequence <- "synonymous"
    } else {
      row$consequence <- "missense"
    }
  }

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

  # Quick frameshift check
  if (del_size %% 3 == 0) {
    row$consequence <- "inframe_deletion"
  } else {
    row$consequence <- "frameshift"
  }

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

  row$aa_alt <- aa_alt
  row
}


#' Fast INS Enrichment
#'
#' @keywords internal
.pf_enrich_ins_fast <- function(gd, row, dna_ref, aa_ref, gene_meta, quiet) {
  # Parse insertion info
  ins_info <- .pm_parse_insertion_annotation(as.character(row$annotation))

  if (is.na(ins_info$position) || is.na(ins_info$sequence)) {
    return(row)
  }

  # Apply insertion
  before_ins <- substr(dna_ref, 1, ins_info$position)
  after_ins <- substr(dna_ref, ins_info$position + 1, nchar(dna_ref))
  dna_alt <- paste0(before_ins, ins_info$sequence, after_ins)
  row$dna_alt <- dna_alt

  # Quick frameshift check
  ins_size <- nchar(ins_info$sequence)
  if (ins_size %% 3 == 0) {
    row$consequence <- "inframe_insertion"
  } else {
    row$consequence <- "frameshift"
  }

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
  row$codon_new <- NA_character_
  row$consequence <- NA_character_

  variant_type <- toupper(as.character(row$type))
  row$qc_note <- sprintf("Structural variant (%s) - alternate sequences not computed", variant_type)

  row
}


#' Fast Intergenic Enrichment
#'
#' @keywords internal
.pf_enrich_intergenic <- function(gd, row, flank, quiet) {
  # Extract DNA window
  pos <- as.integer(gsub(",", "", row$position))
  seq_id <- as.character(row$seq_id)

  start_pos <- max(1L, pos - flank)
  end_pos <- pos + flank

  dna_window <- tryCatch(
    get_roi_dna(gd, chrom = seq_id, start = start_pos, end = end_pos, strand = "+"),
    error = function(e) NA_character_
  )

  row$dna_ref <- dna_window
  row$region <- "intergenic"

  # Apply mutation to window if SNP
  if (toupper(row$type) == "SNP") {
    mutation_str <- as.character(row$mutation)
    mut_parts <- strsplit(gsub("→|->|&gt;", "|", mutation_str), "\\|")[[1]]
    if (length(mut_parts) == 2 && !is.na(dna_window)) {
      window_offset <- pos - start_pos + 1L
      alt_base <- mut_parts[2]
      dna_alt <- dna_window
      substr(dna_alt, window_offset, window_offset) <- alt_base
      row$dna_alt <- dna_alt
    }
  }

  row
}

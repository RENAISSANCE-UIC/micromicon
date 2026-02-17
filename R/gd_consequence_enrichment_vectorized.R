#' Ultra-Fast Vectorized Mutation Consequence Enrichment
#'
#' @description
#' Fully vectorized version with:
#' - Vectorized position → CDS mapping (all mutations in gene at once)
#' - Vectorized codon extraction
#' - Batch codon translation
#' - Vectorized string operations
#'
#' @param gd A GenomeData object
#' @param pm_tbl A data.frame from `predict_mutations()`
#' @param flank Integer, flanking bases for intergenic regions (default: 50)
#' @param quiet Logical, suppress messages (default: FALSE)
#'
#' @return Enriched data.frame with consequence columns
#' @keywords internal
pm_enrich_consequences_vectorized <- function(gd, pm_tbl, flank = 50L, quiet = FALSE) {
  gd_assert(gd, "gd")
  stopifnot(is.data.frame(pm_tbl))

  if (!quiet) {
    cli::cli_inform(c(
      "i" = "pm_enrich_consequences_vectorized(): enriching {nrow(pm_tbl)} mutation{?s}"
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

  # Global caches
  gene_dna_cache <- new.env(parent = emptyenv())
  gene_aa_cache <- new.env(parent = emptyenv())
  gene_meta_cache <- new.env(parent = emptyenv())

  # VECTORIZATION: Parse ALL positions and mutations upfront
  out$..genomic_pos <- as.integer(gsub(",|:.*", "", out$position))
  out$..ref_base <- gsub("(→|->|&gt;).*", "", out$mutation)
  out$..alt_base <- gsub(".*?(→|->|&gt;)", "", out$mutation)
  out$..row_id <- seq_len(n_rows)

  # Clean gene names (vectorized)
  out$..gene_clean <- gsub("\\s*[→←]\\s*$", "", out$gene)
  out$..gene_clean <- trimws(out$..gene_clean)

  # Group by gene for batch processing
  gene_groups <- split(out, out$..gene_clean, drop = FALSE)

  for (gene in names(gene_groups)) {
    group <- gene_groups[[gene]]

    if (is.na(gene) || !nzchar(gene) || nrow(group) == 0) {
      next
    }

    # Process entire gene batch with vectorization
    enriched_group <- .pv_enrich_gene_batch_vectorized(gd, group, flank, quiet,
                                                        gene_dna_cache, gene_aa_cache, gene_meta_cache)

    # Write back
    for (col in c("dna_ref", "dna_alt", "aa_ref", "aa_alt",
                  "codon_ref", "codon_alt", "codon_new", "consequence",
                  "region", "strand", "qc_note")) {
      out[enriched_group$..row_id, col] <- enriched_group[[col]]
    }
  }

  # Handle intergenic
  intergenic_idx <- which(is.na(out$gene) | !nzchar(out$gene) | out$region == "intergenic")
  if (length(intergenic_idx) > 0) {
    for (i in intergenic_idx) {
      out[i, ] <- .pf_enrich_intergenic(gd, out[i, , drop = FALSE], flank, quiet)
    }
  }

  # Clean up temp columns
  out$..genomic_pos <- NULL
  out$..ref_base <- NULL
  out$..alt_base <- NULL
  out$..row_id <- NULL
  out$..gene_clean <- NULL

  if (!quiet) {
    n_coding <- sum(out$region == "coding", na.rm = TRUE)
    n_intergenic <- sum(out$region == "intergenic", na.rm = TRUE)
    cli::cli_inform(c(
      "v" = "Enriched {n_coding} coding and {n_intergenic} mutation{?s}"
    ))
  }

  out
}


#' Vectorized Gene Batch Processing
#'
#' @keywords internal
.pv_enrich_gene_batch_vectorized <- function(gd, group, flank, quiet,
                                              gene_dna_cache, gene_aa_cache, gene_meta_cache) {
  gene <- as.character(group$..gene_clean[1])

  # Get gene metadata once
  gene_meta <- .pf_get_gene_metadata(gd, gene, gene_meta_cache)

  if (is.null(gene_meta) || is.na(gene_meta$strand)) {
    return(group)
  }

  # Get ref sequences once
  dna_ref <- .pm_cached_get_gene_dna(gd, gene, gene_dna_cache)
  aa_ref <- .pm_cached_translate_gene(dna_ref, gene, gene_aa_cache)

  if (is.na(dna_ref) || is.na(aa_ref)) {
    return(group)
  }

  # Fill shared fields (vectorized!)
  group$dna_ref <- dna_ref
  group$aa_ref <- aa_ref
  group$region <- "coding"
  group$strand <- gene_meta$strand

  # VECTORIZATION: Map ALL genomic positions to CDS positions at once
  genomic_positions <- group$..genomic_pos

  if (gene_meta$strand == "-") {
    # Minus strand: count from end
    cds_positions <- gene_meta$end - genomic_positions + 1L
  } else {
    # Plus strand: count from start
    cds_positions <- genomic_positions - gene_meta$start + 1L
  }

  # Filter out-of-range positions
  in_range <- !is.na(cds_positions) & cds_positions >= 1 & cds_positions <= nchar(dna_ref)
  group$..cds_pos <- ifelse(in_range, cds_positions, NA_integer_)

  # VECTORIZATION: Compute ALL codon boundaries at once
  aa_indices <- ceiling(group$..cds_pos / 3)
  codon_starts <- (aa_indices - 1) * 3 + 1
  codon_ends <- pmin(codon_starts + 2, nchar(dna_ref))

  # VECTORIZATION: Extract ALL reference codons at once (vectorized substr!)
  group$codon_ref <- ifelse(
    !is.na(group$..cds_pos),
    substr(rep(dna_ref, nrow(group)), codon_starts, codon_ends),
    NA_character_
  )

  # Process by mutation type
  snp_idx <- which(toupper(group$type) == "SNP" & !is.na(group$..cds_pos))
  del_idx <- which(toupper(group$type) == "DEL")
  ins_idx <- which(toupper(group$type) == "INS")
  structural_idx <- which(toupper(group$type) %in% c("MOB", "AMP", "CON", "INV"))

  # VECTORIZED SNP PROCESSING
  if (length(snp_idx) > 0) {
    snp_group <- group[snp_idx, ]

    # Vectorized: offset within codon
    codon_offsets <- snp_group$..cds_pos - codon_starts[snp_idx] + 1

    # Vectorized: apply complement for minus strand
    alt_bases_cds <- ifelse(
      gene_meta$strand == "-",
      sapply(snp_group$..alt_base, gd_complement_base),
      snp_group$..alt_base
    )

    # Build alternate codons (handle properly in vectorized way)
    codons_alt <- character(length(snp_group$codon_ref))
    for (i in seq_along(snp_group$codon_ref)) {
      codon_ref_i <- snp_group$codon_ref[i]
      if (!is.na(codon_ref_i) && !is.na(codon_offsets[i]) && !is.na(alt_bases_cds[i])) {
        codon_alt_i <- codon_ref_i
        substr(codon_alt_i, codon_offsets[i], codon_offsets[i]) <- alt_bases_cds[i]
        codons_alt[i] <- codon_alt_i
      } else {
        codons_alt[i] <- NA_character_
      }
    }

    snp_group$codon_alt <- codons_alt
    snp_group$codon_new <- codons_alt

    # Translate each codon individually (vectorized with sapply)
    snp_group$..aa_alt_char <- sapply(codons_alt, function(codon) {
      if (is.na(codon) || nchar(codon) != 3) return(NA_character_)
      aa <- tryCatch(
        translate_dna(codon, frame = 1, genetic_code = "11", .internal = TRUE),
        error = function(e) NA_character_
      )
      if (!is.na(aa) && nchar(aa) > 0) substr(aa, 1, 1) else NA_character_
    })

    # Vectorized: build full mutated CDS for each SNP
    for (i in seq_len(nrow(snp_group))) {
      if (!is.na(snp_group$..cds_pos[i])) {
        dna_alt_i <- dna_ref
        substr(dna_alt_i, snp_group$..cds_pos[i], snp_group$..cds_pos[i]) <- alt_bases_cds[i]
        snp_group$dna_alt[i] <- dna_alt_i

        # Build aa_alt by copying aa_ref and changing one position
        aa_alt_i <- aa_ref
        aa_idx <- aa_indices[snp_idx[i]]
        if (!is.na(snp_group$..aa_alt_char[i])) {
          substr(aa_alt_i, aa_idx, aa_idx) <- snp_group$..aa_alt_char[i]
        }
        snp_group$aa_alt[i] <- aa_alt_i

        # Determine consequence
        aa_ref_char <- substr(aa_ref, aa_idx, aa_idx)
        aa_alt_char <- snp_group$..aa_alt_char[i]

        if (!is.na(aa_alt_char)) {
          if (aa_alt_char == "*") {
            snp_group$consequence[i] <- "nonsense"
          } else if (aa_ref_char == aa_alt_char) {
            snp_group$consequence[i] <- "synonymous"
          } else {
            snp_group$consequence[i] <- "missense"
          }
        }
      }
    }

    # Copy back only the enrichment columns
    for (col in c("codon_alt", "codon_new", "dna_alt", "aa_alt", "consequence")) {
      group[snp_idx, col] <- snp_group[[col]]
    }
  }

  # DEL/INS processing (still per-mutation due to variable-length sequences)
  for (i in del_idx) {
    group[i, ] <- .pf_enrich_del_fast(gd, group[i, , drop = FALSE],
                                      dna_ref, aa_ref, gene_meta, quiet)
  }

  for (i in ins_idx) {
    group[i, ] <- .pf_enrich_ins_fast(gd, group[i, , drop = FALSE],
                                      dna_ref, aa_ref, gene_meta, quiet)
  }

  for (i in structural_idx) {
    group[i, ] <- .pf_enrich_structural_fast(group[i, , drop = FALSE])
  }

  # Clean temp columns
  group$..cds_pos <- NULL
  group$..aa_alt_char <- NULL

  group
}

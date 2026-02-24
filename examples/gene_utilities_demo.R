#!/usr/bin/env Rscript
#
# Demo: Gene Utility Functions in micromicon
#
# This script demonstrates the new gene utility functions for working with
# genome_entity objects.

library(micromicon)

# Assuming you have a genome loaded (GenBank or GFF3+FASTA)
# entity <- read_genome("path/to/genome.gbk")
# OR
# entity <- read_genome(gff = "path/to/genome.gff3", fasta = "path/to/genome.fasta")

# For this demo, we'll use a test file if available
if (file.exists("inst/extdata/test.gbk")) {
  entity <- read_genome("inst/extdata/test.gbk")
} else {
  stop("Test file not found. Please provide a genome file.")
}

cat("=== Gene Utility Functions Demo ===\n\n")

# Get all CDS features
cds_features <- features(entity, type = "CDS")
cat(sprintf("Found %d CDS features\n\n", nrow(cds_features)))

# Select first CDS for demo
if (nrow(cds_features) > 0) {
  first_cds <- cds_features[1, ]
  cat(sprintf("Working with CDS: %s\n", first_cds$locus_tag %||% first_cds$gene %||% "unnamed"))
  cat(sprintf("  Location: %s:%d-%d (%s strand)\n\n",
              first_cds$seqname, first_cds$start, first_cds$end, first_cds$strand))

  # 1. Get DNA sequence (5'-to-3' canonical orientation)
  cat("1. Get gene DNA sequence:\n")
  dna <- get_gene_dna(entity, 1)
  cat(sprintf("   Length: %d bp\n", nchar(dna)))
  cat(sprintf("   First 60 bp: %s...\n\n", substr(dna, 1, 60)))

  # 2. Get amino acid sequence (translated)
  cat("2. Get gene amino acid sequence:\n")
  aa <- get_gene_aa(entity, 1)
  cat(sprintf("   Length: %d aa\n", nchar(aa)))
  cat(sprintf("   First 20 aa: %s...\n\n", substr(aa, 1, 20)))

  # 3. Get arbitrary region of DNA
  cat("3. Get DNA for a region of interest (ROI):\n")
  roi_start <- first_cds$start
  roi_end <- first_cds$start + 99
  roi <- get_roi_dna(entity, contig = first_cds$seqname,
                     start = roi_start, end = roi_end, strand = "+")
  cat(sprintf("   Region: %s:%d-%d\n", first_cds$seqname, roi_start, roi_end))
  cat(sprintf("   Sequence: %s\n\n", roi))

  # 4. Translate DNA with different frames
  cat("4. Translate DNA with different reading frames:\n")
  test_dna <- "ATGAAATTTCCCGGGTTTAAATAG"
  cat(sprintf("   Test DNA: %s\n", test_dna))
  for (frame in 1:3) {
    aa_frame <- translate_dna(test_dna, frame = frame, genetic_code = "11")
    cat(sprintf("   Frame %d: %s\n", frame, aa_frame))
  }
  cat("\n")

  # 5. Map genomic position to CDS position
  cat("5. Map genomic position to CDS position:\n")
  genomic_pos <- first_cds$start + 10
  cds_pos <- map_genomic_to_cds(entity, 1, genomic_pos)
  cat(sprintf("   Genomic position %d -> CDS position %d\n", genomic_pos, cds_pos))

  # 6. Map CDS position back to genomic position
  genomic_back <- map_cds_to_genomic(entity, 1, cds_pos)
  cat(sprintf("   CDS position %d -> Genomic position %d\n\n", cds_pos, genomic_back))

  # 7. Get gene information
  cat("7. Get gene information:\n")
  info <- gene_info(entity, 1)
  cat(sprintf("   Chromosome: %s\n", info$chrom))
  cat(sprintf("   Coordinates: %d-%d (%s strand)\n", info$start, info$end, info$strand))
  cat(sprintf("   Length: %d bp\n", info$length))
  cat(sprintf("   Type: %s\n", info$type))
  cat(sprintf("   Frame: %d\n", info$frame))
  if (!is.na(info$gene)) {
    cat(sprintf("   Gene: %s\n", info$gene))
  }
  if (!is.na(info$locus_tag)) {
    cat(sprintf("   Locus tag: %s\n", info$locus_tag))
  }
  if (!is.na(info$product)) {
    cat(sprintf("   Product: %s\n", info$product))
  }
  cat("\n")

  # 8. Validate a variant
  cat("8. Validate variant in gene:\n")
  seqs <- get_contig_sequences(entity)
  ref_base <- substr(seqs[[first_cds$seqname]], first_cds$start, first_cds$start)
  result <- validate_variant_in_gene(entity, 1, first_cds$start, ref_base)
  cat(sprintf("   Position: %d, Ref: %s\n", first_cds$start, ref_base))
  cat(sprintf("   Valid: %s\n", result))
  cat(sprintf("   Message: %s\n\n", attr(result, "message")))

  # Invalid variant (wrong ref base)
  wrong_base <- if (ref_base == "A") "T" else "A"
  result_invalid <- validate_variant_in_gene(entity, 1, first_cds$start, wrong_base)
  cat(sprintf("   Position: %d, Ref: %s (wrong)\n", first_cds$start, wrong_base))
  cat(sprintf("   Valid: %s\n", result_invalid))
  cat(sprintf("   Message: %s\n\n", attr(result_invalid, "message")))
}

cat("=== Demo Complete ===\n")

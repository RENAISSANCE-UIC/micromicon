#!/usr/bin/env Rscript

# Example: Enriching Mutation Consequences
# ==========================================
# This example demonstrates how to use pm_enrich_consequences() to add
# molecular consequence annotations to predicted mutations.

library(micromicon)

# ============================================================================
# Step 1: Load GenomeData object with breseq results
# ============================================================================
# Assuming you have a .gd file and reference genome:
#
# gd <- read_genomedata(
#   gd_path = "path/to/output/output.gd",
#   ref_dir = "path/to/reference/"
# )


# ============================================================================
# Step 2: Get predicted mutations table
# ============================================================================
# mutations <- predict_mutations(gd)
#
# This returns a table like:
# | evidence | type | seq_id    | position | mutation | freq  | annotation | gene  | description        |
# |----------|------|-----------|----------|----------|-------|------------|-------|--------------------|
# | RA       | SNP  | NC_007795 | 1,234    | A→C      | 100%  | coding     | ampC  | beta-lactamase     |
# | RA       | SNP  | NC_007795 | 5,678    | G→T      | 45.2% | intergenic | –     | glnS → / ← ileS    |


# ============================================================================
# Step 3: Enrich with consequences
# ============================================================================
# enriched <- pm_enrich_consequences(gd, mutations, flank = 50)
#
# This adds the following columns:
# - dna_ref: Reference DNA sequence (full CDS for coding; window for intergenic)
# - dna_alt: Alternate DNA sequence with mutation applied
# - aa_ref: Reference amino acid sequence (coding only)
# - aa_alt: Alternate amino acid sequence (coding only)
# - codon_ref: Reference codon at mutation site (coding only)
# - codon_alt: Alternate codon at mutation site (coding only)
# - consequence: "synonymous", "missense", or "nonsense" (coding only)
# - region: "coding" or "intergenic"
# - strand: Gene strand ("+" or "-", coding only)


# ============================================================================
# Step 4: Filter and analyze consequences
# ============================================================================
# # Get all missense mutations
# missense <- subset(enriched, consequence == "missense")
#
# # Get all nonsense mutations
# nonsense <- subset(enriched, consequence == "nonsense")
#
# # Get all synonymous mutations
# synonymous <- subset(enriched, consequence == "synonymous")
#
# # Check amino acid changes
# missense[, c("gene", "position", "mutation", "codon_ref", "codon_alt",
#              "aa_ref", "aa_alt", "consequence")]


# ============================================================================
# Example Output
# ============================================================================
# For a coding SNP (A→C at position 1234 in gene ampC):
#
# | column      | value                                           |
# |-------------|-------------------------------------------------|
# | dna_ref     | ATGGCACAAGTCAT... (full CDS, 1200 bp)          |
# | dna_alt     | ATGGCACAAGTCCT... (mutated CDS)                |
# | aa_ref      | MAQVI... (protein sequence, 400 aa)             |
# | aa_alt      | MAQVP... (mutated protein)                      |
# | codon_ref   | ATG                                             |
# | codon_alt   | CTG                                             |
# | consequence | missense                                        |
# | region      | coding                                          |
# | strand      | +                                               |
#
# For an intergenic SNP (G→T at position 5678):
#
# | column      | value                                           |
# |-------------|-------------------------------------------------|
# | dna_ref     | ...ATCGATCG[G]TCGATCG... (100 bp window)       |
# | dna_alt     | ...ATCGATCG[T]TCGATCG... (mutated window)      |
# | aa_ref      | NA                                              |
# | aa_alt      | NA                                              |
# | codon_ref   | NA                                              |
# | codon_alt   | NA                                              |
# | consequence | NA                                              |
# | region      | intergenic                                      |
# | strand      | NA                                              |


# ============================================================================
# Implementation Details
# ============================================================================
# The function leverages existing micromicon utilities:
# - gd_verify_snp() - Computes SNP effects for coding regions
# - gd_locate() - Classifies regions as coding vs intergenic
# - get_gene_dna(), get_gene_aa() - Extract gene sequences
# - get_roi_dna() - Extract arbitrary genomic windows
# - translate_dna() - Translate CDS to protein
# - reverse_complement() - Handle strand orientation
#
# This ensures consistency with other micromicon functions and minimizes
# code duplication.


# ============================================================================
# Future Extensions
# ============================================================================
# Currently supports SNP mutations only. Future versions will handle:
# - DEL (deletions)
# - INS (insertions)
# - SUB (substitutions)
# - Complex mutations

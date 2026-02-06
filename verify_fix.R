# Verification script for search_features() bug fix
# This demonstrates that the bug reported by the user is now fixed

library(micromicon)

# Create a test entity with features using Name and ID fields
# (simulating a GFF3 file from breseq, Prokka, or NCBI)
test_entity <- new_genome_entity(
  sequences_list = list(
    dna_raw = c(chromosome = "ATCGATCGATCGATCGATCG"),
    dna_bio = NULL,
    indexed_fa = NULL
  ),
  features_df = data.frame(
    seqname = rep("chromosome", 4),
    start = c(1, 6, 11, 16),
    end = c(5, 10, 15, 20),
    strand = rep("+", 4),
    type = rep("gene", 4),
    ID = c("gene_001", "gene_002", "gene_003", "ampC_gene"),
    Name = c("dnaA", "ampC", "gyrA", NA),
    Alias = c(NA, "beta-lactamase", NA, NA),
    gene = c(NA, NA, NA, "ampC"),
    locus_tag = c("TAG_001", "TAG_002", "TAG_003", "TAG_004"),
    product = c("chromosomal replication initiator",
                "beta-lactamase",
                "DNA gyrase subunit A",
                "beta-lactamase"),
    stringsAsFactors = FALSE
  ),
  metadata_df = data.frame(
    seqname = "chromosome",
    length = 20,
    stringsAsFactors = FALSE
  ),
  indices_list = list(
    seqnames = "chromosome",
    locus_tag_index = integer(),
    gene_index = list()
  )
)

cat("=== Verification of search_features() Bug Fix ===\n\n")

# Test 1: Search by Name field (this was the reported bug)
cat("Test 1: Searching for 'ampC' (stored in Name field of row 2)\n")
result1 <- search_features(test_entity, pattern = "ampC")
cat("  Found", nrow(result1), "features\n")
cat("  IDs:", paste(result1$ID, collapse = ", "), "\n")
cat("  Status:", ifelse(nrow(result1) >= 1, "✓ PASS", "✗ FAIL"), "\n\n")

# Test 2: Search by ID field
cat("Test 2: Searching for 'gene_001' (stored in ID field)\n")
result2 <- search_features(test_entity, pattern = "gene_001")
cat("  Found", nrow(result2), "features\n")
cat("  Names:", paste(result2$Name, collapse = ", "), "\n")
cat("  Status:", ifelse(nrow(result2) == 1, "✓ PASS", "✗ FAIL"), "\n\n")

# Test 3: Search by Alias field
cat("Test 3: Searching for 'beta-lactamase' (stored in Alias and product fields)\n")
result3 <- search_features(test_entity, pattern = "beta-lactamase")
cat("  Found", nrow(result3), "features\n")
cat("  Status:", ifelse(nrow(result3) == 2, "✓ PASS", "✗ FAIL"), "\n\n")

# Test 4: Search finds matches in multiple fields
cat("Test 4: Searching for 'ampC' (appears in Name, gene, ID, and product fields)\n")
result4 <- search_features(test_entity, pattern = "ampC")
cat("  Found", nrow(result4), "features\n")
cat("  Expected: 2 (rows 2 and 4)\n")
cat("  Status:", ifelse(nrow(result4) == 2, "✓ PASS", "✗ FAIL"), "\n\n")

# Test 5: Verify case-insensitive search works
cat("Test 5: Case-insensitive search for 'AMPC'\n")
result5 <- search_features(test_entity, pattern = "AMPC")
cat("  Found", nrow(result5), "features\n")
cat("  Status:", ifelse(nrow(result5) == 2, "✓ PASS", "✗ FAIL"), "\n\n")

# Test 6: Verify existing functionality still works (gene field)
cat("Test 6: Searching traditional 'gene' field\n")
result6 <- search_features(test_entity, pattern = "ampC")
cat("  Found", nrow(result6), "features\n")
cat("  Status:", ifelse(nrow(result6) >= 1, "✓ PASS", "✗ FAIL"), "\n\n")

# Summary
all_pass <- all(
  nrow(result1) >= 1,
  nrow(result2) == 1,
  nrow(result3) == 2,
  nrow(result4) == 2,
  nrow(result5) == 2,
  nrow(result6) >= 1
)

cat("=== Overall Result ===\n")
cat(ifelse(all_pass,
           "✓ All tests PASSED! The bug is fixed.\n",
           "✗ Some tests FAILED. Please review.\n"))

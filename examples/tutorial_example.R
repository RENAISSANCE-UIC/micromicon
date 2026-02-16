#!/usr/bin/env Rscript
#
# Tutorial Example: Using pm_enrich_consequences() and %||%
#
# This script demonstrates the new features with working examples.
#

library(micromicon)

cat("=== Micromicon New Features Tutorial ===\n\n")

# ==============================================================================
# Part 1: The %||% Operator (pronounced "or-else")
# ==============================================================================

cat("PART 1: The %||% Operator\n")
cat("-------------------------\n\n")

# Example 1: Basic usage
cat("Example 1: Basic NULL handling\n")
x <- NULL
result1 <- x %||% "default"
cat("  NULL %||% 'default' =", result1, "\n")

x <- "value"
result2 <- x %||% "default"
cat("  'value' %||% 'default' =", result2, "\n\n")

# Example 2: Function parameters
cat("Example 2: Function with optional parameters\n")
greet <- function(name = NULL, greeting = NULL) {
  name <- name %||% "friend"
  greeting <- greeting %||% "Hello"
  paste(greeting, name, "!")
}

cat("  greet() =", greet(), "\n")
cat("  greet('Alice') =", greet("Alice"), "\n")
cat("  greet('Bob', 'Hi') =", greet("Bob", "Hi"), "\n\n")

# Example 3: List extraction with fallbacks
cat("Example 3: Extracting from lists with defaults\n")
config <- list(threads = 4, memory = NULL, timeout = 30)

threads <- config$threads %||% 1
memory <- config$memory %||% "8GB"
debug <- config$debug %||% FALSE

cat("  threads:", threads, "\n")
cat("  memory:", memory, "\n")
cat("  debug:", debug, "\n\n")

# Example 4: Chaining multiple fallbacks
cat("Example 4: Chaining fallbacks (first non-NULL wins)\n")
feature <- list(
  gene = NULL,
  locus_tag = NULL,
  protein_id = "YP_001234",
  product = "hypothetical protein"
)

gene_name <- feature$gene %||%
             feature$locus_tag %||%
             feature$protein_id %||%
             "unknown"

cat("  Best identifier found:", gene_name, "\n\n")

# ==============================================================================
# Part 2: Mutation Consequence Enrichment
# ==============================================================================

cat("PART 2: Mutation Consequence Enrichment\n")
cat("----------------------------------------\n\n")

# Create a minimal mock dataset for demonstration
cat("Creating mock mutation data for demonstration...\n\n")

# Note: In real usage, you would do:
# gd <- read_genomedata("file.gd", ref_dir = "reference/")
# mutations <- predict_mutations(gd)
# enriched <- pm_enrich_consequences(gd, mutations)

# For this tutorial, we'll show the expected structure:
cat("Expected workflow:\n")
cat("  1. gd <- read_genomedata('data.gd', ref_dir = 'reference/')\n")
cat("  2. mutations <- predict_mutations(gd)\n")
cat("  3. enriched <- pm_enrich_consequences(gd, mutations)\n\n")

cat("Expected input (mutations table from predict_mutations):\n")
example_mutations <- data.frame(
  type = c("SNP", "SNP", "SNP", "SNP"),
  seq_id = c("NC_000913", "NC_000913", "NC_000913", "NC_000913"),
  position = c("12345", "67890", "45678:1", "11111"),
  mutation = c("A→C", "G->T", "C>A", "T-&gt;G"),
  gene = c("dnaK", "rpoB", "gyrA", "ftsZ"),
  annotation = c(
    "coding (423/1800 nt)",
    "coding (1234/3537 nt)",
    NA,
    "coding (1/2592 nt)"
  ),
  stringsAsFactors = FALSE
)
print(example_mutations)

cat("\n\nExpected output columns after enrichment:\n")
output_cols <- c(
  "dna_ref", "dna_alt",           # DNA sequences
  "aa_ref", "aa_alt",             # Amino acid sequences
  "codon_ref", "codon_alt",       # Codons at mutation site
  "codon_new",                    # Alias for codon_alt
  "consequence",                  # synonymous/missense/nonsense
  "region",                       # coding/intergenic
  "strand",                       # +/-
  "qc_note"                       # Quality control notes
)
cat(paste("  -", output_cols, collapse = "\n"), "\n\n")

# ==============================================================================
# Part 3: Understanding the New Features
# ==============================================================================

cat("PART 3: New Features Explained\n")
cat("------------------------------\n\n")

cat("✨ Feature 1: QC Notes\n")
cat("  The qc_note column tracks quality control issues:\n")
cat("  - 'SNP has non-zero offset (1) - unusual'\n")
cat("  - 'Normalized alt start GTG->M'\n")
cat("  - 'Used gd_locate fallback'\n")
cat("  - 'Position parse failed'\n\n")

cat("✨ Feature 2: Alt-Start Normalization\n")
cat("  Alternative start codons are normalized to M:\n")
cat("  - GTG (Valine) → M (Methionine)\n")
cat("  - TTG (Leucine) → M\n")
cat("  - CTG (Leucine) → M\n")
cat("  This matches how bacteria actually translate these!\n\n")

cat("✨ Feature 3: HTML-Escaped Formats\n")
cat("  Now supports arrows from web/database exports:\n")
cat("  - A→C (Unicode)\n")
cat("  - A->C (ASCII)\n")
cat("  - A>C (simple)\n")
cat("  - A-&gt;C (HTML escaped)\n")
cat("  - A&gt;C (HTML escaped)\n\n")

cat("✨ Feature 4: Annotation-First Geometry\n")
cat("  Performance optimization:\n")
cat("  - Parses 'coding (i/L nt)' from annotation field first\n")
cat("  - Only calls gd_locate() if annotation is missing\n")
cat("  - Much faster for large mutation tables!\n\n")

cat("✨ Feature 5: Codon Naming Consistency\n")
cat("  Both codon_alt and codon_new are provided:\n")
cat("  - codon_alt: Original column name\n")
cat("  - codon_new: Alias matching gd_verify_snp() output\n")
cat("  They contain identical values!\n\n")

# ==============================================================================
# Part 4: Practical Use Cases
# ==============================================================================

cat("PART 4: Practical Use Cases\n")
cat("---------------------------\n\n")

cat("Use Case 1: Finding Resistance Mutations\n")
cat("  resistance_genes <- c('rpoB', 'gyrA', 'gyrB')\n")
cat("  candidates <- enriched[\n")
cat("    enriched$gene %in% resistance_genes &\n")
cat("    enriched$consequence == 'missense',\n")
cat("  ]\n\n")

cat("Use Case 2: Quality Control Check\n")
cat("  # Find mutations with QC issues\n")
cat("  issues <- enriched[!is.na(enriched$qc_note), ]\n")
cat("  table(issues$qc_note)\n\n")

cat("Use Case 3: Comparing Synonymous vs Non-Synonymous\n")
cat("  # Count consequence types\n")
cat("  table(enriched$consequence)\n")
cat("  \n")
cat("  # Calculate dN/dS ratio\n")
cat("  n_nonsyn <- sum(enriched$consequence %in% c('missense', 'nonsense'))\n")
cat("  n_syn <- sum(enriched$consequence == 'synonymous')\n")
cat("  dn_ds_ratio <- n_nonsyn / n_syn\n\n")

cat("Use Case 4: Exporting for Further Analysis\n")
cat("  # Save enriched mutations\n")
cat("  write.csv(enriched, 'mutations_enriched.csv', row.names = FALSE)\n")
cat("  \n")
cat("  # Or filter to interesting ones\n")
cat("  interesting <- enriched[enriched$consequence == 'nonsense', ]\n")
cat("  write.csv(interesting, 'nonsense_mutations.csv', row.names = FALSE)\n\n")

# ==============================================================================
# Part 5: Common Patterns with %||%
# ==============================================================================

cat("PART 5: Common Patterns with %||%\n")
cat("----------------------------------\n\n")

cat("Pattern 1: Optional flank parameter\n")
my_extract <- function(gd, gene, flank = NULL) {
  flank <- flank %||% 50  # Default to 50 if not provided
  cat("  Using flank =", flank, "\n")
  # ... rest of function
}
my_extract(NULL, "dnaK")           # Uses default 50
my_extract(NULL, "dnaK", 100)      # Uses provided 100
cat("\n")

cat("Pattern 2: Flexible output control\n")
my_analyze <- function(data, quiet = NULL, verbose = NULL) {
  # Convert to boolean with defaults
  quiet <- quiet %||% FALSE
  verbose <- verbose %||% !quiet  # Verbose is opposite of quiet by default

  cat("  quiet =", quiet, ", verbose =", verbose, "\n")
}
my_analyze(NULL)                    # quiet=FALSE, verbose=TRUE
my_analyze(NULL, quiet = TRUE)      # quiet=TRUE, verbose=FALSE
my_analyze(NULL, verbose = FALSE)   # quiet=FALSE, verbose=FALSE
cat("\n")

cat("Pattern 3: Configuration with fallbacks\n")
run_analysis <- function(config = list()) {
  threads <- config$threads %||% 1
  memory <- config$memory %||% "4GB"
  timeout <- config$timeout %||% 3600

  cat("  Configuration:\n")
  cat("    threads:", threads, "\n")
  cat("    memory:", memory, "\n")
  cat("    timeout:", timeout, "seconds\n")
}
run_analysis()                                    # All defaults
run_analysis(list(threads = 8, memory = "16GB"))  # Partial override
cat("\n")

# ==============================================================================
# Summary
# ==============================================================================

cat("=== Summary ===\n\n")
cat("You've learned:\n")
cat("  ✓ How to use the %||% operator (say 'or-else'!)\n")
cat("  ✓ What pm_enrich_consequences() does\n")
cat("  ✓ The new QC tracking features\n")
cat("  ✓ How alt-start codons are normalized\n")
cat("  ✓ Common patterns and use cases\n\n")

cat("Next steps:\n")
cat("  1. Read the full tutorial: TUTORIAL_NEW_FEATURES.md\n")
cat("  2. Try with your own data\n")
cat("  3. Check the documentation: ?pm_enrich_consequences\n")
cat("  4. Use %||% in your own functions!\n\n")

cat("Remember: It's 'or-else', NOT 'grapes-or-or-grapes'! 🍇🚫\n")

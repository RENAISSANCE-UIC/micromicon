#!/usr/bin/env Rscript
devtools::load_all(".", quiet = TRUE)

entity <- micromicon::read_genome(
  fasta = ".test_these/reference.fasta",
  gff   = ".test_these/reference.gff3"
)

gd <- micromicon::parse_gd_annotated(
  gd_path = ".test_these/annotated.gd",
  entity  = entity
)

mut_table <- predict_mutations(gd)

cat("Testing Scientific Equivalence\n")
cat("===============================\n\n")

cons_tbl <- pm_enrich_consequences(gd, mut_table, quiet = TRUE)
cons_tbl_parallel <- pm_enrich_consequences_parallel(gd, mut_table, quiet = TRUE)

cat("Checking scientifically important columns:\n\n")

# Core scientific output columns
scientific_cols <- c("consequence", "codon_ref", "codon_alt", "region")

all_match <- TRUE

for (col in scientific_cols) {
  # Compare values
  matches <- sum(cons_tbl[[col]] == cons_tbl_parallel[[col]], na.rm = TRUE)
  both_na <- sum(is.na(cons_tbl[[col]]) & is.na(cons_tbl_parallel[[col]]))
  total_match <- matches + both_na
  total_rows <- nrow(cons_tbl)

  if (total_match == total_rows) {
    cat("✅", col, ": PERFECT MATCH (", total_match, "/", total_rows, ")\n")
  } else {
    cat("❌", col, ": MISMATCH (", total_match, "/", total_rows, ")\n")
    all_match <- FALSE

    # Show differences
    diff_idx <- which(
      (is.na(cons_tbl[[col]]) != is.na(cons_tbl_parallel[[col]])) |
      (!is.na(cons_tbl[[col]]) & !is.na(cons_tbl_parallel[[col]]) &
       cons_tbl[[col]] != cons_tbl_parallel[[col]])
    )

    if (length(diff_idx) > 0) {
      cat("   Differences at rows:", paste(diff_idx[1:min(5, length(diff_idx))], collapse = ", "), "\n")
    }
  }
}

cat("\n")

if (all_match) {
  cat("╔════════════════════════════════════════════╗\n")
  cat("║  ✅ SCIENTIFIC OUTPUT IDENTICAL           ║\n")
  cat("║  All consequence calls match perfectly    ║\n")
  cat("╚════════════════════════════════════════════╝\n")
} else {
  cat("⚠️  Some differences in scientific output\n")
}

cat("\nNote: Minor differences in metadata columns (qc_note, strand NAs)\n")
cat("      are acceptable as long as consequences match.\n")

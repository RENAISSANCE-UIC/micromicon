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

cat("Testing Exact Equivalence\n")
cat("==========================\n\n")

cat("Running pm_enrich_consequences()...\n")
cons_tbl <- pm_enrich_consequences(gd, mut_table, quiet = TRUE)

cat("Running pm_enrich_consequences_parallel()...\n")
cons_tbl_parallel <- pm_enrich_consequences_parallel(gd, mut_table, quiet = TRUE)

cat("\nComparing results with all.equal()...\n\n")

result <- all.equal(cons_tbl, cons_tbl_parallel)

if (isTRUE(result)) {
  cat("✅ SUCCESS: all.equal() returns TRUE\n")
  cat("   Both functions produce identical output\n\n")
} else {
  cat("❌ FAILURE: Differences found:\n")
  print(result)
  cat("\n")

  # Show detailed comparison
  cat("Detailed comparison:\n")
  cat("  Dimensions match:", identical(dim(cons_tbl), dim(cons_tbl_parallel)), "\n")
  cat("  Column names match:", identical(names(cons_tbl), names(cons_tbl_parallel)), "\n")
  cat("  Row order match:", identical(cons_tbl$position, cons_tbl_parallel$position), "\n")

  # Check consequence column specifically
  cat("\nConsequence comparison:\n")
  cat("  Regular:\n")
  print(table(cons_tbl$consequence, useNA = "ifany"))
  cat("  Parallel:\n")
  print(table(cons_tbl_parallel$consequence, useNA = "ifany"))
}

cat("\nDone.\n")

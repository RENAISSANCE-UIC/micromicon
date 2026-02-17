#!/usr/bin/env Rscript
devtools::load_all(".", quiet = TRUE)

cat("Testing with .now_test_these dataset\n")
cat("=====================================\n\n")

entity <- micromicon::read_genome(
  fasta = ".now_test_these/reference.fasta",
  gff   = ".now_test_these/reference.gff3"
)

gd <- micromicon::parse_gd_annotated(
  gd_path = ".now_test_these/annotated.gd",
  entity  = entity
)

mut_table <- predict_mutations(gd)
cat("Mutations:", nrow(mut_table), "\n\n")

cat("Running pm_enrich_consequences()...\n")
cons_tbl <- pm_enrich_consequences(gd, mut_table, quiet = TRUE)

cat("Running pm_enrich_consequences_parallel()...\n")
cons_tbl_parallel <- pm_enrich_consequences_parallel(gd, mut_table, quiet = TRUE)

cat("\n")
cat("Results:\n")
cat("--------\n")
cat("Regular:  ", sum(cons_tbl$region == "coding", na.rm = TRUE), " coding, ",
    sum(cons_tbl$region == "intergenic", na.rm = TRUE), " intergenic\n", sep = "")
cat("Parallel: ", sum(cons_tbl_parallel$region == "coding", na.rm = TRUE), " coding, ",
    sum(cons_tbl_parallel$region == "intergenic", na.rm = TRUE), " intergenic\n\n", sep = "")

# Test all.equal
cat("Testing all.equal()...\n")
result <- all.equal(cons_tbl, cons_tbl_parallel)

if (isTRUE(result)) {
  cat("\n✅ SUCCESS: all.equal() returns TRUE\n")
} else {
  cat("\n❌ DIFFERENCES FOUND:\n")
  print(result)
  cat("\n")
}

# Check scientific columns
cat("\nScientific columns:\n")
scientific_cols <- c("consequence", "codon_ref", "codon_alt", "region")

for (col in scientific_cols) {
  matches <- sum(cons_tbl[[col]] == cons_tbl_parallel[[col]], na.rm = TRUE)
  both_na <- sum(is.na(cons_tbl[[col]]) & is.na(cons_tbl_parallel[[col]]))
  total_match <- matches + both_na

  if (total_match == nrow(cons_tbl)) {
    cat("  ✅ ", col, ": ", total_match, "/", nrow(cons_tbl), "\n", sep = "")
  } else {
    cat("  ❌ ", col, ": ", total_match, "/", nrow(cons_tbl), "\n", sep = "")
  }
}

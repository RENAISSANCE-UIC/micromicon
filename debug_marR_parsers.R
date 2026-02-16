#!/usr/bin/env Rscript
# Debug why marR deletion is being classified as intergenic

library(micromicon)

cat("=== Debugging marR Deletion Parsers ===\n\n")

# Test input
annotation <- "coding (152-278 / 435 nt)"
mutation <- "Δ127 bp"
position <- "639,834"

cat("Input:\n")
cat("  annotation:", annotation, "\n")
cat("  mutation:", mutation, "\n")
cat("  position:", position, "\n\n")

# Test annotation geometry parser
cat("Step 1: Parse annotation geometry\n")
ann_geo <- micromicon:::.pm_parse_annotation_geometry(annotation)
cat("  region:", ann_geo$region, "\n")
cat("  coding_pos:", ann_geo$coding_pos, "\n")
cat("  coding_len:", ann_geo$coding_len, "\n")

if (is.na(ann_geo$region)) {
  cat("  ✗ PROBLEM: region is NA!\n\n")
} else if (ann_geo$region == "coding") {
  cat("  ✓ Correctly identified as coding\n\n")
} else {
  cat("  ✗ PROBLEM: region is", ann_geo$region, "not 'coding'\n\n")
}

# Test deletion annotation parser
cat("Step 2: Parse deletion range\n")
del_info <- micromicon:::.pm_parse_deletion_annotation(annotation)
cat("  start:", del_info$start, "\n")
cat("  end:", del_info$end, "\n")
cat("  total_len:", del_info$total_len, "\n")

if (is.na(del_info$start) || is.na(del_info$end)) {
  cat("  ✗ PROBLEM: deletion range not parsed!\n\n")
} else {
  cat("  ✓ Deletion range parsed\n")
  cat("  Deletion size:", del_info$end - del_info$start + 1, "bp\n\n")
}

# Test condition that determines coding vs intergenic
cat("Step 3: Check condition for coding deletion\n")
condition <- !is.na(ann_geo$region) && ann_geo$region == "coding" &&
             !is.na(del_info$start) && !is.na(del_info$end)
cat("  Condition result:", condition, "\n")

if (condition) {
  cat("  ✓ Should be treated as CODING deletion\n")
} else {
  cat("  ✗ PROBLEM: Will be treated as INTERGENIC deletion\n")
  cat("  Breakdown:\n")
  cat("    !is.na(ann_geo$region):", !is.na(ann_geo$region), "\n")
  cat("    ann_geo$region == 'coding':", ann_geo$region == "coding", "\n")
  cat("    !is.na(del_info$start):", !is.na(del_info$start), "\n")
  cat("    !is.na(del_info$end):", !is.na(del_info$end), "\n")
}

cat("\n=== Parser Diagnostic Complete ===\n")

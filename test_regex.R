#!/usr/bin/env Rscript
# Test the regex pattern

annotation <- "coding (152-278 / 435 nt)"

cat("Testing regex pattern:\n")
cat("Input:", annotation, "\n\n")

# Test the new pattern
pattern <- "coding\\s*\\(\\s*([\\d\\-]+)\\s*/\\s*([^\\s)]*?)\\s*nt\\s*\\)"
cat("Pattern:", pattern, "\n\n")

coding_match <- regexec(pattern, annotation)
matches <- regmatches(annotation, coding_match)[[1]]

cat("Matches:\n")
print(matches)
cat("\nLength:", length(matches), "\n")

if (length(matches) >= 2) {
  cat("✓ Pattern matched!\n")
  cat("  Full match:", matches[1], "\n")
  cat("  Position/range:", matches[2], "\n")
  if (length(matches) >= 3) {
    cat("  Length:", matches[3], "\n")
  }

  # Parse position
  pos_str <- matches[2]
  if (grepl("-", pos_str)) {
    coding_pos <- as.integer(strsplit(pos_str, "-")[[1]][1])
    cat("  Parsed position:", coding_pos, "\n")
  }
} else {
  cat("✗ Pattern did NOT match\n")
}

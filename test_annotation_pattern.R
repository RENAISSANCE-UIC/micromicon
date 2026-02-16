#!/usr/bin/env Rscript
# Test different patterns on the failing test case

annotation <- "coding ( 123 / 999 nt )"
cat("Input:", annotation, "\n\n")

# Pattern 1: Fixed spaces (no \s*)
p1 <- "coding \\(([0-9\\-]+) / ([0-9]+) nt\\)"
m1 <- regexec(p1, annotation)
cat("Pattern 1 (fixed spaces):\n")
cat("  Pattern:", p1, "\n")
cat("  Result:", regmatches(annotation, m1)[[1]], "\n\n")

# Pattern 2: Use \\s+ (one or more spaces)
p2 <- "coding\\s+\\(\\s+([0-9\\-]+)\\s+/\\s+([0-9]+)\\s+nt\\s+\\)"
m2 <- regexec(p2, annotation)
cat("Pattern 2 (\\s+ for spaces):\n")
cat("  Pattern:", p2, "\n")
result2 <- regmatches(annotation, m2)[[1]]
if (length(result2) > 0) {
  cat("  ✓ Matched:", result2, "\n\n")
} else {
  cat("  ✗ No match\n\n")
}

# Pattern 3: Use \\s* (zero or more spaces)
p3 <- "coding\\s*\\(\\s*([0-9\\-]+)\\s*/\\s*([0-9]+)\\s*nt\\s*\\)"
m3 <- regexec(p3, annotation)
cat("Pattern 3 (\\s* for optional spaces):\n")
cat("  Pattern:", p3, "\n")
result3 <- regmatches(annotation, m3)[[1]]
if (length(result3) > 0) {
  cat("  ✓ Matched:", result3, "\n\n")
} else {
  cat("  ✗ No match\n\n")
}

# Test the patterns on the user's input too
annotation2 <- "coding (152-278 / 435 nt)"
cat("\nTesting on user's input:", annotation2, "\n\n")

for (i in 1:3) {
  p <- get(paste0("p", i))
  m <- regexec(p, annotation2)
  result <- regmatches(annotation2, m)[[1]]
  cat("Pattern", i, ":")
  if (length(result) > 0) {
    cat(" ✓ Matched:", result[2], "/", result[3], "\n")
  } else {
    cat(" ✗ No match\n")
  }
}

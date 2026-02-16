#!/usr/bin/env Rscript
# Test regex patterns incrementally

annotation <- "coding (152-278 / 435 nt)"
cat("Input:", annotation, "\n\n")

# Test 1: Just "coding"
p1 <- "coding"
m1 <- grepl(p1, annotation)
cat("Pattern 1 'coding':", m1, "\n")

# Test 2: "coding ("
p2 <- "coding \\("
m2 <- grepl(p2, annotation)
cat("Pattern 2 'coding \\(':", m2, "\n")

# Test 3: "coding (152-278"
p3 <- "coding \\(152-278"
m3 <- grepl(p3, annotation)
cat("Pattern 3 'coding \\(152-278':", m3, "\n")

# Test 4: Capture numbers and hyphen
p4 <- "coding \\(([0-9\\-]+)"
m4 <- regexec(p4, annotation)
cat("\nPattern 4 'coding \\(([0-9\\-]+)':\n")
print(regmatches(annotation, m4)[[1]])

# Test 5: Add the / part
p5 <- "coding \\(([0-9\\-]+) / "
m5 <- regexec(p5, annotation)
cat("\nPattern 5 'coding \\(([0-9\\-]+) / ':\n")
print(regmatches(annotation, m5)[[1]])

# Test 6: Add number after /
p6 <- "coding \\(([0-9\\-]+) / ([0-9]+)"
m6 <- regexec(p6, annotation)
cat("\nPattern 6 'coding \\(([0-9\\-]+) / ([0-9]+)':\n")
print(regmatches(annotation, m6)[[1]])

# Test 7: Full pattern
p7 <- "coding \\(([0-9\\-]+) / ([0-9]+) nt\\)"
m7 <- regexec(p7, annotation)
cat("\nPattern 7 'coding \\(([0-9\\-]+) / ([0-9]+) nt\\)':\n")
print(regmatches(annotation, m7)[[1]])

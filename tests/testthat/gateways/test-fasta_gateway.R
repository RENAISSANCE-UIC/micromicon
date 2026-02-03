test_that("fasta_gateway reads FASTA file (simple parser)", {
  # Arrange
  gateway <- create_fasta_gateway(use_bioconductor = FALSE)
  test_file <- "inst/extdata/test.fasta"

  skip_if(!file.exists(test_file), "Test file not found")

  # Act
  sequences <- gateway$read(test_file)

  # Assert
  expect_type(sequences, "character")
  expect_true(length(sequences) == 2)
  expect_equal(names(sequences), c("chr1", "chr2"))
  expect_equal(sequences[["chr1"]], "ATCGATCGATCGATCG")
  expect_equal(sequences[["chr2"]], "GCGCGCGCGCGC")
})

test_that("fasta_gateway writes FASTA file (simple writer)", {
  # Arrange
  gateway <- create_fasta_gateway(use_bioconductor = FALSE)
  sequences <- c(
    test1 = "ATCGATCG",
    test2 = "GCGCGCGC"
  )
  temp_file <- tempfile(fileext = ".fasta")
  on.exit(unlink(temp_file))

  # Act
  gateway$write(sequences, temp_file, wrap_width = 80)

  # Assert
  expect_true(file.exists(temp_file))

  # Read back and verify
  lines <- readLines(temp_file)
  expect_true(any(grepl("^>test1", lines)))
  expect_true(any(grepl("^>test2", lines)))
  expect_true(any(grepl("ATCGATCG", lines)))
  expect_true(any(grepl("GCGCGCGC", lines)))
})

test_that("fasta_gateway wraps long sequences", {
  # Arrange
  gateway <- create_fasta_gateway(use_bioconductor = FALSE)
  long_seq <- paste(rep("ATCG", 50), collapse = "")  # 200 bp
  sequences <- c(long = long_seq)
  temp_file <- tempfile(fileext = ".fasta")
  on.exit(unlink(temp_file))

  # Act
  gateway$write(sequences, temp_file, wrap_width = 80)

  # Assert
  lines <- readLines(temp_file)
  seq_lines <- lines[!grepl("^>", lines)]

  # Each line should be <= 80 characters
  expect_true(all(nchar(seq_lines) <= 80))
})

test_that("fasta_gateway handles missing file gracefully", {
  # Arrange
  gateway <- create_fasta_gateway(use_bioconductor = FALSE)

  # Act & Assert
  expect_error(
    gateway$read("nonexistent.fasta"),
    "not found"
  )
})

test_that("wrap_sequence works correctly", {
  # Test the helper function
  seq <- "ATCGATCGATCGATCG"  # 16 bp

  # Wrap at 8
  wrapped <- wrap_sequence(seq, 8)
  expect_equal(length(wrapped), 2)
  expect_equal(wrapped[1], "ATCGATCG")
  expect_equal(wrapped[2], "ATCGATCG")

  # Wrap at 10
  wrapped <- wrap_sequence(seq, 10)
  expect_equal(length(wrapped), 2)
  expect_equal(wrapped[1], "ATCGATCGAT")
  expect_equal(wrapped[2], "CGATCG")

  # Wrap at 20 (no wrapping needed)
  wrapped <- wrap_sequence(seq, 20)
  expect_equal(length(wrapped), 1)
  expect_equal(wrapped[1], seq)
})

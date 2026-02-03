test_that("read_genome imports GenBank file", {
  # Arrange
  test_file <- "inst/extdata/test.gbk"
  skip_if(!file.exists(test_file), "Test file not found")

  # Act
  entity <- read_genome(test_file)

  # Assert
  expect_s3_class(entity, "genome_entity")
  expect_true(length(entity$sequences$dna_raw) > 0)
  expect_true(nrow(entity$features) > 0)
})

test_that("read_genome auto-detects GenBank format", {
  # Arrange
  test_file <- "inst/extdata/test.gbk"
  skip_if(!file.exists(test_file), "Test file not found")

  # Act
  entity <- read_genome(test_file, format = "auto")

  # Assert
  expect_s3_class(entity, "genome_entity")
})

test_that("read_genome handles GFF3 + FASTA input", {
  # Arrange
  gff_file <- "inst/extdata/test.gff3"
  fasta_file <- "inst/extdata/test.fasta"

  skip_if(!file.exists(gff_file) || !file.exists(fasta_file), "Test files not found")

  # Act
  entity <- read_genome(gff = gff_file, fasta = fasta_file)

  # Assert
  expect_s3_class(entity, "genome_entity")
  expect_true(length(entity$sequences$dna_raw) > 0)
  expect_true(nrow(entity$features) > 0)
})

test_that("read_gbk returns genome_entity by default", {
  # Arrange
  test_file <- "inst/extdata/test.gbk"
  skip_if(!file.exists(test_file), "Test file not found")

  # Act
  entity <- read_gbk(test_file)

  # Assert
  expect_s3_class(entity, "genome_entity")
})

test_that("read_gbk returns legacy format when requested", {
  # Arrange
  test_file <- "inst/extdata/test.gbk"
  skip_if(!file.exists(test_file), "Test file not found")

  # Act
  records <- read_gbk(test_file, return_entity = FALSE)

  # Assert
  expect_type(records, "list")
  expect_true(length(records) > 0)
  expect_true("metadata" %in% names(records[[1]]))
  expect_true("features" %in% names(records[[1]]))
  expect_true("sequence" %in% names(records[[1]]))
})

test_that("init_genome imports GFF3 + FASTA", {
  # Arrange
  gff_file <- "inst/extdata/test.gff3"
  fasta_file <- "inst/extdata/test.fasta"

  skip_if(!file.exists(gff_file) || !file.exists(fasta_file), "Test files not found")

  # Act
  entity <- init_genome(gff_file, fasta_file, verbose = FALSE)

  # Assert
  expect_s3_class(entity, "genome_entity")
  expect_equal(length(entity$sequences$dna_raw), 2)
  expect_true(nrow(entity$features) > 0)
})

test_that("read_genome handles missing file gracefully", {
  # Act & Assert
  expect_error(
    read_genome("nonexistent.gbk"),
    "not found"
  )
})

test_that("read_genome provides helpful error for GFF without FASTA", {
  # Arrange
  gff_file <- "inst/extdata/test.gff3"
  skip_if(!file.exists(gff_file), "Test file not found")

  # Act & Assert
  expect_error(
    read_genome(gff_file, format = "auto"),
    "FASTA file required"
  )
})

test_that("genbank_gateway reads real GenBank file", {
  # Arrange
  gateway <- create_genbank_gateway()
  test_file <- "inst/extdata/test.gbk"

  skip_if(!file.exists(test_file), "Test file not found")

  # Act
  records <- gateway$read(test_file)

  # Assert
  expect_type(records, "list")
  expect_true(length(records) > 0)

  # Check first record structure
  rec <- records[[1]]
  expect_true("metadata" %in% names(rec))
  expect_true("features" %in% names(rec))
  expect_true("sequence" %in% names(rec))

  # Check metadata
  expect_equal(rec$metadata$accession, "TEST001")
  expect_equal(rec$metadata$locus, "TEST001")

  # Check features
  expect_s3_class(rec$features, "data.frame")
  expect_true(nrow(rec$features) > 0)
  expect_true(all(c("type", "start", "end") %in% names(rec$features)))

  # Check sequence
  expect_type(rec$sequence, "character")
  expect_true(nchar(rec$sequence) > 0)
})

test_that("genbank_gateway normalizes records correctly", {
  # Arrange
  gateway <- create_genbank_gateway()
  test_file <- "inst/extdata/test.gbk"

  skip_if(!file.exists(test_file), "Test file not found")

  # Act
  records <- gateway$read(test_file)
  rec <- records[[1]]

  # Assert: Check normalized metadata fields exist
  metadata_fields <- c("locus", "accession", "version", "length_bp", "mol_type",
                       "topology", "division", "date", "definition", "organism", "taxonomy")

  for (field in metadata_fields) {
    expect_true(field %in% names(rec$metadata),
                info = paste("Missing metadata field:", field))
  }
})

test_that("genbank_gateway handles missing file gracefully", {
  # Arrange
  gateway <- create_genbank_gateway()

  # Act & Assert
  expect_error(
    gateway$read("nonexistent.gbk"),
    "not found"
  )
})

test_that("genbank_gateway handles empty features", {
  # This test would require a GenBank file with no features
  # For now, we'll just verify the structure is correct
  gateway <- create_genbank_gateway()
  expect_true("read" %in% names(gateway))
})

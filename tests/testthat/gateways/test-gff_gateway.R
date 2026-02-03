test_that("gff_gateway reads GFF3 file (simple parser)", {
  # Arrange
  gateway <- create_gff_gateway(use_bioconductor = FALSE)
  test_file <- "inst/extdata/test.gff3"

  skip_if(!file.exists(test_file), "Test file not found")

  # Act
  features <- gateway$read(test_file)

  # Assert
  expect_s3_class(features, "data.frame")
  expect_true(nrow(features) == 3)

  # Check required columns
  required_cols <- c("seqname", "start", "end", "strand", "type")
  expect_true(all(required_cols %in% names(features)))

  # Check data
  expect_equal(features$seqname, c("chr1", "chr1", "chr2"))
  expect_equal(features$type, c("gene", "CDS", "gene"))
  expect_equal(features$start, c(1, 1, 1))
  expect_equal(features$end, c(10, 9, 12))
})

test_that("gff_gateway parses attributes correctly", {
  # Arrange
  gateway <- create_gff_gateway(use_bioconductor = FALSE)
  test_file <- "inst/extdata/test.gff3"

  skip_if(!file.exists(test_file), "Test file not found")

  # Act
  features <- gateway$read(test_file)

  # Assert: Check parsed attributes
  expect_true("gene" %in% names(features))
  expect_true("locus_tag" %in% names(features))

  # Check specific values
  expect_equal(features$gene[1], "testA")
  expect_equal(features$locus_tag[1], "GENE001")
  expect_equal(features$gene[3], "testB")
  expect_equal(features$locus_tag[3], "GENE002")
})

test_that("gff_gateway handles missing file gracefully", {
  # Arrange
  gateway <- create_gff_gateway(use_bioconductor = FALSE)

  # Act & Assert
  expect_error(
    gateway$read("nonexistent.gff3"),
    "not found"
  )
})

test_that("parse_gff_attributes extracts common fields", {
  # Arrange
  df <- data.frame(
    seqname = "chr1",
    start = 1,
    end = 10,
    strand = "+",
    type = "gene",
    attributes = "ID=gene1;Name=testA;gene=testA;locus_tag=GENE001;product=test protein",
    stringsAsFactors = FALSE
  )

  # Act
  result <- parse_gff_attributes(df)

  # Assert
  expect_equal(result$gene, "testA")
  expect_equal(result$locus_tag, "GENE001")
  expect_equal(result$product, "test protein")
  expect_equal(result$ID, "gene1")
  expect_equal(result$Name, "testA")
})

test_that("parse_gff_attributes handles URL encoding", {
  # Arrange
  df <- data.frame(
    seqname = "chr1",
    start = 1,
    end = 10,
    strand = "+",
    type = "gene",
    attributes = "product=test%20protein%20A",  # Space encoded as %20
    stringsAsFactors = FALSE
  )

  # Act
  result <- parse_gff_attributes(df)

  # Assert
  expect_equal(result$product, "test protein A")
})

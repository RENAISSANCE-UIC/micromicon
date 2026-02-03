test_that("import_genbank transforms data correctly", {
  # Arrange: Create mock gateway
  mock_gateway <- list(
    read = function(path) {
      list(
        list(
          metadata = list(
            locus = "TESTLOC",
            accession = "TEST001",
            length_bp = 10,
            mol_type = "DNA",
            topology = "linear",
            division = "BCT",
            date = "01-JAN-2020",
            definition = "Test record",
            version = "TEST001.1",
            organism = "Test organism",
            taxonomy = "Bacteria"
          ),
          features = data.frame(
            type = c("gene", "CDS"),
            start = c(1, 1),
            end = c(10, 9),
            strand = c("+", "+"),
            gene = c("testA", "testA"),
            locus_tag = c("GENE001", "GENE001"),
            product = c(NA, "Test protein"),
            stringsAsFactors = FALSE
          ),
          sequence = "ATCGATCGAT"
        )
      )
    }
  )

  # Act
  entity <- execute_import_genbank(mock_gateway, "fake.gbk")

  # Assert
  expect_s3_class(entity, "genome_entity")
  expect_equal(length(entity$sequences$dna_raw), 1)
  expect_equal(entity$sequences$dna_raw[["TEST001"]], "ATCGATCGAT")
  expect_equal(nrow(entity$features), 2)
  expect_equal(entity$features$seqname[1], "TEST001")
  expect_equal(entity$features$type[1], "gene")
  expect_equal(nrow(entity$metadata), 1)
  expect_equal(entity$metadata$seqname[1], "TEST001")
  expect_equal(entity$metadata$topology[1], "linear")
})

test_that("import_genbank handles multiple records", {
  # Arrange
  mock_gateway <- list(
    read = function(path) {
      list(
        list(
          metadata = list(locus = "REC1", accession = "ACC1", length_bp = 8),
          features = data.frame(
            type = "gene",
            start = 1,
            end = 8,
            strand = "+",
            gene = "geneA",
            stringsAsFactors = FALSE
          ),
          sequence = "ATCGATCG"
        ),
        list(
          metadata = list(locus = "REC2", accession = "ACC2", length_bp = 6),
          features = data.frame(
            type = "gene",
            start = 1,
            end = 6,
            strand = "+",
            gene = "geneB",
            stringsAsFactors = FALSE
          ),
          sequence = "GCGCGC"
        )
      )
    }
  )

  # Act
  entity <- execute_import_genbank(mock_gateway, "multi.gbk")

  # Assert
  expect_equal(length(entity$sequences$dna_raw), 2)
  expect_equal(names(entity$sequences$dna_raw), c("ACC1", "ACC2"))
  expect_equal(nrow(entity$features), 2)
  expect_equal(nrow(entity$metadata), 2)
  expect_true("ACC1" %in% entity$indices$seqnames)
  expect_true("ACC2" %in% entity$indices$seqnames)
})

test_that("import_genbank builds indices correctly", {
  # Arrange
  mock_gateway <- list(
    read = function(path) {
      list(
        list(
          metadata = list(locus = "TEST", accession = "TEST001"),
          features = data.frame(
            type = c("gene", "CDS"),
            start = c(1, 1),
            end = c(10, 9),
            strand = c("+", "+"),
            gene = c("geneA", "geneA"),
            locus_tag = c("GENE001", "GENE001"),
            stringsAsFactors = FALSE
          ),
          sequence = "ATCGATCGAT"
        )
      )
    }
  )

  # Act
  entity <- execute_import_genbank(mock_gateway, "test.gbk")

  # Assert: Check indices
  expect_true("GENE001" %in% names(entity$indices$locus_tag_index))
  expect_equal(entity$indices$locus_tag_index[["GENE001"]], 1)  # First row

  expect_true("geneA" %in% names(entity$indices$gene_index))
  expect_equal(length(entity$indices$gene_index[["geneA"]]), 2)  # Both features
})

test_that("import_genbank validates entity after creation", {
  # Arrange: Mock gateway returns invalid data (coordinates out of bounds)
  bad_gateway <- list(
    read = function(path) {
      list(
        list(
          metadata = list(locus = "TEST", accession = "TEST001"),
          features = data.frame(
            type = "gene",
            start = 1,
            end = 100,  # Beyond sequence length
            strand = "+",
            stringsAsFactors = FALSE
          ),
          sequence = "ATCG"  # Only 4 bp
        )
      )
    }
  )

  # Act & Assert
  expect_error(
    execute_import_genbank(bad_gateway, "bad.gbk"),
    "exceed"
  )
})

test_that("import_genbank handles empty features", {
  # Arrange
  mock_gateway <- list(
    read = function(path) {
      list(
        list(
          metadata = list(locus = "TEST", accession = "TEST001"),
          features = data.frame(),  # No features
          sequence = "ATCGATCG"
        )
      )
    }
  )

  # Act
  entity <- execute_import_genbank(mock_gateway, "no_features.gbk")

  # Assert
  expect_equal(nrow(entity$features), 0)
  expect_equal(length(entity$sequences$dna_raw), 1)
})

test_that("import_genbank uses locus when accession missing", {
  # Arrange
  mock_gateway <- list(
    read = function(path) {
      list(
        list(
          metadata = list(
            locus = "MYLOCUS",
            accession = NA,  # No accession
            length_bp = 8
          ),
          features = data.frame(),
          sequence = "ATCGATCG"
        )
      )
    }
  )

  # Act
  entity <- execute_import_genbank(mock_gateway, "test.gbk")

  # Assert
  expect_equal(names(entity$sequences$dna_raw)[1], "MYLOCUS")
  expect_equal(entity$metadata$seqname[1], "MYLOCUS")
})

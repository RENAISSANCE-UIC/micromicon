test_that("import_gff_fasta combines GFF and FASTA correctly", {
  # Arrange: Mock gateways
  mock_gff <- list(
    read = function(path) {
      data.frame(
        seqname = c("chr1", "chr1"),
        start = c(1, 5),
        end = c(4, 8),
        strand = c("+", "+"),
        type = c("gene", "CDS"),
        gene = c("geneA", "geneA"),
        locus_tag = c("GENE001", "GENE001"),
        stringsAsFactors = FALSE
      )
    }
  )

  mock_fasta <- list(
    read = function(path) {
      c(chr1 = "ATCGATCG")
    }
  )

  # Act
  entity <- execute_import_gff_fasta(
    mock_gff, mock_fasta,
    "test.gff3", "test.fasta",
    options = list(verbose = FALSE)
  )

  # Assert
  expect_s3_class(entity, "genome_entity")
  expect_equal(length(entity$sequences$dna_raw), 1)
  expect_equal(entity$sequences$dna_raw[["chr1"]], "ATCGATCG")
  expect_equal(nrow(entity$features), 2)
  expect_equal(entity$features$seqname[1], "chr1")
  expect_equal(nrow(entity$metadata), 1)
  expect_equal(entity$metadata$length_bp[1], 8)
})

test_that("import_gff_fasta harmonizes seqnames", {
  # Arrange: GFF has "1", FASTA has "chr1"
  mock_gff <- list(
    read = function(path) {
      data.frame(
        seqname = "1",  # Numeric seqname
        start = 1,
        end = 10,
        strand = "+",
        type = "gene",
        stringsAsFactors = FALSE
      )
    }
  )

  mock_fasta <- list(
    read = function(path) {
      c(chr1 = "ATCGATCGATCG")  # With "chr" prefix
    }
  )

  # Act
  entity <- execute_import_gff_fasta(
    mock_gff, mock_fasta,
    "test.gff3", "test.fasta",
    options = list(auto_harmonize = TRUE, verbose = FALSE)
  )

  # Assert: Should harmonize "1" -> "chr1"
  expect_equal(nrow(entity$features), 1)
  expect_equal(entity$features$seqname[1], "chr1")
})

test_that("import_gff_fasta filters orphaned features", {
  # Arrange: Feature on seqname not in FASTA
  mock_gff <- list(
    read = function(path) {
      data.frame(
        seqname = c("chr1", "chr2", "chr3"),  # chr3 not in FASTA
        start = c(1, 1, 1),
        end = c(10, 10, 10),
        strand = c("+", "+", "+"),
        type = c("gene", "gene", "gene"),
        stringsAsFactors = FALSE
      )
    }
  )

  mock_fasta <- list(
    read = function(path) {
      c(chr1 = "ATCGATCGATCG", chr2 = "GCGCGCGCGCGC")  # No chr3
    }
  )

  # Act
  expect_warning(
    entity <- execute_import_gff_fasta(
      mock_gff, mock_fasta,
      "test.gff3", "test.fasta",
      options = list(auto_harmonize = FALSE, verbose = FALSE)
    ),
    "mismatch"
  )

  # Assert: chr3 feature should be filtered out
  expect_equal(nrow(entity$features), 2)
  expect_false("chr3" %in% entity$features$seqname)
})

test_that("import_gff_fasta builds indices", {
  # Arrange
  mock_gff <- list(
    read = function(path) {
      data.frame(
        seqname = "chr1",
        start = 1,
        end = 10,
        strand = "+",
        type = "gene",
        gene = "geneA",
        locus_tag = "GENE001",
        stringsAsFactors = FALSE
      )
    }
  )

  mock_fasta <- list(
    read = function(path) {
      c(chr1 = "ATCGATCGATCG")
    }
  )

  # Act
  entity <- execute_import_gff_fasta(
    mock_gff, mock_fasta,
    "test.gff3", "test.fasta",
    options = list(verbose = FALSE)
  )

  # Assert
  expect_true("GENE001" %in% names(entity$indices$locus_tag_index))
  expect_true("geneA" %in% names(entity$indices$gene_index))
})

test_that("import_gff_fasta validates entity", {
  # Arrange: Invalid coordinates
  mock_gff <- list(
    read = function(path) {
      data.frame(
        seqname = "chr1",
        start = 1,
        end = 100,  # Beyond sequence
        strand = "+",
        type = "gene",
        stringsAsFactors = FALSE
      )
    }
  )

  mock_fasta <- list(
    read = function(path) {
      c(chr1 = "ATCG")  # Only 4 bp
    }
  )

  # Act & Assert
  expect_error(
    execute_import_gff_fasta(
      mock_gff, mock_fasta,
      "bad.gff3", "bad.fasta",
      options = list(verbose = FALSE)
    ),
    "exceed"
  )
})

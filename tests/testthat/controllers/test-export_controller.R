test_that("write_fasta writes genome_entity to file", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = c(test1 = "ATCGATCG", test2 = "GCGCGCGC"),
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(),
    metadata_df = data.frame(seqname = c("test1", "test2"), stringsAsFactors = FALSE),
    indices_list = list(seqnames = c("test1", "test2"))
  )

  temp_file <- tempfile(fileext = ".fasta")
  on.exit(unlink(temp_file))

  # Act
  result <- write_fasta(entity, temp_file)

  # Assert
  expect_true(file.exists(temp_file))
  expect_equal(result, temp_file)

  # Read back and verify
  lines <- readLines(temp_file)
  expect_true(any(grepl("^>test1", lines)))
  expect_true(any(grepl("ATCGATCG", lines)))
})

test_that("write_fasta writes character vector to file", {
  # Arrange
  seqs <- c(seq1 = "ATCG", seq2 = "GCGC")
  temp_file <- tempfile(fileext = ".fasta")
  on.exit(unlink(temp_file))

  # Act
  write_fasta(seqs, temp_file)

  # Assert
  expect_true(file.exists(temp_file))

  lines <- readLines(temp_file)
  expect_true(any(grepl("^>seq1", lines)))
  expect_true(any(grepl("ATCG", lines)))
})

test_that("write_gff3 writes features to GFF3 file", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(dna_raw = c(chr1 = "ATCGATCG"), dna_bio = NULL, indexed_fa = NULL),
    features_df = data.frame(
      seqname = "chr1",
      start = 1,
      end = 8,
      strand = "+",
      type = "gene",
      gene = "testA",
      locus_tag = "GENE001",
      stringsAsFactors = FALSE
    ),
    metadata_df = data.frame(seqname = "chr1", stringsAsFactors = FALSE),
    indices_list = list(seqnames = "chr1")
  )

  temp_file <- tempfile(fileext = ".gff3")
  on.exit(unlink(temp_file))

  # Act
  result <- write_gff3(entity, temp_file)

  # Assert
  expect_true(file.exists(temp_file))
  expect_equal(result, temp_file)

  # Read back and verify
  lines <- readLines(temp_file)
  expect_true(any(grepl("##gff-version 3", lines)))
  expect_true(any(grepl("chr1", lines)))
  expect_true(any(grepl("gene=testA", lines)))
  expect_true(any(grepl("locus_tag=GENE001", lines)))
})

test_that("write_gbk_fasta is a wrapper for write_fasta", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = c(test = "ATCG"),
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(),
    metadata_df = data.frame(seqname = "test", stringsAsFactors = FALSE),
    indices_list = list(seqnames = "test")
  )

  temp_file <- tempfile(fileext = ".fasta")
  on.exit(unlink(temp_file))

  # Act
  write_gbk_fasta(entity, temp_file)

  # Assert
  expect_true(file.exists(temp_file))
})

test_that("write_gbk_gff3 is a wrapper for write_gff3", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(dna_raw = c(chr1 = "ATCG"), dna_bio = NULL, indexed_fa = NULL),
    features_df = data.frame(
      seqname = "chr1",
      start = 1,
      end = 4,
      strand = "+",
      type = "gene",
      stringsAsFactors = FALSE
    ),
    metadata_df = data.frame(seqname = "chr1", stringsAsFactors = FALSE),
    indices_list = list(seqnames = "chr1")
  )

  temp_file <- tempfile(fileext = ".gff3")
  on.exit(unlink(temp_file))

  # Act
  write_gbk_gff3(entity, temp_file)

  # Assert
  expect_true(file.exists(temp_file))
})

test_that("write_gff3 handles empty features", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(dna_raw = c(chr1 = "ATCG"), dna_bio = NULL, indexed_fa = NULL),
    features_df = data.frame(),  # No features
    metadata_df = data.frame(seqname = "chr1", stringsAsFactors = FALSE),
    indices_list = list(seqnames = "chr1")
  )

  temp_file <- tempfile(fileext = ".gff3")
  on.exit(unlink(temp_file))

  # Act & Assert
  expect_warning(
    write_gff3(entity, temp_file),
    "No features"
  )
})

test_that("build_gff3_attributes creates valid attributes string", {
  # Arrange
  feat <- data.frame(
    seqname = "chr1",
    start = 1,
    end = 10,
    type = "gene",
    gene = "testA",
    locus_tag = "GENE001",
    product = "test protein",
    stringsAsFactors = FALSE
  )

  # Act
  attrs <- build_gff3_attributes(feat)

  # Assert
  expect_type(attrs, "character")
  expect_true(grepl("gene=testA", attrs))
  expect_true(grepl("locus_tag=GENE001", attrs))
  expect_true(grepl("product=test", attrs))
})

test_that("build_gff3_attributes handles missing fields", {
  # Arrange
  feat <- data.frame(
    seqname = "chr1",
    start = 1,
    end = 10,
    type = "gene",
    stringsAsFactors = FALSE
  )

  # Act
  attrs <- build_gff3_attributes(feat)

  # Assert: Should return "." when no attributes
  expect_equal(attrs, ".")
})

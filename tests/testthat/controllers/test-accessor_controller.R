test_that("sequences returns character vector", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = c(chr1 = "ATCG", chr2 = "GCGC"),
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(),
    metadata_df = data.frame(seqname = c("chr1", "chr2"), stringsAsFactors = FALSE),
    indices_list = list(seqnames = c("chr1", "chr2"))
  )

  # Act
  seqs <- sequences(entity)

  # Assert
  expect_type(seqs, "character")
  expect_equal(length(seqs), 2)
  expect_equal(names(seqs), c("chr1", "chr2"))
  expect_equal(seqs[["chr1"]], "ATCG")
})

test_that("features returns data.frame", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(dna_raw = c(chr1 = "ATCG"), dna_bio = NULL, indexed_fa = NULL),
    features_df = data.frame(
      seqname = "chr1",
      start = 1,
      end = 4,
      type = "gene",
      gene = "testA",
      stringsAsFactors = FALSE
    ),
    metadata_df = data.frame(seqname = "chr1", stringsAsFactors = FALSE),
    indices_list = list(seqnames = "chr1")
  )

  # Act
  feats <- features(entity)

  # Assert
  expect_s3_class(feats, "data.frame")
  expect_equal(nrow(feats), 1)
  expect_equal(feats$gene[1], "testA")
})

test_that("features filters by type", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(dna_raw = c(chr1 = "ATCG"), dna_bio = NULL, indexed_fa = NULL),
    features_df = data.frame(
      seqname = c("chr1", "chr1"),
      start = c(1, 1),
      end = c(4, 3),
      type = c("gene", "CDS"),
      stringsAsFactors = FALSE
    ),
    metadata_df = data.frame(seqname = "chr1", stringsAsFactors = FALSE),
    indices_list = list(seqnames = "chr1")
  )

  # Act
  genes <- features(entity, type = "gene")

  # Assert
  expect_equal(nrow(genes), 1)
  expect_equal(genes$type[1], "gene")
})

test_that("metadata returns data.frame", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(dna_raw = c(chr1 = "ATCG"), dna_bio = NULL, indexed_fa = NULL),
    features_df = data.frame(),
    metadata_df = data.frame(
      seqname = "chr1",
      length = 4,
      topology = "linear",
      stringsAsFactors = FALSE
    ),
    indices_list = list(seqnames = "chr1")
  )

  # Act
  meta <- metadata(entity)

  # Assert
  expect_s3_class(meta, "data.frame")
  expect_equal(nrow(meta), 1)
  expect_equal(meta$seqname[1], "chr1")
  expect_equal(meta$topology[1], "linear")
})

test_that("seqnames returns character vector", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = c(chr1 = "ATCG", chr2 = "GCGC"),
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(),
    metadata_df = data.frame(seqname = c("chr1", "chr2"), stringsAsFactors = FALSE),
    indices_list = list(seqnames = c("chr1", "chr2"))
  )

  # Act
  names <- seqnames(entity)

  # Assert
  expect_type(names, "character")
  expect_equal(length(names), 2)
  expect_equal(names, c("chr1", "chr2"))
})

test_that("accessor functions validate input", {
  # Arrange
  not_entity <- list(foo = "bar")

  # Act & Assert
  expect_error(sequences(not_entity), "genome_entity")
  expect_error(features(not_entity), "genome_entity")
  expect_error(metadata(not_entity), "genome_entity")
  expect_error(seqnames(not_entity), "genome_entity")
})

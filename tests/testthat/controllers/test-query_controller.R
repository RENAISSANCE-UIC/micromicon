test_that("extract_sequences_by_name finds features", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = c(chr1 = "ATCGATCGATCGATCG"),
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(
      seqname = "chr1",
      start = 1,
      end = 4,
      strand = "+",
      type = "gene",
      gene = "dnaA",
      locus_tag = "GENE001",
      stringsAsFactors = FALSE
    ),
    metadata_df = data.frame(seqname = "chr1", stringsAsFactors = FALSE),
    indices_list = list(seqnames = "chr1")
  )

  # Act
  seqs <- extract_sequences_by_name(entity, "dnaA")

  # Assert
  expect_equal(length(seqs), 1)
  expect_equal(as.character(seqs[1]), "ATCG")
  expect_true(grepl("dnaA", names(seqs)[1]))
})

test_that("extract_sequences_by_coords extracts region", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = c(chr1 = "ATCGATCGATCGATCG"),
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(),
    metadata_df = data.frame(seqname = "chr1", stringsAsFactors = FALSE),
    indices_list = list(seqnames = "chr1")
  )

  # Act
  seq <- extract_sequences_by_coords(entity, "chr1", 1, 4)

  # Assert
  expect_equal(length(seq), 1)
  expect_equal(as.character(seq), "ATCG")
  expect_equal(names(seq), "chr1:1-4")
})

test_that("search_features filters by type", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(dna_raw = c(chr1 = "ATCGATCGATCG"), dna_bio = NULL, indexed_fa = NULL),
    features_df = data.frame(
      seqname = c("chr1", "chr1", "chr1"),
      start = c(1, 1, 5),
      end = c(4, 3, 8),
      type = c("gene", "CDS", "gene"),
      gene = c("testA", "testA", "testB"),
      strand = c("+", "+", "+"),
      stringsAsFactors = FALSE
    ),
    metadata_df = data.frame(seqname = "chr1", stringsAsFactors = FALSE),
    indices_list = list(seqnames = "chr1")
  )

  # Act
  genes <- search_features(entity, type = "gene")

  # Assert
  expect_equal(nrow(genes), 2)
  expect_true(all(genes$type == "gene"))
})

test_that("search_features filters by pattern", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(dna_raw = c(chr1 = "ATCGATCGATCG"), dna_bio = NULL, indexed_fa = NULL),
    features_df = data.frame(
      seqname = c("chr1", "chr1"),
      start = c(1, 5),
      end = c(4, 8),
      type = c("gene", "gene"),
      gene = c("dnaA", "gyrA"),
      strand = c("+", "+"),
      stringsAsFactors = FALSE
    ),
    metadata_df = data.frame(seqname = "chr1", stringsAsFactors = FALSE),
    indices_list = list(seqnames = "chr1")
  )

  # Act
  dna_genes <- search_features(entity, pattern = "dna")

  # Assert
  expect_equal(nrow(dna_genes), 1)
  expect_equal(dna_genes$gene[1], "dnaA")
})

test_that("search_features filters by seqname", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = c(chr1 = "ATCG", chr2 = "GCGC"),
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(
      seqname = c("chr1", "chr2"),
      start = c(1, 1),
      end = c(4, 4),
      type = c("gene", "gene"),
      strand = c("+", "+"),
      stringsAsFactors = FALSE
    ),
    metadata_df = data.frame(seqname = c("chr1", "chr2"), stringsAsFactors = FALSE),
    indices_list = list(seqnames = c("chr1", "chr2"))
  )

  # Act
  chr1_genes <- search_features(entity, seqname = "chr1")

  # Assert
  expect_equal(nrow(chr1_genes), 1)
  expect_equal(chr1_genes$seqname[1], "chr1")
})

test_that("search_features filters by coordinates", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(dna_raw = c(chr1 = "ATCGATCGATCG"), dna_bio = NULL, indexed_fa = NULL),
    features_df = data.frame(
      seqname = "chr1",
      start = c(1, 5, 9),
      end = c(4, 8, 12),
      type = c("gene", "gene", "gene"),
      strand = c("+", "+", "+"),
      stringsAsFactors = FALSE
    ),
    metadata_df = data.frame(seqname = "chr1", stringsAsFactors = FALSE),
    indices_list = list(seqnames = "chr1")
  )

  # Act: Features starting at position 5 or later
  late_features <- search_features(entity, start = 5)

  # Assert
  expect_equal(nrow(late_features), 2)
  expect_true(all(late_features$start >= 5))
})

test_that("query functions validate input", {
  # Arrange
  not_entity <- list(foo = "bar")

  # Act & Assert
  expect_error(extract_sequences_by_name(not_entity, "test"), "genome_entity")
  expect_error(extract_sequences_by_coords(not_entity, "chr1", 1, 10), "genome_entity")
  expect_error(search_features(not_entity, type = "gene"), "genome_entity")
})

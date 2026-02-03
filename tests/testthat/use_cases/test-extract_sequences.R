test_that("extract_sequences_by_name finds features by gene name", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = c(chr1 = "ATCGATCGATCGATCG"),
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(
      seqname = c("chr1", "chr1"),
      start = c(1, 5),
      end = c(4, 8),
      strand = c("+", "+"),
      type = c("gene", "CDS"),
      gene = c("dnaA", "dnaA"),
      locus_tag = c("GENE001", "GENE001"),
      stringsAsFactors = FALSE
    ),
    metadata_df = data.frame(
      seqname = "chr1",
      length_bp = 16,
      stringsAsFactors = FALSE
    ),
    indices_list = list(
      seqnames = "chr1"
    )
  )

  # Act
  seqs <- execute_extract_sequences_by_name(entity, "dnaA")

  # Assert
  expect_equal(length(seqs), 2)
  expect_equal(seqs[1], c("chr1:1-4 dnaA" = "ATCG"))
  expect_equal(seqs[2], c("chr1:5-8 dnaA" = "ATCG"))
})

test_that("extract_sequences_by_name filters by type", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = c(chr1 = "ATCGATCGATCGATCG"),
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(
      seqname = c("chr1", "chr1"),
      start = c(1, 5),
      end = c(4, 8),
      strand = c("+", "+"),
      type = c("gene", "CDS"),
      gene = c("dnaA", "dnaA"),
      stringsAsFactors = FALSE
    ),
    metadata_df = data.frame(
      seqname = "chr1",
      length_bp = 16,
      stringsAsFactors = FALSE
    ),
    indices_list = list(seqnames = "chr1")
  )

  # Act: Filter to only CDS
  seqs <- execute_extract_sequences_by_name(
    entity, "dnaA",
    options = list(type_filter = "CDS")
  )

  # Assert: Only CDS feature
  expect_equal(length(seqs), 1)
  expect_true(grepl("5-8", names(seqs)[1]))
})

test_that("extract_sequences_by_name handles minus strand", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = c(chr1 = "ATCGATCG"),
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(
      seqname = "chr1",
      start = 1,
      end = 4,
      strand = "-",  # Minus strand
      type = "gene",
      gene = "testGene",
      stringsAsFactors = FALSE
    ),
    metadata_df = data.frame(
      seqname = "chr1",
      length_bp = 8,
      stringsAsFactors = FALSE
    ),
    indices_list = list(seqnames = "chr1")
  )

  # Act
  seqs <- execute_extract_sequences_by_name(entity, "testGene")

  # Assert: Should be reverse complement of ATCG = CGAT
  expect_equal(as.character(seqs[1]), "CGAT")
})

test_that("extract_sequences_by_coords extracts single region", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = c(chr1 = "ATCGATCGATCGATCG"),
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(),
    metadata_df = data.frame(
      seqname = "chr1",
      length_bp = 16,
      stringsAsFactors = FALSE
    ),
    indices_list = list(seqnames = "chr1")
  )

  # Act
  seq <- execute_extract_sequences_by_coords(entity, "chr1", 1, 4)

  # Assert
  expect_equal(length(seq), 1)
  expect_equal(as.character(seq), "ATCG")
  expect_equal(names(seq), "chr1:1-4")
})

test_that("extract_sequences_by_coords extracts multiple regions", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = c(chr1 = "ATCGATCGATCGATCG", chr2 = "GCGCGCGC"),
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(),
    metadata_df = data.frame(
      seqname = c("chr1", "chr2"),
      length_bp = c(16, 8),
      stringsAsFactors = FALSE
    ),
    indices_list = list(seqnames = c("chr1", "chr2"))
  )

  # Act
  seqs <- execute_extract_sequences_by_coords(
    entity,
    seqname = c("chr1", "chr2"),
    start = c(1, 1),
    end = c(4, 4)
  )

  # Assert
  expect_equal(length(seqs), 2)
  expect_equal(as.character(seqs[1]), "ATCG")
  expect_equal(as.character(seqs[2]), "GCGC")
})

test_that("extract_sequences_by_coords handles reverse complement", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = c(chr1 = "ATCGATCG"),
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(),
    metadata_df = data.frame(
      seqname = "chr1",
      length_bp = 8,
      stringsAsFactors = FALSE
    ),
    indices_list = list(seqnames = "chr1")
  )

  # Act
  seq <- execute_extract_sequences_by_coords(
    entity, "chr1", 1, 4,
    options = list(strand = "-")
  )

  # Assert: Reverse complement of ATCG = CGAT
  expect_equal(as.character(seq), "CGAT")
})

test_that("extract_sequences_by_coords validates coordinates", {
  # Arrange
  entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = c(chr1 = "ATCG"),
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(),
    metadata_df = data.frame(
      seqname = "chr1",
      length_bp = 4,
      stringsAsFactors = FALSE
    ),
    indices_list = list(seqnames = "chr1")
  )

  # Act & Assert: Out of bounds
  expect_error(
    execute_extract_sequences_by_coords(entity, "chr1", 1, 100),
    "exceed"
  )

  # Unknown sequence
  expect_error(
    execute_extract_sequences_by_coords(entity, "chr2", 1, 4),
    "not found"
  )
})

test_that("reverse_complement works correctly", {
  expect_equal(reverse_complement("ATCG"), "CGAT")
  expect_equal(reverse_complement("AAAA"), "TTTT")
  expect_equal(reverse_complement("GCGC"), "GCGC")  # Palindrome
  expect_equal(reverse_complement("atcg"), "cgat")  # Lowercase
})

test_that("translate_dna works correctly", {
  # ATG = M (Met, start codon)
  expect_equal(translate_dna("ATG"), "M")

  # ATG GCT TAA = M A * (stop)
  expect_equal(translate_dna("ATGGCTTAA"), "MA*")

  # Incomplete codon (only 2 bases) - should not translate last 2
  expect_equal(translate_dna("ATGGC"), "M")

  # Empty or too short
  expect_equal(translate_dna("AT"), "")
  expect_equal(translate_dna(""), "")
})

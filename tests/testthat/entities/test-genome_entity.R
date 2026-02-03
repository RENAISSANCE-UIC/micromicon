test_that("new_genome_entity creates valid entity", {
  # Arrange
  sequences_list <- list(
    dna_raw = c(chr1 = "ATCGATCG", chr2 = "GCGCGC"),
    dna_bio = NULL,
    indexed_fa = NULL
  )

  features_df <- data.frame(
    seqname = c("chr1", "chr1", "chr2"),
    start = c(1, 5, 1),
    end = c(4, 8, 6),
    strand = c("+", "+", "+"),
    type = c("gene", "CDS", "gene"),
    locus_tag = c("GENE001", "GENE001", "GENE002"),
    stringsAsFactors = FALSE
  )

  metadata_df <- data.frame(
    seqname = c("chr1", "chr2"),
    length = c(8, 6),
    topology = c("linear", "circular"),
    molecule_type = c("DNA", "DNA"),
    stringsAsFactors = FALSE
  )

  indices_list <- list(
    seqnames = c("chr1", "chr2"),
    locus_tag_index = c(GENE001 = 1, GENE002 = 3),
    gene_index = list()
  )

  # Act
  entity <- new_genome_entity(sequences_list, features_df, metadata_df, indices_list)

  # Assert
  expect_s3_class(entity, "genome_entity")
  expect_true(all(c("sequences", "features", "metadata", "indices") %in% names(entity)))
  expect_equal(length(entity$sequences$dna_raw), 2)
  expect_equal(nrow(entity$features), 3)
})

test_that("validate_genome_entity accepts valid entity", {
  # Arrange
  valid_entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = c(chr1 = "ATCGATCG"),
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(
      seqname = "chr1",
      start = 1,
      end = 4,
      strand = "+",
      type = "gene",
      stringsAsFactors = FALSE
    ),
    metadata_df = data.frame(
      seqname = "chr1",
      length = 8,
      stringsAsFactors = FALSE
    ),
    indices_list = list(
      seqnames = "chr1",
      locus_tag_index = integer(),
      gene_index = list()
    )
  )

  # Act & Assert
  expect_true(validate_genome_entity(valid_entity))
})

test_that("validate_genome_entity rejects non-genome_entity object", {
  # Arrange
  not_entity <- list(foo = "bar")

  # Act & Assert
  expect_error(validate_genome_entity(not_entity), "not a genome_entity")
})

test_that("validate_genome_entity rejects missing components", {
  # Arrange
  incomplete_entity <- structure(
    list(sequences = list(), features = data.frame()),  # Missing metadata and indices
    class = "genome_entity"
  )

  # Act & Assert
  expect_error(validate_genome_entity(incomplete_entity), "Missing required components")
})

test_that("validate_genome_entity rejects invalid coordinates", {
  # Arrange
  invalid_entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = c(chr1 = "ATCG"),  # Only 4 bp
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(
      seqname = "chr1",
      start = 1,
      end = 10,  # Beyond sequence length
      strand = "+",
      stringsAsFactors = FALSE
    ),
    metadata_df = data.frame(
      seqname = "chr1",
      length = 4,
      stringsAsFactors = FALSE
    ),
    indices_list = list(
      seqnames = "chr1"
    )
  )

  # Act & Assert
  expect_error(validate_genome_entity(invalid_entity), "exceed")
})

test_that("validate_genome_entity rejects invalid strand", {
  # Arrange
  invalid_entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = c(chr1 = "ATCG"),
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(
      seqname = "chr1",
      start = 1,
      end = 4,
      strand = "INVALID",  # Invalid strand
      stringsAsFactors = FALSE
    ),
    metadata_df = data.frame(
      seqname = "chr1",
      length = 4,
      stringsAsFactors = FALSE
    ),
    indices_list = list(
      seqnames = "chr1"
    )
  )

  # Act & Assert
  expect_error(validate_genome_entity(invalid_entity), "Invalid strand")
})

test_that("validate_genome_entity rejects feature referencing unknown sequence", {
  # Arrange
  invalid_entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = c(chr1 = "ATCG"),
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(
      seqname = "chr2",  # Sequence doesn't exist
      start = 1,
      end = 4,
      strand = "+",
      stringsAsFactors = FALSE
    ),
    metadata_df = data.frame(
      seqname = "chr1",
      length = 4,
      stringsAsFactors = FALSE
    ),
    indices_list = list(
      seqnames = "chr1"
    )
  )

  # Act & Assert
  expect_error(validate_genome_entity(invalid_entity), "unknown sequence")
})

test_that("genome_entity handles empty sequences and features", {
  # Arrange
  empty_entity <- new_genome_entity(
    sequences_list = list(
      dna_raw = character(),
      dna_bio = NULL,
      indexed_fa = NULL
    ),
    features_df = data.frame(),
    metadata_df = data.frame(),
    indices_list = list(
      seqnames = character()
    )
  )

  # Act & Assert
  expect_true(validate_genome_entity(empty_entity))
})

test_that("has_bioconductor works", {
  result <- has_bioconductor()
  expect_type(result, "logical")
  expect_length(result, 1)
})

test_that("sequences accessor works with character format", {
  # Create test entity with sequences
  dna <- c(seq1 = "ATCGATCG", seq2 = "GGCCTTAA")
  entity <- new_genome_entity(
    sequences_list = list(dna_raw = dna, dna_bio = NULL, indexed_fa = NULL)
  )

  seqs <- sequences(entity, format = "character")
  expect_type(seqs, "character")
  expect_equal(names(seqs), c("seq1", "seq2"))
  expect_equal(seqs["seq1"], c(seq1 = "ATCGATCG"))
})

test_that("sequences accessor errors without Bioconductor", {
  entity <- new_genome_entity()

  if (!has_bioconductor()) {
    expect_error(
      sequences(entity, format = "DNAStringSet"),
      "DNAStringSet format requires Bioconductor"
    )
  }
})

test_that("features accessor works with data.frame format", {
  # Create test entity with features
  feat_df <- data.frame(
    seqname = "seq1",
    start = 100,
    end = 200,
    strand = 1,
    type = "CDS",
    stringsAsFactors = FALSE
  )
  entity <- new_genome_entity(
    features_df = feat_df
  )

  feats <- features(entity, format = "data.frame")
  expect_s3_class(feats, "data.frame")
  expect_equal(nrow(feats), 1)
  expect_equal(feats$type[1], "CDS")
})

test_that("metadata accessor works", {
  entity <- new_genome_entity()

  # Get all metadata (Clean Architecture: metadata is a data.frame)
  meta <- metadata(entity)
  expect_s3_class(meta, "data.frame")
})

test_that("seqnames accessor works", {
  entity <- new_genome_entity(
    indices = list(
      seqnames = c("chr1", "chr2"),
      seqname_to_record = c(chr1 = 1, chr2 = 2),
      locus_tag_index = integer(),
      gene_index = list()
    )
  )

  names <- seqnames(entity)
  expect_type(names, "character")
  expect_equal(names, c("chr1", "chr2"))
})

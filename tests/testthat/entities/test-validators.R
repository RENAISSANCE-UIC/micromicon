test_that("validate_genomic_coordinates accepts valid coordinates", {
  expect_true(validate_genomic_coordinates(1, 100, "chr1"))
  expect_true(validate_genomic_coordinates(50, 150, "chr1", max_length = 200))
  expect_true(validate_genomic_coordinates(50, 50, "chr1"))  # Single position
})

test_that("validate_genomic_coordinates rejects negative start", {
  expect_error(
    validate_genomic_coordinates(-1, 100, "chr1"),
    "positive"
  )

  expect_error(
    validate_genomic_coordinates(0, 100, "chr1"),
    "positive"
  )
})

test_that("validate_genomic_coordinates rejects start > end", {
  expect_error(
    validate_genomic_coordinates(100, 50, "chr1"),
    "Start coordinate must be"
  )
})

test_that("validate_genomic_coordinates rejects out-of-bounds coordinates", {
  expect_error(
    validate_genomic_coordinates(1, 200, "chr1", max_length = 100),
    "exceed"
  )
})

test_that("validate_strand accepts valid strands", {
  expect_true(validate_strand("+"))
  expect_true(validate_strand("-"))
  expect_true(validate_strand("."))
  expect_true(validate_strand("*"))
  expect_true(validate_strand(NA_character_))
})

test_that("validate_strand rejects invalid strands", {
  expect_error(validate_strand("F"), "Invalid strand")
  expect_error(validate_strand("forward"), "Invalid strand")
  expect_error(validate_strand("1"), "Invalid strand")
})

test_that("validate_dna_sequence accepts valid sequences", {
  expect_true(validate_dna_sequence("ATCG"))
  expect_true(validate_dna_sequence("ATCGNNNATCG"))  # With ambiguity codes
  expect_true(validate_dna_sequence("ATCG--ATCG"))   # With gaps
  expect_true(validate_dna_sequence("atcg"))  # Lowercase (gets uppercased)
  expect_true(validate_dna_sequence("RYSWKM"))  # Ambiguity codes
})

test_that("validate_dna_sequence rejects invalid characters", {
  expect_error(validate_dna_sequence("ATCGX"), "Invalid nucleotide")
  expect_error(validate_dna_sequence("ATCG123"), "Invalid nucleotide")
  expect_error(validate_dna_sequence("ATCG "), "Invalid nucleotide")  # Space
})

test_that("validate_dna_sequence respects allow_gaps parameter", {
  expect_true(validate_dna_sequence("ATCG--ATCG", allow_gaps = TRUE))
  expect_error(validate_dna_sequence("ATCG--ATCG", allow_gaps = FALSE), "Invalid nucleotide")
})

test_that("validate_feature_type accepts recognized types", {
  expect_true(validate_feature_type("gene"))
  expect_true(validate_feature_type("CDS"))
  expect_true(validate_feature_type("mRNA"))
  expect_true(validate_feature_type("tRNA"))
  expect_true(validate_feature_type("rRNA"))
})

test_that("validate_feature_type warns on unrecognized types when allow_custom=TRUE", {
  expect_warning(
    validate_feature_type("my_custom_feature", allow_custom = TRUE),
    "Unrecognized"
  )
})

test_that("validate_feature_type rejects unrecognized types when allow_custom=FALSE", {
  expect_error(
    validate_feature_type("my_custom_feature", allow_custom = FALSE),
    "Unrecognized"
  )
})

test_that("validate_locus_tag accepts valid locus tags", {
  expect_true(validate_locus_tag("GENE001"))
  expect_true(validate_locus_tag("ABC_12345"))
  expect_true(validate_locus_tag("locus_tag_v2"))
  expect_true(validate_locus_tag("A"))  # Single character
})

test_that("validate_locus_tag rejects invalid characters", {
  expect_error(validate_locus_tag("GENE 001"), "Invalid locus tag")  # Space
  expect_error(validate_locus_tag("GENE-001"), "Invalid locus tag")  # Hyphen
  expect_error(validate_locus_tag("GENE.001"), "Invalid locus tag")  # Dot
})

test_that("validate_locus_tag rejects too-long tags", {
  long_tag <- paste(rep("A", 100), collapse = "")
  expect_error(validate_locus_tag(long_tag, max_length = 50), "too long")
})

test_that("validate_sequence_topology accepts valid topologies", {
  expect_true(validate_sequence_topology("linear"))
  expect_true(validate_sequence_topology("circular"))
  expect_true(validate_sequence_topology(NA_character_))
})

test_that("validate_sequence_topology rejects invalid topologies", {
  expect_error(validate_sequence_topology("supercoiled"), "Invalid topology")
  expect_error(validate_sequence_topology("relaxed"), "Invalid topology")
})

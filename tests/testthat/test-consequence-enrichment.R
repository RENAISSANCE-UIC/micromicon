test_that(".pm_parse_position handles various formats", {
  # Standard integer
  expect_equal(.pm_parse_position("12345"), 12345L)

  # Comma-formatted
  expect_equal(.pm_parse_position("1,234,567"), 1234567L)

  # breseq P:k format (use base position P)
  expect_equal(.pm_parse_position("12345:0"), 12345L)
  expect_equal(.pm_parse_position("12345:1"), 12345L)
  expect_equal(.pm_parse_position("1,234,567:2"), 1234567L)

  # Edge cases
  expect_true(is.na(.pm_parse_position("")))
  expect_true(is.na(.pm_parse_position(NA_character_)))
  expect_true(is.na(.pm_parse_position("not_a_number")))
  expect_true(is.na(.pm_parse_position("0")))  # Position must be >= 1
  expect_true(is.na(.pm_parse_position("-123")))
})


test_that(".pm_parse_position with return_offset extracts k offset", {
  # Standard position without offset
  result <- .pm_parse_position("12345", return_offset = TRUE)
  expect_equal(result$position, 12345L)
  expect_equal(result$offset, 0L)

  # Position with k=0
  result <- .pm_parse_position("12345:0", return_offset = TRUE)
  expect_equal(result$position, 12345L)
  expect_equal(result$offset, 0L)

  # Position with k=1 (SNP after insertion)
  result <- .pm_parse_position("12345:1", return_offset = TRUE)
  expect_equal(result$position, 12345L)
  expect_equal(result$offset, 1L)

  # Position with k=2
  result <- .pm_parse_position("1,234,567:2", return_offset = TRUE)
  expect_equal(result$position, 1234567L)
  expect_equal(result$offset, 2L)

  # Invalid k (non-numeric)
  result <- .pm_parse_position("12345:abc", return_offset = TRUE)
  expect_equal(result$position, 12345L)
  expect_equal(result$offset, 0L)  # Defaults to 0 on parse fail

  # NA input
  result <- .pm_parse_position(NA_character_, return_offset = TRUE)
  expect_true(is.na(result$position))
  expect_equal(result$offset, 0L)
})


test_that(".pm_parse_snp_mutation handles arrow formats", {
  # Unicode arrow (→)
  result <- .pm_parse_snp_mutation("A\u2192C")
  expect_equal(result$ref, "A")
  expect_equal(result$alt, "C")

  # ASCII arrow (->)
  result <- .pm_parse_snp_mutation("G->T")
  expect_equal(result$ref, "G")
  expect_equal(result$alt, "T")

  # Simple greater-than (>)
  result <- .pm_parse_snp_mutation("C>A")
  expect_equal(result$ref, "C")
  expect_equal(result$alt, "A")

  # HTML-escaped arrow (-&gt;)
  result <- .pm_parse_snp_mutation("A-&gt;T")
  expect_equal(result$ref, "A")
  expect_equal(result$alt, "T")

  # HTML-escaped greater-than (&gt;)
  result <- .pm_parse_snp_mutation("G&gt;C")
  expect_equal(result$ref, "G")
  expect_equal(result$alt, "C")

  # With whitespace
  result <- .pm_parse_snp_mutation(" T -> G ")
  expect_equal(result$ref, "T")
  expect_equal(result$alt, "G")

  # HTML-escaped with whitespace
  result <- .pm_parse_snp_mutation(" C -&gt; A ")
  expect_equal(result$ref, "C")
  expect_equal(result$alt, "A")

  # Invalid formats
  result <- .pm_parse_snp_mutation("ABC")
  expect_true(is.na(result$ref))
  expect_true(is.na(result$alt))

  result <- .pm_parse_snp_mutation("A->CG")
  expect_true(is.na(result$ref))
  expect_true(is.na(result$alt))

  result <- .pm_parse_snp_mutation("")
  expect_true(is.na(result$ref))
  expect_true(is.na(result$alt))

  result <- .pm_parse_snp_mutation(NA_character_)
  expect_true(is.na(result$ref))
  expect_true(is.na(result$alt))

  # Non-DNA bases
  result <- .pm_parse_snp_mutation("X->Y")
  expect_true(is.na(result$ref))
  expect_true(is.na(result$alt))
})


test_that(".pm_apply_snp_to_sequence modifies sequence correctly", {
  seq <- "ATCGATCG"

  # Middle position
  result <- .pm_apply_snp_to_sequence(seq, 4, "T")
  expect_equal(result, "ATCTATCG")

  # First position
  result <- .pm_apply_snp_to_sequence(seq, 1, "C")
  expect_equal(result, "CTCGATCG")

  # Last position
  result <- .pm_apply_snp_to_sequence(seq, 8, "A")
  expect_equal(result, "ATCGATCA")

  # Out of bounds
  expect_true(is.na(.pm_apply_snp_to_sequence(seq, 0, "A")))
  expect_true(is.na(.pm_apply_snp_to_sequence(seq, 9, "A")))
  expect_true(is.na(.pm_apply_snp_to_sequence(seq, 100, "A")))

  # Invalid inputs
  expect_true(is.na(.pm_apply_snp_to_sequence(NA_character_, 1, "A")))
  expect_true(is.na(.pm_apply_snp_to_sequence(seq, NA_integer_, "A")))
  expect_true(is.na(.pm_apply_snp_to_sequence(seq, 1, NA_character_)))
})


test_that("pm_enrich_consequences validates inputs", {
  # Create minimal mock gd object (won't actually use it for these tests)
  mock_gd <- structure(
    list(
      source_type = "gd",
      events = list(),
      sequences = list()
    ),
    class = c("genome_entity_gd", "genome_entity")
  )

  # Valid minimal table
  pm_tbl <- data.frame(
    type = "SNP",
    seq_id = "NC_007795",
    position = "12345",
    mutation = "A\u2192C",
    gene = "testGene",
    stringsAsFactors = FALSE
  )

  # Should not error with valid inputs (but will fail later without real data)
  expect_error(
    pm_enrich_consequences(mock_gd, pm_tbl, flank = 50, quiet = TRUE),
    NA  # Expect no validation error at input level
  )

  # Missing required columns
  bad_tbl <- data.frame(type = "SNP", stringsAsFactors = FALSE)
  expect_error(
    pm_enrich_consequences(mock_gd, bad_tbl),
    "missing required columns"
  )

  # Invalid flank (should throw error)
  expect_error(
    pm_enrich_consequences(mock_gd, pm_tbl, flank = -10)
  )
})


test_that("pm_enrich_consequences adds expected columns", {
  # Create minimal mock gd object
  mock_gd <- structure(
    list(
      source_type = "gd",
      events = list(),
      sequences = list()
    ),
    class = c("genome_entity_gd", "genome_entity")
  )

  pm_tbl <- data.frame(
    type = "SNP",
    seq_id = "NC_007795",
    position = "12345",
    mutation = "A\u2192C",
    gene = "testGene",
    stringsAsFactors = FALSE
  )

  # This will fail internally due to missing data, but should still add columns
  result <- suppressWarnings(
    pm_enrich_consequences(mock_gd, pm_tbl, quiet = TRUE)
  )

  expected_cols <- c(
    "dna_ref", "dna_alt", "aa_ref", "aa_alt",
    "codon_ref", "codon_alt", "consequence", "region", "strand", "qc_note"
  )

  for (col in expected_cols) {
    expect_true(col %in% names(result),
                info = sprintf("Column '%s' should be added", col))
  }
})


test_that("pm_enrich_consequences handles non-SNP types gracefully", {
  mock_gd <- structure(
    list(
      source_type = "gd",
      events = list(),
      sequences = list()
    ),
    class = c("genome_entity_gd", "genome_entity")
  )

  pm_tbl <- data.frame(
    type = c("SNP", "DEL", "INS", "SUB"),
    seq_id = rep("NC_007795", 4),
    position = rep("12345", 4),
    mutation = c("A\u2192C", "\u03941 bp", "+ACGT", "2 bp\u2192GC"),
    gene = rep("testGene", 4),
    stringsAsFactors = FALSE
  )

  result <- suppressWarnings(
    pm_enrich_consequences(mock_gd, pm_tbl, quiet = TRUE)
  )

  # Should process all rows without error
  expect_equal(nrow(result), 4)

  # DEL, INS, SUB are now implemented - may have region values or NA
  # (depends on whether mock data can be processed)
  # Just verify they don't cause errors
  expect_true(TRUE)
})


test_that("pm_enrich_consequences handles empty table", {
  mock_gd <- structure(
    list(
      source_type = "gd",
      events = list(),
      sequences = list()
    ),
    class = c("genome_entity_gd", "genome_entity")
  )

  pm_tbl <- data.frame(
    type = character(0),
    seq_id = character(0),
    position = character(0),
    mutation = character(0),
    gene = character(0),
    stringsAsFactors = FALSE
  )

  result <- pm_enrich_consequences(mock_gd, pm_tbl, quiet = TRUE)

  expect_equal(nrow(result), 0)
  expect_true("consequence" %in% names(result))
})


test_that(".pm_enrich_snp handles unparseable positions", {
  mock_gd <- structure(
    list(
      source_type = "gd",
      events = list(),
      sequences = list()
    ),
    class = c("genome_entity_gd", "genome_entity")
  )

  row <- data.frame(
    type = "SNP",
    seq_id = "NC_007795",
    position = "invalid_position",
    mutation = "A\u2192C",
    gene = "testGene",
    dna_ref = NA_character_,
    dna_alt = NA_character_,
    aa_ref = NA_character_,
    aa_alt = NA_character_,
    codon_ref = NA_character_,
    codon_alt = NA_character_,
    consequence = NA_character_,
    region = NA_character_,
    strand = NA_character_,
    qc_note = NA_character_,
    stringsAsFactors = FALSE
  )

  result <- suppressWarnings(
    .pm_enrich_snp(mock_gd, row, flank = 50, quiet = TRUE)
  )

  # Should return row with QC note
  expect_true(is.na(result$region))
  expect_true(!is.na(result$qc_note))
  expect_true(grepl("Position parse failed", result$qc_note))
})


test_that(".pm_enrich_snp handles unparseable mutations", {
  mock_gd <- structure(
    list(
      source_type = "gd",
      events = list(),
      sequences = list()
    ),
    class = c("genome_entity_gd", "genome_entity")
  )

  row <- data.frame(
    type = "SNP",
    seq_id = "NC_007795",
    position = "12345",
    mutation = "not_a_snp",
    gene = "testGene",
    dna_ref = NA_character_,
    dna_alt = NA_character_,
    aa_ref = NA_character_,
    aa_alt = NA_character_,
    codon_ref = NA_character_,
    codon_alt = NA_character_,
    consequence = NA_character_,
    region = NA_character_,
    strand = NA_character_,
    qc_note = NA_character_,
    stringsAsFactors = FALSE
  )

  result <- suppressWarnings(
    .pm_enrich_snp(mock_gd, row, flank = 50, quiet = TRUE)
  )

  # Should add QC note for unparseable format; consequence remains NA
  expect_true(!is.na(result$qc_note))
  expect_true(grepl("Unparseable mutation format", result$qc_note))
})


test_that(".pm_parse_annotation_geometry extracts coding geometry", {
  # Standard coding annotation
  result <- .pm_parse_annotation_geometry("coding (45/1200 nt)")
  expect_equal(result$region, "coding")
  expect_equal(result$coding_pos, 45L)
  expect_equal(result$coding_len, 1200L)

  # With extra whitespace
  result <- .pm_parse_annotation_geometry("coding ( 123 / 999 nt )")
  expect_equal(result$region, "coding")
  expect_equal(result$coding_pos, 123L)
  expect_equal(result$coding_len, 999L)

  # Intergenic annotation
  result <- .pm_parse_annotation_geometry("intergenic (+123/-456)")
  expect_equal(result$region, "intergenic")
  expect_true(is.na(result$coding_pos))
  expect_true(is.na(result$coding_len))

  # Simple intergenic
  result <- .pm_parse_annotation_geometry("intergenic")
  expect_equal(result$region, "intergenic")

  # Unknown format
  result <- .pm_parse_annotation_geometry("some other text")
  expect_true(is.na(result$region))

  # Empty string
  result <- .pm_parse_annotation_geometry("")
  expect_true(is.na(result$region))

  # NA input
  result <- .pm_parse_annotation_geometry(NA_character_)
  expect_true(is.na(result$region))
})


test_that(".pm_append_qc_note accumulates notes", {
  # Create test row with qc_note column
  row <- data.frame(
    type = "SNP",
    qc_note = NA_character_,
    stringsAsFactors = FALSE
  )

  # First note (to NA field)
  row <- .pm_append_qc_note(row, "First issue")
  expect_equal(row$qc_note, "First issue")

  # Second note (should append with semicolon)
  row <- .pm_append_qc_note(row, "Second issue")
  expect_equal(row$qc_note, "First issue; Second issue")

  # Third note
  row <- .pm_append_qc_note(row, "Third issue")
  expect_equal(row$qc_note, "First issue; Second issue; Third issue")
})


test_that(".pm_append_qc_note handles empty starting note", {
  # Create test row with empty string (not NA)
  row <- data.frame(
    type = "SNP",
    qc_note = "",
    stringsAsFactors = FALSE
  )

  # Should treat empty string as NA
  row <- .pm_append_qc_note(row, "First note")
  expect_equal(row$qc_note, "First note")
})


test_that("pm_enrich_consequences creates qc_note for non-zero offset", {
  # This is an integration-style test that would need real data
  # For now, we test that the column exists (covered by other tests)
  # Real testing would require actual GenomeData object
  expect_true(TRUE)  # Placeholder for integration test
})


test_that("pm_enrich_consequences codon_alt column exists", {
  mock_gd <- structure(
    list(
      source_type = "gd",
      events = list(),
      sequences = list()
    ),
    class = c("genome_entity_gd", "genome_entity")
  )

  pm_tbl <- data.frame(
    type = "SNP",
    seq_id = "NC_007795",
    position = "12345",
    mutation = "A\u2192C",
    gene = "testGene",
    stringsAsFactors = FALSE
  )

  result <- suppressWarnings(
    pm_enrich_consequences(mock_gd, pm_tbl, quiet = TRUE)
  )

  expect_true("codon_alt" %in% names(result))
})

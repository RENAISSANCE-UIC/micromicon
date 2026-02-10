test_that("get_gene_dna extracts DNA sequence correctly", {
  skip_if(!file.exists("inst/extdata/test.gbk"), "Test file not found")

  entity <- read_genome("inst/extdata/test.gbk")

  # Get first CDS feature
  cds_features <- features(entity, type = "CDS")
  skip_if(nrow(cds_features) == 0, "No CDS features in test data")

  # Test by feature index
  dna1 <- get_gene_dna(entity, 1)
  expect_type(dna1, "character")
  expect_true(nchar(dna1) > 0)
  expect_true(grepl("^[ACGTN]+$", dna1, ignore.case = TRUE))

  # Test by locus_tag if available
  if ("locus_tag" %in% names(cds_features) && !is.na(cds_features$locus_tag[1])) {
    dna2 <- get_gene_dna(entity, cds_features$locus_tag[1])
    expect_type(dna2, "character")
    expect_true(nchar(dna2) > 0)
  }

  # Test by gene name if available
  if ("gene" %in% names(cds_features) && !is.na(cds_features$gene[1])) {
    dna3 <- get_gene_dna(entity, cds_features$gene[1])
    expect_type(dna3, "character")
    expect_true(nchar(dna3) > 0)
  }

  # Test that sequence length matches coordinates
  feat <- cds_features[1, ]
  expected_len <- feat$end - feat$start + 1
  expect_equal(nchar(dna1), expected_len)
})

test_that("get_gene_aa translates CDS to amino acids", {
  skip_if(!file.exists("inst/extdata/test.gbk"), "Test file not found")

  entity <- read_genome("inst/extdata/test.gbk")

  # Get first CDS feature
  cds_features <- features(entity, type = "CDS")
  skip_if(nrow(cds_features) == 0, "No CDS features in test data")

  # Get amino acid sequence
  aa <- get_gene_aa(entity, 1)
  expect_type(aa, "character")
  expect_true(nchar(aa) > 0)

  # Check that it's amino acids (standard one-letter codes)
  expect_true(grepl("^[ACDEFGHIKLMNPQRSTVWY*X]+$", aa))

  # Check length ratio (DNA/3 ≈ AA length)
  dna <- get_gene_dna(entity, 1)
  expect_equal(nchar(aa), floor(nchar(dna) / 3))

  # Test with custom genetic code
  aa_custom <- get_gene_aa(entity, 1, genetic_code = "1")
  expect_type(aa_custom, "character")
  expect_true(nchar(aa_custom) > 0)
})

test_that("get_roi_dna extracts arbitrary regions", {
  skip_if(!file.exists("inst/extdata/test.gbk"), "Test file not found")

  entity <- read_genome("inst/extdata/test.gbk")
  seqs <- sequences(entity)
  seqname <- names(seqs)[1]
  seq_len <- nchar(seqs[1])

  # Extract a region
  roi <- get_roi_dna(entity, chrom = seqname, start = 1, end = 100, strand = "+")
  expect_type(roi, "character")
  expect_equal(nchar(roi), 100)

  # Verify it matches the source sequence
  expected <- substr(seqs[1], 1, 100)
  expect_equal(roi, expected)

  # Test minus strand (should reverse complement)
  roi_minus <- get_roi_dna(entity, chrom = seqname, start = 1, end = 100, strand = "-")
  expect_type(roi_minus, "character")
  expect_equal(nchar(roi_minus), 100)
  expect_false(roi == roi_minus)  # Should be different
})

test_that("translate_dna translates correctly with different frames", {
  # Test sequence: ATG (M) + 6 codons
  dna <- "ATGAAATTTCCCGGGTTTAAATAG"

  # Frame 1
  aa1 <- translate_dna(dna, frame = 1, genetic_code = "11")
  expect_type(aa1, "character")
  expect_equal(nchar(aa1), 8)  # 24/3 = 8 codons

  # Frame 2 (skip first base)
  aa2 <- translate_dna(dna, frame = 2, genetic_code = "11")
  expect_type(aa2, "character")
  expect_equal(nchar(aa2), 7)  # 23/3 = 7 codons (1 base skipped)

  # Frame 3 (skip first two bases)
  aa3 <- translate_dna(dna, frame = 3, genetic_code = "11")
  expect_type(aa3, "character")
  expect_equal(nchar(aa3), 7)  # 22/3 = 7 codons (2 bases skipped)

  # Different frames should give different results
  expect_false(aa1 == aa2)
  expect_false(aa1 == aa3)

  # Test genetic code 1 vs 11 (should be same for most codons)
  aa_gc1 <- translate_dna(dna, frame = 1, genetic_code = "1")
  aa_gc11 <- translate_dna(dna, frame = 1, genetic_code = "11")
  expect_type(aa_gc1, "character")
  expect_type(aa_gc11, "character")
})

test_that("map_genomic_to_cds and map_cds_to_genomic work correctly", {
  skip_if(!file.exists("inst/extdata/test.gbk"), "Test file not found")

  entity <- read_genome("inst/extdata/test.gbk")

  # Get a CDS on plus strand if available
  cds_features <- features(entity, type = "CDS")
  skip_if(nrow(cds_features) == 0, "No CDS features in test data")

  plus_strand_features <- cds_features[cds_features$strand == "+", ]
  if (nrow(plus_strand_features) > 0) {
    feat <- plus_strand_features[1, ]

    # Test genomic to CDS mapping
    genomic_pos <- feat$start + 10
    cds_pos <- map_genomic_to_cds(entity, 1, genomic_pos)
    expect_type(cds_pos, "integer")
    expect_equal(cds_pos, 11L)  # 10 + 1 (1-based)

    # Test inverse mapping
    genomic_back <- map_cds_to_genomic(entity, 1, cds_pos)
    expect_equal(genomic_back, genomic_pos)

    # Test position outside gene
    outside_pos <- feat$end + 100
    cds_outside <- map_genomic_to_cds(entity, 1, outside_pos)
    expect_true(is.na(cds_outside))
  }

  # Test minus strand if available
  minus_strand_features <- cds_features[cds_features$strand == "-", ]
  if (nrow(minus_strand_features) > 0) {
    idx <- which(entity$features$strand == "-" & entity$features$type == "CDS")[1]
    feat <- entity$features[idx, ]

    # For minus strand, CDS position counts from end
    genomic_pos <- feat$end - 10
    cds_pos <- map_genomic_to_cds(entity, idx, genomic_pos)
    expect_type(cds_pos, "integer")
    expect_equal(cds_pos, 11L)  # 10 + 1

    # Test inverse
    genomic_back <- map_cds_to_genomic(entity, idx, cds_pos)
    expect_equal(genomic_back, genomic_pos)
  }
})

test_that("gene_info returns correct metadata", {
  skip_if(!file.exists("inst/extdata/test.gbk"), "Test file not found")

  entity <- read_genome("inst/extdata/test.gbk")

  # Get first CDS
  cds_features <- features(entity, type = "CDS")
  skip_if(nrow(cds_features) == 0, "No CDS features in test data")

  info <- gene_info(entity, 1)

  # Check required fields
  expect_type(info, "list")
  expect_true("chrom" %in% names(info))
  expect_true("start" %in% names(info))
  expect_true("end" %in% names(info))
  expect_true("strand" %in% names(info))
  expect_true("type" %in% names(info))
  expect_true("length" %in% names(info))
  expect_true("frame" %in% names(info))

  # Check values
  feat <- cds_features[1, ]
  expect_equal(info$chrom, feat$seqname)
  expect_equal(info$start, feat$start)
  expect_equal(info$end, feat$end)
  expect_equal(info$length, feat$end - feat$start + 1)

  # Check optional fields
  if (!is.na(feat$gene)) {
    expect_equal(info$gene, feat$gene)
  }
  if (!is.na(feat$locus_tag)) {
    expect_equal(info$locus_tag, feat$locus_tag)
  }
})

test_that("validate_variant_in_gene validates correctly", {
  skip_if(!file.exists("inst/extdata/test.gbk"), "Test file not found")

  entity <- read_genome("inst/extdata/test.gbk")

  # Get first CDS
  cds_features <- features(entity, type = "CDS")
  skip_if(nrow(cds_features) == 0, "No CDS features in test data")

  feat <- cds_features[1, ]

  # Get the actual reference base at gene start
  seqs <- sequences(entity)
  ref_base <- substr(seqs[[feat$seqname]], feat$start, feat$start)

  # Valid variant (correct position and ref base)
  result <- validate_variant_in_gene(entity, 1, feat$start, ref_base)
  expect_true(result)
  expect_true(grepl("Valid", attr(result, "message")))

  # Invalid: position outside gene
  result_outside <- validate_variant_in_gene(entity, 1, feat$end + 100, "A")
  expect_false(result_outside)
  expect_true(grepl("outside gene bounds", attr(result_outside, "message")))

  # Invalid: wrong reference base
  wrong_base <- if (ref_base == "A") "T" else "A"
  result_wrong_ref <- validate_variant_in_gene(entity, 1, feat$start, wrong_base)
  expect_false(result_wrong_ref)
  expect_true(grepl("mismatch", attr(result_wrong_ref, "message")))
})

test_that("gene resolution works with different identifiers", {
  skip_if(!file.exists("inst/extdata/test.gbk"), "Test file not found")

  entity <- read_genome("inst/extdata/test.gbk")

  cds_features <- features(entity, type = "CDS")
  skip_if(nrow(cds_features) == 0, "No CDS features in test data")

  # Test with integer index
  dna_idx <- get_gene_dna(entity, 1)
  expect_type(dna_idx, "character")

  # Test with data.frame row
  dna_df <- get_gene_dna(entity, cds_features[1, , drop = FALSE])
  expect_equal(dna_idx, dna_df)

  # Test error handling for invalid index
  expect_error(
    get_gene_dna(entity, 99999),
    "out of range"
  )

  # Test error handling for non-existent gene name
  expect_error(
    get_gene_dna(entity, "NONEXISTENT_GENE_12345"),
    "not found"
  )
})

test_that("reverse_complement works correctly", {
  # Test basic complement
  seq <- "ATCG"
  rc <- reverse_complement(seq)
  expect_equal(rc, "CGAT")

  # Test with longer sequence
  seq2 <- "ATCGATCG"
  rc2 <- reverse_complement(seq2)
  expect_equal(rc2, "CGATCGAT")

  # Test with ambiguous bases
  seq3 <- "ATCGN"
  rc3 <- reverse_complement(seq3)
  expect_equal(rc3, "NCGAT")

  # Test double reverse complement (should give original)
  rc_rc <- reverse_complement(rc)
  expect_equal(rc_rc, seq)
})

test_that("genetic code tables work correctly", {
  # Test standard genetic code
  dna <- "ATGTAA"  # ATG = M, TAA = * (stop)
  aa <- translate_dna(dna, genetic_code = "1")
  expect_equal(aa, "M*")

  # Test bacterial genetic code
  aa_bac <- translate_dna(dna, genetic_code = "11")
  expect_equal(aa_bac, "M*")

  # Test error for unsupported genetic code
  expect_error(
    translate_dna(dna, genetic_code = "99"),
    "not implemented"
  )
})

test_that("End-to-end GenBank workflow works", {
  skip_if(!file.exists("inst/extdata/test.gbk"), "Test file not found")

  # Read GenBank file
  entity <- read_genome("inst/extdata/test.gbk")

  # Validate structure
  expect_s3_class(entity, "genome_entity")
  expect_true(length(entity$sequences$dna_raw) > 0)

  # Access data
  seqs <- sequences(entity)
  expect_type(seqs, "character")

  feats <- features(entity)
  expect_s3_class(feats, "data.frame")

  meta <- metadata(entity)
  expect_s3_class(meta, "data.frame")

  names <- seqnames(entity)
  expect_type(names, "character")

  # Print/summary should work
  expect_output(print(entity), "genome_entity")
  expect_output(summary(entity), "genome_entity Summary")

  # Export to FASTA
  temp_fasta <- tempfile(fileext = ".fasta")
  on.exit(unlink(temp_fasta), add = TRUE)

  result <- write_fasta(entity, temp_fasta)
  expect_true(file.exists(temp_fasta))

  # Read exported file
  lines <- readLines(temp_fasta)
  expect_true(any(grepl("^>", lines)))
  expect_true(any(grepl("^[ACGTN]+$", lines)))
})

test_that("End-to-end GFF3+FASTA workflow works", {
  skip_if(!file.exists("inst/extdata/test.gff3"), "Test GFF3 file not found")
  skip_if(!file.exists("inst/extdata/test.fasta"), "Test FASTA file not found")

  # Read GFF3 + FASTA
  entity <- read_genome(
    gff = "inst/extdata/test.gff3",
    fasta = "inst/extdata/test.fasta"
  )

  # Validate structure
  expect_s3_class(entity, "genome_entity")
  expect_true(length(entity$sequences$dna_raw) > 0)
  expect_true(nrow(entity$features) > 0)

  # Query features
  if (nrow(entity$features) > 0 && "type" %in% names(entity$features)) {
    genes <- search_features(entity, type = "gene")
    expect_s3_class(genes, "data.frame")
  }
})

test_that("Clean Architecture layers work independently", {
  # Test Entity layer
  entity <- new_genome_entity()
  expect_s3_class(entity, "genome_entity")
  expect_silent(validate_genome_entity(entity))

  # Test Gateway layer
  gbk_gateway <- create_genbank_gateway()
  expect_type(gbk_gateway, "list")
  expect_true("read" %in% names(gbk_gateway))

  fasta_gateway <- create_fasta_gateway(use_bioconductor = FALSE)
  expect_type(fasta_gateway, "list")
  expect_true("read" %in% names(fasta_gateway))
  expect_true("write" %in% names(fasta_gateway))

  # Test Use Case layer (with mock gateway)
  mock_gateway <- list(
    read = function(path) {
      list(
        list(
          metadata = list(
            locus = "TEST001",
            accession = "TEST001",
            length_bp = 16,
            mol_type = "DNA",
            topology = "linear"
          ),
          features = data.frame(
            type = "gene",
            start = 1,
            end = 10,
            strand = "+",
            gene = "testGene",
            locus_tag = "TEST_001",
            stringsAsFactors = FALSE
          ),
          sequence = "ATCGATCGATCGATCG"
        )
      )
    }
  )

  result <- execute_import_genbank(mock_gateway, "dummy.gbk")
  expect_s3_class(result, "genome_entity")
})

test_that("Backward compatibility is maintained", {
  skip_if(!file.exists("inst/extdata/test.gbk"), "Test file not found")

  # Old API: read_gbk with return_entity = FALSE
  gbk_list <- read_gbk("inst/extdata/test.gbk", return_entity = FALSE)
  expect_type(gbk_list, "list")
  expect_true(length(gbk_list) > 0)

  # Should have old structure
  expect_true("metadata" %in% names(gbk_list[[1]]))
  expect_true("features" %in% names(gbk_list[[1]]))
  expect_true("sequence" %in% names(gbk_list[[1]]))

  # Convert to entity
  entity <- gbk_to_entity(gbk_list)
  expect_s3_class(entity, "genome_entity")

  # Convert back to legacy format
  legacy <- entity_to_gbk_list(entity)
  expect_type(legacy, "list")
})

test_that("Error handling works correctly", {
  # Invalid file path
  expect_error(
    read_genome("nonexistent.gbk"),
    "not found"
  )

  # Invalid entity
  bad_entity <- list(foo = "bar")
  class(bad_entity) <- "genome_entity"

  expect_error(
    validate_genome_entity(bad_entity),
    "Missing required components"
  )

  # Invalid input to accessors
  expect_error(
    sequences(list(not_an_entity = TRUE)),
    "not implemented for class"
  )
})

test_that("Dependency injection pattern works", {
  # Create a custom gateway
  custom_gateway <- list(
    read = function(path) {
      list(
        list(
          metadata = list(locus = "CUSTOM"),
          features = data.frame(),
          sequence = "ATCG"
        )
      )
    }
  )

  # Use case should work with custom gateway
  result <- execute_import_genbank(custom_gateway, "test.gbk")
  expect_s3_class(result, "genome_entity")
})

test_that("Query functions work with use cases", {
  # Create a test entity
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
      gene = "testGene",
      locus_tag = "TEST_001",
      stringsAsFactors = FALSE
    ),
    metadata_df = data.frame(
      seqname = "chr1",
      length = 16,
      stringsAsFactors = FALSE
    ),
    indices_list = list(
      seqnames = "chr1",
      locus_tag_index = integer(),
      gene_index = list()
    )
  )

  # Test use case layer directly
  result <- execute_extract_sequences_by_coords(
    entity,
    seqname = "chr1",
    start = 1,
    end = 4,
    list()
  )
  expect_type(result, "character")
  expect_equal(as.character(result), "ATCG")

  # Test use case for name extraction
  result2 <- execute_extract_sequences_by_name(
    entity,
    pattern = "testGene",
    fields = c("gene", "locus_tag", "product"),
    list()
  )
  expect_type(result2, "character")
  expect_equal(length(result2), 1)

  # Search features works with data.frame
  genes <- search_features(entity, type = "gene")
  expect_s3_class(genes, "data.frame")
  expect_equal(nrow(genes), 1)
})

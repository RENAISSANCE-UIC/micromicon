test_that("write_gff3 warns for GenBank-sourced data", {
  # Mock genome_entity with GenBank source
  genome <- list(
    sequences = list(dna_raw = c(seq1 = "ATCG")),
    features = data.frame(
      seqname = "seq1",
      type = "CDS",
      start = 1,
      end = 4,
      strand = "+",
      stringsAsFactors = FALSE
    ),
    metadata = data.frame(
      seqname = "seq1",
      length_bp = 4,
      stringsAsFactors = FALSE
    ),
    indices = list(seqnames = "seq1")
  )
  class(genome) <- "genome_entity"
  attr(genome, "import_source") <- "genbank"

  # Expect warning when exporting
  expect_message(
    write_gff3(genome, tempfile()),
    "GenBank-sourced data to GFF3 LOSES metadata"
  )
})

test_that("write_fasta warns for GenBank-sourced data", {
  # Mock genome_entity with GenBank source
  genome <- list(
    sequences = list(dna_raw = c(seq1 = "ATCG")),
    features = data.frame(
      seqname = "seq1",
      type = "CDS",
      start = 1,
      end = 4,
      strand = "+",
      stringsAsFactors = FALSE
    ),
    metadata = data.frame(
      seqname = "seq1",
      length_bp = 4,
      stringsAsFactors = FALSE
    ),
    indices = list(seqnames = "seq1")
  )
  class(genome) <- "genome_entity"
  attr(genome, "import_source") <- "genbank"

  # Expect warning when exporting
  expect_message(
    write_fasta(genome, tempfile()),
    "GenBank-sourced data to FASTA LOSES metadata"
  )
})

test_that("write_gff3 does not warn for GFF3-sourced data", {
  # Mock genome_entity with GFF3 source
  genome <- list(
    sequences = list(dna_raw = c(seq1 = "ATCG")),
    features = data.frame(
      seqname = "seq1",
      type = "CDS",
      start = 1,
      end = 4,
      strand = "+",
      stringsAsFactors = FALSE
    ),
    metadata = data.frame(
      seqname = "seq1",
      length_bp = 4,
      stringsAsFactors = FALSE
    ),
    indices = list(seqnames = "seq1")
  )
  class(genome) <- "genome_entity"
  attr(genome, "import_source") <- "gff3_fasta"

  # Capture all output
  output <- capture.output(
    result <- write_gff3(genome, tempfile()),
    type = "message"
  )

  # Should not contain warning about metadata loss
  expect_false(any(grepl("GenBank-sourced data to GFF3 LOSES metadata", output)))
})

test_that("write_fasta does not warn for GFF3-sourced data", {
  # Mock genome_entity with GFF3 source
  genome <- list(
    sequences = list(dna_raw = c(seq1 = "ATCG")),
    features = data.frame(
      seqname = "seq1",
      type = "CDS",
      start = 1,
      end = 4,
      strand = "+",
      stringsAsFactors = FALSE
    ),
    metadata = data.frame(
      seqname = "seq1",
      length_bp = 4,
      stringsAsFactors = FALSE
    ),
    indices = list(seqnames = "seq1")
  )
  class(genome) <- "genome_entity"
  attr(genome, "import_source") <- "gff3_fasta"

  # Capture all output
  output <- capture.output(
    result <- write_fasta(genome, tempfile()),
    type = "message"
  )

  # Should not contain warning about metadata loss
  expect_false(any(grepl("GenBank-sourced data to FASTA LOSES metadata", output)))
})

test_that("warnings can be suppressed with option", {
  # Mock genome_entity with GenBank source
  genome <- list(
    sequences = list(dna_raw = c(seq1 = "ATCG")),
    features = data.frame(
      seqname = "seq1",
      type = "CDS",
      start = 1,
      end = 4,
      strand = "+",
      stringsAsFactors = FALSE
    ),
    metadata = data.frame(
      seqname = "seq1",
      length_bp = 4,
      stringsAsFactors = FALSE
    ),
    indices = list(seqnames = "seq1")
  )
  class(genome) <- "genome_entity"
  attr(genome, "import_source") <- "genbank"

  # Suppress warnings
  withr::with_options(
    list(micromicon.warn_export = FALSE),
    {
      # Capture all output
      output <- capture.output(
        result <- write_gff3(genome, tempfile()),
        type = "message"
      )

      # Should not contain warning
      expect_false(any(grepl("GenBank-sourced data to GFF3 LOSES metadata", output)))
    }
  )
})

test_that("write_gff3 handles missing source attribute gracefully", {
  # Mock genome_entity without source attribute
  genome <- list(
    sequences = list(dna_raw = c(seq1 = "ATCG")),
    features = data.frame(
      seqname = "seq1",
      type = "CDS",
      start = 1,
      end = 4,
      strand = "+",
      stringsAsFactors = FALSE
    ),
    metadata = data.frame(
      seqname = "seq1",
      length_bp = 4,
      stringsAsFactors = FALSE
    ),
    indices = list(seqnames = "seq1")
  )
  class(genome) <- "genome_entity"
  # No import_source attribute set

  # Should not error
  expect_no_error(write_gff3(genome, tempfile()))
})

test_that("write_fasta handles missing source attribute gracefully", {
  # Mock genome_entity without source attribute
  genome <- list(
    sequences = list(dna_raw = c(seq1 = "ATCG")),
    features = data.frame(
      seqname = "seq1",
      type = "CDS",
      start = 1,
      end = 4,
      strand = "+",
      stringsAsFactors = FALSE
    ),
    metadata = data.frame(
      seqname = "seq1",
      length_bp = 4,
      stringsAsFactors = FALSE
    ),
    indices = list(seqnames = "seq1")
  )
  class(genome) <- "genome_entity"
  # No import_source attribute set

  # Should not error
  expect_no_error(write_fasta(genome, tempfile()))
})

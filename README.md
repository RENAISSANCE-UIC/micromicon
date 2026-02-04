# micRomicon

**A format-agnostic microbial genomics toolkit for R**

`micromicon` is a clean-architecture toolkit for reading, representing, and examining microbial genomes in R. Whether a genome arrives through GenBank, GFF3+FASTA, or another representation—breseq outputs, variant calls, curated mutation tables, or other analyses that encode sequence change—`micromicon` resolves each into a unified `genome_entity` for downstream interrogation, including future support for assessing the functional consequences of observed mutations.

## Why micRomicon?

We wanted a free, open-source toolkit that worked naturally for R users and lowered the friction of moving among file formats commonly used in microbial genomics. GenBank, GFF3+FASTA, and mutation-oriented outputs each bring their own structural hurdles, and the parsing logic for these is often scattered across different packages and domains. After decades of doing this the old way, we wanted a dedicated system to do the parsing and formatting for us, so we could reroute cognitive bandwidth to doing the actual science. 

`micromicon` consolidates these inputs into a single, stable representation (`genome_entity`) so that import, storage, query, and export operations follow the same patterns regardless of where the data originated. What began as a collection of convenience wrappers has grown into a format-agnostic foundation for routine bacterial genome analysis, with space reserved for future tooling, including variant-aware workflows and functional consequence interpretation.

## Installation

```r
# Install from GitHub
devtools::install_github("RENAISSANCE-UIC/micromicon")

# Or from local source
devtools::install_local("/path/to/micromicon")
```

### Optional (But Recommended): Bioconductor Integration

For advanced features (GRanges, DNAStringSet):

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GenomicRanges", "Biostrings", "rtracklayer"))
```

## Quick Start

```r
library(micromicon)

# Read any format
genome <- read_genome("ecoli.gbk")  # or GFF3+FASTA

# Explore
features(genome, type = "CDS")
sequences(genome)
genome_metadata(genome)

# Query
ampC <- search_features(genome, pattern = "ampC", type = "CDS")

# Extract
protein <- ampC$translation
dna <- extract_by_coords(genome, ampC$seqname, ampC$start, ampC$end)

# Export
write_gff3(genome, "output.gff3")
write_fasta(genome, "output.fasta")
```

## Core Features

### Format-Agnostic Input

Read genomes from multiple formats into a unified representation:

```r
# GenBank
genome1 <- read_genome("ecoli.gbk")

# GFF3 + FASTA
genome2 <- read_genome(gff = "ecoli.gff3", fasta = "ecoli.fna")

# All become genome_entity
class(genome1)  # "genome_entity"
class(genome2)  # "genome_entity"
```

### Search and Filter

Find features by type, name, or location:

```r
# All CDS features
cds <- features(genome, type = "CDS")

# Search by gene name
dna_genes <- search_features(genome, pattern = "dna")

# Features in a region
region <- search_features(genome,
  seqname = "chr1",
  start = 1000,
  end = 5000
)

# Helper for complex filters
filtered <- feat_filter(features(genome),
  type = "CDS",
  name = "lac"
)
```

### Extract Sequences

Pull sequences by coordinates or feature names:

```r
# By coordinates
seq <- extract_by_coords(genome, "chr1", 1000, 2000)

# By gene name
gene_seq <- extract_by_name(genome, "ampC")

# With translation
protein <- extract_by_name(genome, "ampC",
  type = "CDS",
  translate = TRUE
)

# Genomic regions with context
roi <- roi_coords("chr1", 1000, 2000, "+")
```

### Export to Standard Formats

Round-trip between formats:

```r
# Read GenBank
genome <- read_genome("input.gbk")

# Export to GFF3 + FASTA
write_gff3(genome, "output.gff3")
write_fasta(genome, "output.fasta")

# Re-import (lossless where specs allow)
genome2 <- read_genome(gff = "output.gff3", fasta = "output.fasta")
```

### Local BLASTP Integration (PROVISIONAL)

Compare proteins against local databases:

```r
# Extract protein sequence
ampC <- search_features(genome, pattern = "ampC", type = "CDS")
protein <- ampC$translation

# BLASTP against local SwissProt
hits <- blast_protein(protein, database = "swissprot")

# Filter high-quality matches
top_hits <- reduce_hits(hits,
  min_qcov = 80,
  min_pident = 50,
  besthit = FALSE,
  max_per_query = 5
)
```

**Requirements**: Local BLAST+ and database. See [BLAST Setup Guide](BLAST_SETUP.md).

## Design Philosophy

### Clean Architecture

`micromicon` follows Clean Architecture principles:

- **Entities**: `genome_entity` as the core domain object
- **Use Cases**: Pure functions for genome operations
- **Controllers**: S3 generics for extensibility
- **Gateways**: Format-specific parsing (GenBank, GFF3, etc.)

This separation ensures:
- Format-agnostic operations
- Easy testing and validation
- Extensibility for future formats
- No framework lock-in

### Functional Style

All operations return new data—never mutate inputs:

```r
# Functional pipeline
genome |>
  search_features(pattern = "lac", type = "CDS") |>
  extract_by_name(genome, ., translate = TRUE)
```

### S3 Generics for Extensibility

Core accessors use S3 dispatch, enabling custom implementations:

```r
# Built-in: in-memory genome
genome <- read_genome("ecoli.gbk")
features(genome)  # Returns data.frame

# Future: remote database genome
genome_remote <- connect_to_db(...)
features(genome_remote)  # Queries database, same interface
```

## Documentation

- **Getting Started**: See the [vignette](vignettes/micromicon-intro.Rmd)
- **BLAST Setup**: [BLAST_SETUP.md](BLAST_SETUP.md)
- **Function Reference**: `?read_genome`, `?search_features`, etc.

## Extensibility

`micromicon` is designed to be extended:

### Custom Genome Sources

Implement methods for new `genome_entity` subclasses:

```r
# Define a remote genome class
genome_entity_remote <- structure(
  list(connection = db_conn, cache = list()),
  class = c("genome_entity_remote", "genome_entity")
)

# Implement methods
features.genome_entity_remote <- function(x, ...) {
  query_database(x$connection, "SELECT * FROM features")
}

# Works seamlessly
features(remote_genome)  # Dispatches to your method
```

### Future Directions

`micromicon` is built to enable support:
- **Variant-aware workflows**: Track mutations across genomes
- **Functional consequence prediction**: Assess mutation impacts
- **breseq integration**: Import mutation calls
- **Comparative genomics**: Multi-genome operations

## Dependencies

**Required**:
- R ≥ 4.0
- cli, dplyr, readr, tibble

**Optional** (but recommended):
- Bioconductor packages (GenomicRanges, Biostrings) for advanced features
- BLAST+ for local BLASTP

## Contributing

We welcome contributions! Areas of interest:
- Additional format parsers (EMBL, PTT, etc.)
- Variant call integration (VCF, breseq)
- Functional annotation tools
- Performance optimizations

## Citation

If you use `micromicon` in your research, please cite:

```
[Citation information - to be added]
```

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Acknowledgments

`micromicon` builds on the shoulders of:
- Bioconductor ecosystem (GenomicRanges, Biostrings, rtracklayer)
- NCBI BLAST+ toolkit
- Clean Architecture principles by Robert C. Martin

---

**Questions?** Open an [issue](https://github.com/your-username/micromicon/issues) or check the [vignette](vignettes/micromicon-intro.Rmd).

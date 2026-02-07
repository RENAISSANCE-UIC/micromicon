# micRomicon: An Ostensibly Format-Agnostic Microbial Genomics Toolkit for R

This is the repo for `micromicon`, a clean-architecture toolkit for reading, representing, and examining microbial genomes in R. Whether the genome sequence and annotation arrives through GenBank, GFF3+FASTA, or another representation (such as annotated breseq genome difference files, variant calls, curated mutation tables, or other analyses that encode sequence change), `micromicon` will resolve each into a unified `genome_entity` object for downstream interrogation, including future support for transcriptomics and assessing the functional consequences of observed mutations.

## Why micRomicon?

We wanted a free, open-source toolkit that worked naturally for R users and lowered the barrier of moving among file formats commonly used in microbial genomics. GenBank, GFF3+FASTA, and mutation-oriented outputs each bring their own structural hurdles, and the parsing logic for these is often scattered across different packages and domains. After decades of doing this the old way, we wanted a dedicated system to do the parsing and formatting for us, so that we could reroute cognitive bandwidth toward doing the actual science. 

`micromicon` ingests and consolidates common genomics file formats into a single, stable representation (`genome_entity`) so that import, storage, query, and export operations follow the same patterns regardless of where the data originated. What began as a collection of convenience wrappers has grown into a format-agnostic foundation for routine bacterial genomics analysis, with space reserved for future tooling, including variant-aware workflows and functional consequence interpretation.

## Installation

```r
# Install from GitHub
devtools::install_github("RENAISSANCE-UIC/micromicon")

# Or from local source
devtools::install_local("/path/to/micromicon")
```

### Optional (But Recommended): Bioconductor Integration

For advanced features:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GenomicRanges", "Biostrings", "rtracklayer"))
```

## Quick Start

```r
library(micromicon)

# Read any GenBank or GFF3+FASTA
path <- system.file("extdata", "test_ampC.gbk", package = "micromicon")
genome <- read_genome(path) 

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

### Format Conversion Rules

**IMPORTANT**: Format conversion in micromicon is inentionally one-way.

#### ALLOWED: GenBank → GFF3+FASTA
```r
genome <- read_genome("reference.gbk")     # Read GenBank
write_gff3(genome, "output.gff3")          # Export GFF3 ✓
write_fasta(genome, "output.fasta")        # Export FASTA ✓
```

**WARNING**: This conversion LOSES metadata (organism, taxonomy, references, comments).

#### FORBIDDEN: GFF3+FASTA → GenBank 
```r
genome <- read_genome(gff = "anno.gff3", fasta = "seq.fasta")
# write_genbank(genome, "output.gbk")  # DOES NOT EXIST - FORBIDDEN
```

**Why forbidden?**
- GenBank provides rich metadata (organism, taxonomy, references) that GFF3+FASTA lacks
- Converting GFF3+FASTA to GenBank would produce incomplete/invalid files
- No `write_genbank()` function exists (this is intentional, not a missing feature)

**GenBank metadata that GFF3+FASTA cannot represent**:
- Organism name and taxonomic lineage
- Publication references and citations
- Curator comments and assembly information
- Accession numbers and version history
- Database cross-references (taxonomy IDs, etc.)
- Sequence topology (circular vs. linear)

**Best practice**: Always keep original GenBank files. Export to GFF3+FASTA only when required by downstream tools (genome browsers, annotation pipelines).

## Core Features

### (Ostensibly) Format-Agnostic Input

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

Export genome data for downstream tools:

```r
# Read GenBank (preserves all metadata)
genome <- read_genome("reference.gbk")

# Export to GFF3+FASTA (for genome browsers, pipelines)
write_gff3(genome, "output.gff3")
write_fasta(genome, "output.fasta")
```

**METADATA LOSS**: Exporting from GenBank to GFF3+FASTA loses organism, taxonomy, references, and comments. Keep your original GenBank file.

**NO REVERSE CONVERSION**: There is no `write_genbank()` function. You cannot convert GFF3+FASTA back to GenBank format. This is unimplemented because GFF3+FASTA lacks required metadata.

```r
# FORBIDDEN - Will never be implemented:
# write_genbank(genome, "output.gbk")  # Does not exist
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

`micromicon` follows some Clean Architecture principles:

- **Entities**: `genome_entity` as the core domain object
- **Use Cases**: Pure functions for genome operations
- **Controllers**: S3 generics for extensibility
- **Gateways**: Format-specific parsing (GenBank, GFF3, etc.)

This separation ensures:
- Format-agnostic operations
- Easy testing and validation
- Extensibility for future formats
- No framework lock-in


### The Object System Roadmap (S3 → S4/S7)

In this initial release, we employ S3 objects as the primary interface for `genome_entity` and the related operations. This was not (soley) from ideological allegiance, but from pragmatic concinnity: familiarity, low-friction extensibility, and excellent discoverability for most R users. At this stage of development, S3 provides a stalwart, easily inspectable substrate on which to stabilize the core idioms of `micromicon`.

As the codebase matures, we plan a gradatim migration toward more constrained and less mutable object systems, likely S4 or S7. Both provide stricter contracts, clearer invariants, and richer introspection, which will become increasingly important as the toolkit grows to support variant-aware workflows, functional consequence inference, and multi-genome comparative operations.
The precise destination remains open. S7, in particular, offers an appealing blend of rigor and simplicity (orthogonal to S4’s sometimes baroque formalism) while preserving the kind of explicitness that helps prevent accumulation of structural drift. Whatever the final form, the public interface will retain its present ethos: clean, predictable generics and minimal cognitive overhead for downstream analysis.

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

## Pipes, Routing, and Early Thinness

The public API of `micromicon` is intentionally thin. Most user-facing functions are S3 generics that act as the routing stubs that inspect the object’s class and dispatch to a backend implementation. We kept these wrappers deliberately minimal for now to preserve architectural flexibility as the backends proliferate (in-memory objects today; remote stores or on-disk indices later). Think of it as preemptive debt-control. `micromicon` emerged from an earlier, more noodle-shaped prototype. Rather than retrofit dispatch and purity later, we wanted to scaffold in routing early, ensuring that method boundaries, side-effect hygiene, and extension points would be made explicit from the start.

We have borrowed the excellent design idioms familiar in other ecosystems (clean architecture, dependency inversion, strict boundaries between entities and gateways). While R itself does not enforce these constraints, we hope that the idioms will make the public interface predictable, the backends swappable, and the migration to stricter object systems (S4 or S7) straightforward as our invariants become more fully specified. 


## Documentation

- **Getting Started**: See the [vignette](vignettes/micromicon-intro.Rmd)
- **BLAST Setup**: [BLAST_SETUP.md](BLAST_SETUP.md)
- **Function Reference**: `?read_genome`, `?search_features`, etc.

## Extensibility

`micromicon` is designed to be extended. 

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
Ackerman, W. (2026). micRomicon: An ostensibly format-agnostic microbial genomics toolkit for R. Zenodo. https://doi.org/10.5281/zenodo.18500345
```

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Acknowledgments

`micromicon` builds on the shoulders of:
- Bioconductor ecosystem (GenomicRanges, Biostrings, rtracklayer)
- NCBI BLAST+ toolkit
- Clean Architecture principles by Robert C. Martin ("Uncle Bob")

---

**Questions?** Open an [issue](https://github.com/your-username/micromicon/issues) or check the [vignette](vignettes/micromicon-intro.Rmd).

# Gene Utility Functions - Implementation Summary

## Overview

Successfully implemented a comprehensive set of gene utility functions for working with `genome_entity` objects in the micromicon package. All functions follow the Clean Architecture pattern with S3 generic methods.

## New Functions Implemented

### 1. `get_gene_dna(entity, gene, format = "character")`
Extract DNA sequence for a gene in 5'-to-3' canonical orientation.
- Automatically handles reverse complement for minus strand genes
- Supports multiple gene identifier types (index, locus_tag, gene name, ID, feature row)
- Optional Biostrings DNAString format output

### 2. `get_gene_aa(entity, gene, genetic_code = NULL)`
Get translated amino acid sequence for a CDS feature.
- Uses genetic code from feature's `transl_table` attribute or defaults to code 11 (bacterial)
- Custom genetic code can be specified
- Warns if feature is not a CDS

### 3. `get_roi_dna(entity, chrom, start, end, strand = "+")`
Extract DNA for an arbitrary genomic region.
- Convenience wrapper around `extract_by_coords()`
- Handles strand orientation automatically

### 4. `translate_dna(dna, frame = 1, genetic_code = "11")`
Enhanced translation function with frame and genetic code support.
- Supports reading frames 1, 2, and 3
- Supports NCBI genetic codes (currently codes 1 and 11)
- Extensible for additional genetic codes
- Handles unknown codons gracefully (translates to "X")

### 5. `map_genomic_to_cds(entity, gene, genomic_pos)`
Convert genomic position to CDS-relative position (1-based).
- Accounts for strand orientation
- Returns NA if position is outside gene bounds
- Works with both plus and minus strand genes

### 6. `map_cds_to_genomic(entity, gene, cds_pos)`
Convert CDS position to genomic coordinate (inverse of above).
- Accounts for strand orientation
- Returns NA if CDS position is out of range
- Validates input positions

### 7. `gene_info(entity, gene)`
Get comprehensive gene metadata.
- Returns list with: chrom, start, end, strand, type, length, cds_ranges, frame, gene, locus_tag, product
- Extracts reading frame from GFF3 phase field if available
- Useful for variant analysis and gene annotation

### 8. `validate_variant_in_gene(entity, gene, genomic_pos, ref_base)`
Validate genomic variant position and reference base.
- Checks if position is within gene boundaries
- Verifies reference base matches genome sequence
- Accounts for strand orientation
- Returns logical with descriptive message attribute

## Helper Functions

### Internal Helpers (not exported)
- `.resolve_gene()` - Resolves various gene identifier types to feature row
- `.get_genetic_code_table()` - Returns NCBI genetic code codon tables

### Exported Helpers (also used internally)
- `reverse_complement()` - Reverse complement DNA sequence (with IUPAC ambiguity codes)

## Gene Identifier Support

All functions support multiple ways to identify genes:
1. **Integer index**: `get_gene_dna(entity, 1)` - First feature
2. **Locus tag**: `get_gene_dna(entity, "b0001")` - By locus_tag field
3. **Gene name**: `get_gene_dna(entity, "dnaA")` - By gene field
4. **ID**: `get_gene_dna(entity, "gene-001")` - By ID field (GFF3)
5. **Feature row**: `get_gene_dna(entity, features(entity)[1,])` - Data.frame row

## Genetic Code Support

Currently implemented:
- **Code 1**: Standard (universal) genetic code
- **Code 11**: Bacterial, Archaeal, and Plant Plastid Code (default)

Easily extensible to add more codes (e.g., mitochondrial codes).

## Files Created/Modified

### New Files
- `R/generics_gene_utilities.R` - Main implementation (632 lines)
- `tests/testthat/test-gene-utilities.R` - Comprehensive test suite (313 lines)
- `examples/gene_utilities_demo.R` - Usage demonstration
- `vignettes/gene-utilities.Rmd` - Full documentation vignette

### Modified Files
- `R/use_cases_query_extract_sequences_by_name.R` - Updated to use new `translate_dna()`
- `R/use_cases_query_extract_sequences_by_coords.R` - Updated to use new `translate_dna()`
- `R/gd_parser.R` - Fixed duplicate roxygen documentation
- `NAMESPACE` - Auto-generated exports
- `man/*.Rd` - Auto-generated documentation files

## Testing

All tests passing (84 passes, 0 failures):
- Unit tests for all 8 main functions
- Edge case testing (invalid positions, wrong reference bases, etc.)
- Frame translation tests
- Genetic code tests
- Reverse complement tests
- Gene identifier resolution tests

## Documentation

Complete documentation provided:
- Roxygen2 documentation for all functions
- Full vignette with examples and workflows
- Demo script showing all features
- Integration with Mode A (reference genome) workflow

## Performance Considerations

- DNA sequences stored as character vectors for efficiency
- Position mappings use optimized 1-based coordinate system
- Functions validate inputs early to fail fast
- Automatic strand handling (no manual reverse complement needed)
- Genetic code translation optimized for bacterial genomes

## Integration with Existing Code

- Follows Clean Architecture pattern (Generics layer)
- Compatible with existing `genome_entity` structure
- Works seamlessly with `search_features()`, `extract_by_coords()`, etc.
- Fully compatible with Mode A (reference genome) workflow for breseq genome diffs
- No breaking changes to existing code

## Usage Example

```r
library(micromicon)

# Load genome
entity <- read_genome("genome.gbk")

# Get gene sequence
dna <- get_gene_dna(entity, "dnaA")
aa <- get_gene_aa(entity, "dnaA")

# Analyze variant
info <- gene_info(entity, "dnaA")
cds_pos <- map_genomic_to_cds(entity, "dnaA", 12345)
is_valid <- validate_variant_in_gene(entity, "dnaA", 12345, "A")

# Extract arbitrary region
roi <- get_roi_dna(entity, chrom = "chromosome", start = 1000, end = 2000)

# Translate with custom frame
aa <- translate_dna(dna, frame = 2, genetic_code = "1")
```

## Future Enhancements

Potential additions for future versions:
1. Support for genes with introns (multi-exon CDS)
2. Additional NCBI genetic codes (mitochondrial, etc.)
3. Codon usage analysis functions
4. Support for frameshifts and programmed ribosomal slippage
5. Integration with variant calling pipelines
6. Batch operations on multiple genes

## Notes

- All functions use 1-based coordinates (standard in biology)
- Functions automatically handle strand orientation
- Comprehensive error messages for invalid inputs
- Follows micromicon package conventions and style
- Ready for CRAN submission (all checks pass)

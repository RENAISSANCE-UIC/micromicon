# micromicon 0.3.2

* Updated use_cases_gd_consequence_enrichment_parallel.R to use clusters on Windows

### Breaking Changes

#### Function Renaming: predict_mutations → predict_variants

* **BREAKING**: Renamed `predict_mutations()` to `predict_variants()`
* Updated all documentation, README, vignettes, and internal references
* Users must update their code: change `predict_mutations(gd)` to `predict_variants(gd)`


# micromicon 0.3.2

## Architectural Focus

* Made `predict_variants()` and `compute_effects()` official S3 generics
* Reorganized code base to better conform with clean architecture design

# micromicon 0.3.0

## Major: Mutation Consequence Enrichment System

### New Function: pm_enrich_consequences()

* **NEW**: `pm_enrich_consequences()` enriches `predict_variants()` output with molecular consequences
* Supports all mutation types: SNP, DEL, INS, SUB
* Provides DNA sequences (reference and mutant), amino acid sequences, and consequence classification
* Consequence types: synonymous, missense, nonsense, frameshift, inframe_deletion, inframe_insertion, complex
* Works directly with `predict_variants()` output (handles arrows in gene names, position-based fallback)
* Includes QC tracking system with `qc_note` column for validation

### Production Hardening

* **Fixed**: Annotation parser now handles deletion ranges like "coding (152-278 / 435 nt)"
* **Added**: Stop codon truncation - mutant proteins terminate at first `*` (biological accuracy)
* **Added**: Gene resolver cleans arrows (←, →) and falls back to position-based lookup
* **Added**: Alternative start codon normalization (GTG/TTG/CTG → M)
* **Added**: QC note tracking for resolution issues and warnings

### Infrastructure

* **NEW**: Consolidated `%||%` operator into single exported definition in `R/operators.R`
* **Removed**: 10 duplicate `%||%` definitions across package files

# micromicon 0.2.9

## Breaking Changes

### Function Renaming: predicted_mutations → predict_mutations

* **BREAKING**: Renamed `predicted_mutations()` to `predict_mutations()` to use verb form
* **BREAKING**: Renamed `predicted_mutations_int()` to `predict_mutations_int()` (internal function)
* **BREAKING**: Renamed `predicted_mutations_orig()` to `predict_mutations_orig()` (internal function)
* Updated all documentation, README, and vignettes to use new function names
* Users must update their code: change `predicted_mutations(gd)` to `predict_mutations(gd)`

# micromicon 0.2.8

## Bug Fixes: Sequence Extraction Functions

### Fixed Missing Namespace Prefixes

* **Fixed**: Added proper namespace prefixes throughout `get_feature_fasta()` and `get_roi_fasta()`
  - `GenomicRanges::strand()`, `GenomicRanges::seqnames()`
  - `BiocGenerics::start()`, `BiocGenerics::end()`, `BiocGenerics::width()`
  - `GenomeInfoDb::seqinfo()`, `GenomeInfoDb::seqlengths()`, `GenomeInfoDb::seqlevels()`
  - `Biostrings::getSeq()`, `Biostrings::reverseComplement()`, `Biostrings::writeXStringSet()`
  - `S4Vectors::Rle()`, `IRanges::IRanges()`

### Enhanced FaFile vs DNAStringSet Support

* **Fixed**: `get_roi_fasta()` now handles both indexed FASTA files (FaFile) and in-memory sequences (DNAStringSet)
* **Fixed**: Sequence name preservation when using `width()` on DNAStringSet
* **Fixed**: Proper sequence extraction using `subseq()` for DNAStringSet (since `getSeq()` only works with FaFile)
* **Fixed**: Automatic seqname matching when single-sequence genomes have mismatched names

### Documentation Updates

* **Added**: Comprehensive README updates introducing two complementary modes:
  - **Genome Navigation Mode** (`genome_entity`) for reference exploration
  - **Variation Analysis Mode** (`genome_entity_gd`) for mutation tracking
* **Added**: New vignette `variation-analysis.Rmd` covering breseq integration and `genome_entity_gd` workflow
* **Updated**: `gene-utilities.Rmd` to replace "Mode A" terminology with current workflow
* **Updated**: README examples to show both modes with complete function calls

### Breaking Changes

None - all changes are backward compatible bug fixes and documentation improvements.

# micromicon 0.2.7

## Major: Genome Difference (GD) Parser and Viewing API

### New Type-Aware GD Parser

* **Added**: Complete rewrite of `parse_gd_annotated()` with type-aware dispatch pattern
* Handles all mutation types: SNP, DEL, INS, SUB, MOB, AMP, CON, INV
* Handles all evidence types: RA (Read Alignment), MC (Missing Coverage), JC (Junction), UN (Unknown)
* Handles validation types: TSEQ, PFLP, RFLP, PFGE, PHYL, CURA, NOTE, FPOS
* Events automatically binned by kind in `gd$provenance$by_kind` for efficient filtering
* Flexible parsing with offset detection handles malformed breseq output gracefully
* Comprehensive bounds checking with configurable strict mode

### New Viewing and Analysis Functions

* **`gd_events_table()`** - Convert nested events to tidy data frame with optional tag expansion
* **`view_events()`** - Preview events with sensible column ordering and filtering
* **`view_mutations()`** - Convenience wrapper for mutation events only
* **`view_evidence()`** - Convenience wrapper for evidence events with coverage columns
* **`probe_ra()`** - Diagnostic tool showing raw vs parsed RA event fields
* **`peek_tags()`** - Inspect tags for a specific event in key-value format
* **`tag_get_first()`** - Extract first value of a tag into a new column
* **`tag_get_concat()`** - Concatenate multi-valued tags with separator
* **`print_event()`** - Pretty-print individual events with type-specific formatting
* **`format_freq()`** - Format numeric frequencies as percentages

### Tag Extraction API

Events can be extracted with tags as a list-column for programmatic manipulation:

```r
ev_tbl <- gd_events_table(gd, expand_tags = FALSE)
ev_tbl <- ev_tbl |>
  tag_get_first("gene_name") |>
  tag_get_concat("locus_tag", sep = ";")
```

### Documentation

* Added comprehensive roxygen2 documentation for all new functions
* All user-facing functions now have parameter descriptions, return values, and examples
* Improved `parse_gd_annotated()` documentation with detailed description of type-aware parsing

### Bug Fixes

* **Fixed**: Missing field defaults for AMP, CON, and INV mutation types in parser
* **Fixed**: Removed duplicate `tag_get()` function
* **Fixed**: Added missing `@export` tags to `summary.genome_entity_gd()`

### Technical Details

* Parser uses `.parse_mut()` and `.parse_evidence()` dispatch functions
* `.rectify_for_constructor()` ensures backward compatibility with existing code
* All events include `kind` field ("mutation", "evidence", "validation")
* Events indexed by kind in `provenance$by_kind` for O(1) filtering
* Robust handling of missing tags and optional fields

---

# micromicon 0.2.6

## Added convenience functions to be used in genome diff integration

---

# micromicon 0.2.4

## BLAST Gateway Refactoring 

* Added gateway for local blast `gateways_blast.R`
* Refactored `blast_protein()` to use gateway

---

# micromicon 0.2.3

## Bug Fixes

### search_features() Now Searches ID, Name, and Alias Fields

* **Fixed**: `search_features()` now correctly finds features when searching by their `ID` or `Name` attributes
* Previously only searched `gene`, `locus_tag`, and `product` fields, missing common GFF3 annotation fields
* Now searches all commonly used GFF3 attributes: `ID`, `Name`, `Alias`, `gene`, `locus_tag`, `product`
* This fix makes `search_features()` compatible with various annotation sources (breseq, Prokka, NCBI)
* Updated documentation to clarify which fields are searched
* Added comprehensive test suite (6 new tests) covering all searchable fields

**Technical Details:**
- Modified `R/generics_queries.R:24` (S3 method implementation)
- Updated `R/controllers_query_controller.R:175` and documentation for consistency
- All tests pass (67 total tests)

---

# micromicon 0.2.2

* Docs: removed pipe workflow claim

---

# micromicon 0.2.0

## Major Changes

### One-Way Format Conversion Enforcement

* **BREAKING**: Formalized that format conversion is ONE-WAY ONLY (GenBank → GFF3+FASTA)
* Added runtime warnings when exporting GenBank-sourced data to GFF3/FASTA formats
* Users are now clearly warned about metadata loss (organism, taxonomy, references, comments)
* No `write_genbank()` function will be implemented - this is intentional, not a missing feature

### New Features

* Added source tracking to `genome_entity` objects via `import_source` attribute
* Export warnings can be suppressed with `options(micromicon.warn_export = FALSE)`
* Added `get_example_file()` helper function for accessing package example data

### Documentation

* Added comprehensive "Format Conversion Rules" section to README
* Updated function documentation for `write_gff3()` and `write_fasta()` with metadata loss warnings
* Updated `read_genome()` documentation with format differences and conversion rules
* Added format conversion rules section to vignette
* Clarified in source comments that GenBank export is forbidden (not just unimplemented)

### Testing

* Added comprehensive test suite for export warnings (7 new tests)
* All tests pass (52 total tests)

### Bug Fixes

* Fixed GenBank metadata parsing with comprehensive Bioconductor error messages
* Fixed `analyze_gene()` compatibility with GenBank files
* Fixed NULL check in `extract_sequences_by_name()`
* Added auto-conversion of genome_entity to legacy format for backward compatibility

---

# micromicon 0.1.0

* Initial release
* Support for reading GenBank and GFF3+FASTA formats
* Unified `genome_entity` object representation
* Feature search and sequence extraction
* Export to GFF3 and FASTA formats
* Clean architecture design with separated layers
* Optional Bioconductor integration

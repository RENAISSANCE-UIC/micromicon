

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
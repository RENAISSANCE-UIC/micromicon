# TODO List for micromicon

## High Priority

### NAMESPACE Cleanup
- [ ] Review and audit all exported functions
- [ ] Consider making these internal (using `@keywords internal` instead of `@export`):
  - `execute_import_genbank()` - line 53
  - `execute_import_gff_fasta()` - line 54
  - `execute_build_indices()` - line 52
  - `validate_*()` functions (lines 106-112 in NAMESPACE)
    - `validate_dna_sequence()`
    - `validate_feature_type()`
    - `validate_genome_entity()`
    - `validate_genomic_coordinates()`
    - `validate_locus_tag()`
    - `validate_sequence_topology()`
    - `validate_strand()`
  - Converter functions that might be internal:
    - `entity_to_gbk_list()` - line 50
    - `entity_to_legacy_genome_obj()` - line 51
    - `gbk_to_entity()` - line 66
    - `gff_fasta_to_entity()` - line 79

**Decision needed**: Are these meant for advanced users/extensibility, or should they be internal?

## Medium Priority

### Documentation
- [ ] Review all function documentation for consistency
- [ ] Ensure all `@param` tags match actual function signatures
- [ ] Fix remaining roxygen2 warnings from `devtools::check_man()`
- [ ] Add more examples to vignettes

### Testing
- [ ] Verify tests handle both "gene" and "CDS" feature types appropriately
- [ ] Add tests for case-insensitive type matching in `search_features()`
- [ ] Test GFF3 from various sources (breseq, Prokka, NCBI)

## Low Priority

### Code Quality
- [ ] Review use of `@keywords internal` throughout codebase
- [ ] Standardize error messages using cli package
- [ ] Consider deprecation strategy for legacy functions
- [ ] Review and consolidate duplicate functionality (if any)

### Performance
- [ ] Profile common workflows
- [ ] Optimize feature searching for large genomes
- [ ] Consider caching strategies for repeated queries

## Future Features

### Enhancements
- [ ] VCF format support
- [ ] Functional consequence prediction for mutations
- [ ] Multi-genome comparative analysis
- [ ] Time-series mutation tracking
- [ ] Integration with phylogenetic tools

### Infrastructure
- [ ] S4/S7 migration planning (see README roadmap)
- [ ] CI/CD setup for automated testing
- [ ] CRAN submission preparation
- [ ] Benchmarking suite

## Notes

### Recent Fixes (v0.2.8)
- ✅ Fixed namespace prefixes in sequence extraction functions
- ✅ Added DNAStringSet vs FaFile handling
- ✅ Made `search_features()` type parameter case-insensitive
- ✅ Documented feature type variations by annotation source
- ✅ Made `search_features_legacy_internal()` truly internal

### Design Decisions
- Keep backward compatibility where possible
- Prefer clear error messages over silent failures
- Document assumptions about data sources (breseq, Prokka, etc.)
- Maintain clean architecture separation

---

**Last updated**: 2026-02-12

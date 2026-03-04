# TODO List for micromicon

## High Priority

### NAMESPACE Cleanup
- [x] Review and audit all exported functions
- [x] Make non-generics internal (`@keywords internal`) — kept only S3 generics + `read_genome()` + `read_variants()`
- [x] Make `execute_import_genbank()`, `execute_import_gff_fasta()`, `execute_build_indices()` internal
- [x] Make all `validate_*()` functions internal
- [x] Make converter functions internal (`entity_to_gbk_list()`, `entity_to_legacy_genome_obj()`, `gbk_to_entity()`, `gff_fasta_to_entity()`)

## Medium Priority

### Documentation
- [x] Review all function documentation for consistency
- [x] Ensure all `@param` tags match actual function signatures
- [x] Fix remaining roxygen2 warnings from `devtools::check_man()`
- [x] Add more examples to vignettes

### Testing
- [x] Verify tests handle both "gene" and "CDS" feature types appropriately
- [ ] Add tests for case-insensitive type matching in `search_features()`
- [ ] Test GFF3 from various sources (breseq, Prokka, NCBI)

## Low Priority

### Code Quality
- [x] Review use of `@keywords internal` throughout codebase
- [x] Standardize error messages using cli package
- [x] Consider deprecation strategy for legacy functions
- [x] Review and consolidate duplicate functionality (if any)
- [x] Remove commented-out code blocks in `queries.R` (next code wrangling session)
- [ ] Investigate overlap between `generics_export.R` and `generics_export_methods.R`

### Performance
- [x] Profile common workflows
- [ ] Optimize feature searching for large genomes - Maybe use compiled code for some performance boosts? Rust?
- [x] Consider caching strategies for repeated queries

## Future Features

### Enhancements
- [x] `read_variants()` S3 generic — dispatches on `genome_entity`, currently supports `"gd"` (annotated breseq genome diff); extensible to `"vcf"`
- [ ] VCF format support (via `read_variants(..., format = "vcf")`)
- [x] Functional consequence prediction for mutations (`pm_enrich_consequences()` / `annotate_variants()`)
- [ ] Multi-genome comparative analysis - betaverse project
- [ ] Time-series mutation tracking - betaverse project
- [ ] Integration with phylogenetic tools - betaverse project

### Infrastructure
- [ ] S4/S7 migration planning (see README roadmap)
- [ ] CI/CD setup for automated testing
- [ ] CRAN submission preparation - maybe ROpenSci?
- [ ] Benchmarking suite

**Last updated**: 2026-03-03

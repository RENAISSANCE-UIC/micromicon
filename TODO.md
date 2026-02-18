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
- [x] Standardize error messages using cli package
- [ ] Consider deprecation strategy for legacy functions
- [ ] Review and consolidate duplicate functionality (if any)
- [ ] Remove commented-out code blocks in `queries.R` (next code wrangling session)
- [ ] Investigate overlap between `generics_export.R` and `generics_export_methods.R`

### Performance
- [ ] Profile common workflows
- [ ] Optimize feature searching for large genomes
- [ ] Consider caching strategies for repeated queries

## Future Features

### Enhancements
- [x] `read_variants()` S3 generic — dispatches on `genome_entity`, currently supports `"gd"` (annotated breseq genome diff); extensible to `"vcf"`
- [ ] VCF format support (via `read_variants(..., format = "vcf")`)
- [x] Functional consequence prediction for mutations (`pm_enrich_consequences()` / `compute_effects()`)
- [ ] Multi-genome comparative analysis
- [ ] Time-series mutation tracking
- [ ] Integration with phylogenetic tools

### Infrastructure
- [ ] S4/S7 migration planning (see README roadmap)
- [ ] CI/CD setup for automated testing
- [ ] CRAN submission preparation
- [ ] Benchmarking suite

## Notes

### Recent Work (v0.3.x, 2026-02)

**Architecture reorganization:**
- ✅ Migrated `queries.R` functions to correct layers (S3 methods → `generics_queries.R`, helpers → `controllers_query_controller.R`)
- ✅ Fixed alphabetical-load shadow bug (queries.R was overwriting earlier-loading generics/controllers)
- ✅ Restored `use_cases_gd_consequence_enrichment.R` (deleted in prior wrangling, restored from git)
- ✅ Added `compute_effects.genome_entity_gd` method

**Naming consistency:**
- ✅ `genome_metadata()` → `get_genome_metadata()`
- ✅ `genome_seqnames()` → `get_genome_seqnames()`
- ✅ `features()` → `get_features()`

**New generics:**
- ✅ `get_roi_features()` — exported S3 generic wrapping `analyze_roi()`
- ✅ `read_variants()` — exported S3 generic; `read_variants(ref, "annotated.gd")` replaces direct `parse_gd_annotated()` calls

**Bug fixes:**
- ✅ `predict_mutations()` duplicate RA entries for INS mutations (fixed `scalar_or` truncation + added deduplication in `gd_presenters.R`)
- ✅ JC parser negative `alignment_overlap` bug (added `.valid_int()` in `gd_parser.R`)

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
- Public API: 21 exports total — 19 S3 generics + `read_genome()` + `read_variants()`

---

**Last updated**: 2026-02-18

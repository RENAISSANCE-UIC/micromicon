# micromicon Architecture Dictionary

A reference for contributors and Future Claude: what code belongs where.

---

## Naming Convention

Every source file follows the pattern:

```
R/<layer>_<topic>.R
```

The layer prefix encodes *what kind of code lives there*, not the domain. The topic names the domain subject. Multiple topic segments are separated by underscores (e.g., `use_cases_import_import_genbank.R`).

**Exception:** Experimental and prototype code uses a prefix such as `beta_<name>.R` and is exempt from all layer rules below. Beta files may be messy, uncategorized, and unstable. They should not export symbols to NAMESPACE.

---

## Layers

### `controllers_`

**Purpose:** The user-facing API surface. Every function a user calls by name lives here or is dispatched through a generic defined in `generics_`.

**Rules:**
- Functions are exported (`@export`)
- Heavy roxygen documentation: `@param`, `@return`, `@examples`, `@seealso`
- Validates user-supplied arguments (type checks, path existence, etc.)
- Does **not** contain business logic — delegates immediately to a use case or generic
- Does **not** call gateways or parsers directly
- May call other controllers only for formatting/display

**Examples:**
- `controllers_query_controller.R` — `search_features()`, `extract_by_name()`, `gene_info()`
- `controllers_gd_views.R` — `view_mutations()`, `view_evidence()`, `peek_tags()`
- `controllers_gd_mutations.R` — `predict_mutations()`, `compute_effects()`
- `controllers_import_controller.R` — `read_genome()`, `init_genome()`

---

### `generics_`

**Purpose:** S3 generic definitions and their method dispatch. This is the polymorphism layer.

**Rules:**
- Define the generic with `@export` and document the *interface* (not the implementation)
- Define S3 methods (`function.class_name`) in the **same file** as the generic when the method count is small; split into a companion `generics_<topic>_methods.R` if methods grow large
- Method bodies may contain logic when the logic is trivially tied to dispatch; otherwise delegate to a use case
- No gateway or parser calls

**Examples:**
- `generics_queries.R` — `search_features.genome_entity()`, `get_gene_dna.genome_entity()`
- `generics_accessors.R` — `features()`, `sequences()`, `genome_metadata()`
- `generics_print.R` — `print.genome_entity()`, `print.genome_entity_gd()`
- `generics_summary.R` — `summary.genome_entity()`, `summary.genome_entity_gd()`
- `generics_gd.R` — GD-specific generics (`compute_effects`, `predict_mutations`)
- `generics_export.R` — generic definition; `generics_export_methods.R` — S3 methods

---

### `use_cases_`

**Purpose:** Business logic and orchestration. A use case takes validated inputs, calls gateways/parsers/entities to do work, and returns a result. It makes decisions; it does not do I/O directly.

**Rules:**
- Functions are **not** exported as a rule (use `@keywords internal`); exceptions exist for advanced users
- Accepts domain objects and primitive values; returns domain objects or tidy data
- May call: entities, gateways, frameworks parsers, framework utils
- Does **not** call controllers or other use cases in a chain (keep the call graph flat)
- One primary action per file; file name mirrors the action (e.g., `use_cases_import_import_genbank.R`)

**Examples:**
- `use_cases_gd_parse.R` — `parse_gd_annotated()` (orchestrates parser + entity construction)
- `use_cases_import_import_genbank.R` — `execute_import_genbank()`
- `use_cases_query_extract_sequences_by_name.R` — `extract_sequences_by_name()`
- `use_cases_gd_consequence_enrichment_fast.R` — consequence enrichment pipeline

---

### `entities_`

**Purpose:** Core domain data structures: constructors, field schemas, and structural validators.

**Rules:**
- Defines `new_<type>()` constructors and `validate_<type>()` validators
- No I/O, no business logic, no gateway calls
- May call framework utils for pure data manipulation
- Validators throw errors on structural inconsistencies (not semantic errors — those belong in use cases)
- Exported selectively: constructors for extension authors; validators may be exported or internal depending on user value

**Examples:**
- `entities_genome_entity.R` — `new_genome_entity()`, `validate_genome_entity()`
- `entities_genome_entity_gd.R` — `new_genome_entity_gd()`, `validate_genome_entity_gd()`
- `entities_validators.R` — standalone validation helpers (`validate_strand()`, `validate_genomic_coordinates()`, etc.)

---

### `gateways_`

**Purpose:** Integration surfaces to external formats, tools, or data sources. A gateway is the only code that touches the file system, external processes, or external APIs.

**Rules:**
- Returned via a factory function `create_<name>_gateway()` that returns a named list of closures
- Each closure does exactly one external operation and returns an R-native result
- No domain object construction inside the gateway — return lists, data frames, or raw character vectors
- No business logic; the use case decides what to do with the result
- File I/O (read/write), subprocess calls (`system2`), and HTTP requests all belong here

**Examples:**
- `gateways_genbank_gateway.R` — `create_genbank_gateway()`: reads GenBank files
- `gateways_fasta_gateway.R` — `create_fasta_gateway()`: reads/writes FASTA
- `gateways_blast.R` — `create_blast_gateway()`: runs BLASTP subprocesses
- `gateways_gd_annotation.R` — gateway for breseq annotation lookups
- `gateways_gd_reference.R` — gateway for reference manifest construction

---

### `frameworks_parsers_`

**Purpose:** Pure mechanical adapters. Turn raw text or bytes into R-native transport shapes (lists, data frames, named vectors). Zero domain decisions.

**Rules:**
- No entity construction — return plain lists or data frames, never `new_*()` objects
- No I/O — accept a file path or character vector as input, do not open connections themselves (or if they must, they are not responsible for what to do with failures)
- Deterministic: same input → same output, no side effects
- Functions are `@keywords internal`; not exported
- May call framework utils

**Examples:**
- `frameworks_parsers_gd_parser.R` — `parse_gd_raw()`: parses `.gd` text → `list(header, events, kind_vec)`
- `frameworks_parsers_genbank_parser.R` — GenBank file parser → raw record list

---

### `frameworks_utils_`

**Purpose:** Stateless utility functions with no domain knowledge. Shared toolbox for all other layers.

**Rules:**
- Functions are purely functional (no side effects, no I/O)
- No domain objects — works on primitive R types (character, numeric, list, data.frame)
- `@keywords internal`; not exported unless the utility has standalone user value (e.g., `format_freq`)
- Subtopic suffixes help organize: `_general`, `_gd`, `_gd_helpers`, `_string_utils`, `_file_utils`, `_helpers`

**Examples:**
- `frameworks_utils_general.R` — `first_non_na()`, `scalar_or()`, `format_freq()`
- `frameworks_utils_gd.R` — GD-specific pure utilities (hashing, annotation parsing helpers)
- `frameworks_utils_gd_helpers.R` — `is_annotated_gd()`, `canonical_event_hash()`
- `frameworks_utils_string_utils.R` — string manipulation utilities
- `frameworks_utils_helpers.R` — general-purpose R helpers

---

### `frameworks_` (other subtypes)

**Purpose:** Platform-level concerns that don't fit the parser or utils subtypes.

| Subtype | What goes here |
|---|---|
| `frameworks_operators.R` | Custom infix operators (`%\|\|%`, etc.) |
| `frameworks_bioconductor_compat.R` | Bioconductor version detection, conditional loading |
| `frameworks_utils_examples.R` | `get_example_file()` and test fixture helpers |
| `frameworks_utils_legacy_*.R` | Backward-compat shims for old APIs; to be pruned |

---



# Position Mapping Specification

## Overview

The `map_genomic_to_cds()` and `map_cds_to_genomic()` functions provide bijective mappings between genomic coordinates and CDS (coding DNA sequence) positions.

## Requirements

✓ **Respect gene strand** (+ or -)
✓ **Return 1-based positions**
✓ **Strictly monotonic in CDS space**
⚠️ **Respect exons** - NOT YET IMPLEMENTED (single-exon only)

## Formulas

### Plus Strand (+)

**Genomic → CDS:**
```
cds_pos = genomic_pos - gene_start + 1
```

**CDS → Genomic:**
```
genomic_pos = gene_start + cds_pos - 1
```

### Minus Strand (-)

**Genomic → CDS:**
```
cds_pos = gene_end - genomic_pos + 1
```

**CDS → Genomic:**
```
genomic_pos = gene_end - cds_pos + 1
```

## Key Concepts

### 1. CDS Position = Canonical 5'→3' Coding Orientation

- **CDS position 1** is always the first base of the coding sequence in the 5'→3' direction
- For **plus strand**: CDS position 1 is at `gene_start` (numerically smallest)
- For **minus strand**: CDS position 1 is at `gene_end` (numerically largest)

### 2. Strand Matters

**Plus strand gene (301-1434):**
```
Genomic:  301  302  303  ...  1434
CDS:        1    2    3  ...  1134
          (increasing together)
```

**Minus strand gene (100-200):**
```
Genomic:  200  199  198  ...  100
CDS:        1    2    3  ...  101
          (inverse relationship)
```

### 3. Monotonicity

- In **CDS space**, positions are always strictly increasing (1, 2, 3, ...)
- For **plus strand**: CDS increases with genomic coordinate
- For **minus strand**: CDS increases as genomic coordinate decreases
- This ensures CDS position always represents 5'→3' coding direction

### 4. Bijectivity

The mappings are perfectly bijective (one-to-one and onto):
```r
genomic → CDS → genomic  (identity)
CDS → genomic → CDS      (identity)
```

## Test Results

### Plus Strand (ampC: 301-1434, +)

```
Genomic → CDS:
  301 → 1      ✓
  302 → 2      ✓
  351 → 51     ✓
  1434 → 1134  ✓

CDS → Genomic:
  1 → 301      ✓
  2 → 302      ✓
  51 → 351     ✓
  1134 → 1434  ✓

Bijectivity: PASS ✓
Monotonicity: PASS ✓ (strictly increasing)
```

### Minus Strand (test gene: 100-200, -)

```
Genomic → CDS:
  200 → 1    ✓  (5' end)
  199 → 2    ✓
  150 → 51   ✓
  100 → 101  ✓  (3' end)

CDS → Genomic:
  1 → 200    ✓  (5' end)
  2 → 199    ✓
  51 → 150   ✓
  101 → 100  ✓  (3' end)

Bijectivity: PASS ✓
Monotonicity: PASS ✓ (strictly increasing in CDS space)
```

## Implementation Details

### Current Scope: Single-Exon Genes

The current implementation assumes **single-exon CDS** (typical for bacterial genes):
- Uses `feat$start` and `feat$end` directly
- Returns `NA` for positions outside gene boundaries
- Works for both strands

### Future: Multi-Exon Support

For eukaryotic genes with introns, the mapping would need to:
1. Parse exon coordinates from GFF3 (e.g., from `ranges` column)
2. Skip intron regions when calculating CDS position
3. Handle edge cases (positions in introns → NA)

Example for multi-exon gene:
```
Exon 1: 100-200 (101 bp)
Intron: 201-300 (skipped)
Exon 2: 301-400 (100 bp)

Genomic → CDS:
  100 → 1     (exon 1)
  200 → 101   (exon 1)
  250 → NA    (intron)
  301 → 102   (exon 2)
  400 → 201   (exon 2)
```

## Usage Examples

### Example 1: SNP in Plus Strand Gene

```r
entity <- read_genome("genome.gbk")

# ampC gene: 2607270-2608403 on + strand
# SNP at genomic position 2607320

cds_pos <- map_genomic_to_cds(entity, "ampC", 2607320)
# Returns: 51 (51st base of CDS)

codon_num <- ceiling(cds_pos / 3)
codon_pos <- ((cds_pos - 1) %% 3) + 1
# Codon 17, position 3
```

### Example 2: SNP in Minus Strand Gene

```r
# Gene on minus strand: 1000-2000
# SNP at genomic position 1950

cds_pos <- map_genomic_to_cds(entity, "minusGene", 1950)
# Returns: 51 (51st base in 5'→3' CDS direction)

# Inverse mapping
genomic_pos <- map_cds_to_genomic(entity, "minusGene", 51)
# Returns: 1950 (verified)
```

### Example 3: Find Position of Mutation

```r
# You have a CDS mutation at position 100
# What's the genomic coordinate?

genomic_pos <- map_cds_to_genomic(entity, "dnaA", 100)
# Returns the exact genomic coordinate

# Verify reference base
ref_base <- get_roi_dna(entity, info$chrom, genomic_pos, genomic_pos)
```

## Validation

### Boundary Testing
- ✓ Start position (CDS 1)
- ✓ End position (CDS length)
- ✓ Positions outside gene → NA
- ✓ Invalid CDS positions → NA

### Mathematical Properties
- ✓ Bijective (one-to-one)
- ✓ Monotonic (strictly increasing in CDS space)
- ✓ Symmetric (g→c→g = g, c→g→c = c)
- ✓ 1-based indexing throughout

## Error Handling

```r
# Position outside gene
map_genomic_to_cds(entity, "dnaA", 99999)
# Returns: NA

# Invalid CDS position
map_cds_to_genomic(entity, "dnaA", 99999)
# Returns: NA

# Position in intron (future multi-exon support)
map_genomic_to_cds(entity, "eukaryotic_gene", intron_pos)
# Will return: NA (position not in CDS)
```

## See Also

- `gene_info()` - Get gene coordinates and strand
- `validate_variant_in_gene()` - Validate variant at genomic position
- `get_gene_dna()` - Extract gene sequence (handles strand automatically)

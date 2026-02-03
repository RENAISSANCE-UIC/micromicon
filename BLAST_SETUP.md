# Setting Up Local BLAST for micromicon

The `blast_protein()` function requires a **local BLAST database**. Remote BLAST is not supported.

## Requirements

### 1. Install BLAST+ Command-Line Tools

**Option A: Download from NCBI**
```bash
# Download from: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
# Extract and add to PATH
```

**Option B: Install via Conda**
```bash
conda install -c bioconda blast
```

**Option C: Install via Package Manager**
```bash
# Ubuntu/Debian
sudo apt-get install ncbi-blast+

# macOS
brew install blast
```

### 2. Download a Local BLAST Database

**SwissProt (Recommended for testing - smaller size)**
```bash
# Create directory for databases
mkdir -p ~/blastdb
cd ~/blastdb

# Download SwissProt (curated, ~500K proteins, ~300MB compressed)
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz

# Extract
tar -xzf swissprot.tar.gz
```

**Other databases:**
- `nr` - Non-redundant protein database (very large, ~100GB)
- `refseq_protein` - RefSeq proteins
- `pdbaa` - PDB protein sequences

Browse available databases: ftp://ftp.ncbi.nlm.nih.gov/blast/db/

### 3. Set BLASTDB Environment Variable

**In R (temporary, current session only):**
```r
Sys.setenv(BLASTDB = "/path/to/blastdb")
```

**In .Renviron (persistent):**
```r
# Edit .Renviron file
usethis::edit_r_environ()

# Add this line:
BLASTDB=/path/to/blastdb
```

**In shell profile (system-wide):**
```bash
# Add to ~/.bashrc or ~/.zshrc
export BLASTDB="/path/to/blastdb"
```

## Usage

### Basic Example

```r
library(micromicon)

# Set BLASTDB if not already set
Sys.setenv(BLASTDB = "/path/to/blastdb")

# Read genome
genome <- read_genome("ecoli.gbk")

# Get protein sequence
acrB_cds <- search_features(genome, pattern = "acrB", type = "CDS")
acrB_protein <- acrB_cds$translation

# Run BLASTP against local SwissProt
results <- blast_protein(
  sequence = acrB_protein,
  database = "swissprot",
  evalue = 1e-5,
  max_hits = 20
)

# View results
print(results)
```

### Advanced Usage

```r
# Use different database
results <- blast_protein(
  sequence = protein_seq,
  database = "nr",  # requires local nr database
  evalue = 1e-10,
  threads = 8,
  max_hits = 50
)

# Provide explicit database directory (overrides BLASTDB)
results <- blast_protein(
  sequence = protein_seq,
  database = "swissprot",
  dbdir = "/custom/path/to/databases",
  evalue = 1e-5
)

# Filter results
high_quality <- reduce_hits(
  results,
  min_qcov = 80,      # minimum 80% query coverage
  min_pident = 50,    # minimum 50% identity
  besthit = TRUE      # only best hit per query
)
```

## Troubleshooting

### Error: "Could not locate blastp"
- BLAST+ is not installed or not in PATH
- Solution: Install BLAST+ (see above)

### Error: "Cannot find BLAST database files"
- Database not downloaded or extracted
- BLASTDB not set correctly
- Solution: Check that .psq, .phr, .pin files exist in BLASTDB directory

### Error: "blastp exited with status 2"
- Database path contains spaces or special characters
- Solution: Use BLASTDB environment variable (not explicit path)
- Or move database to path without spaces

### Warning: "Using relative database path without explicit dbdir"
- This is OK if BLASTDB is set
- The function will find the database via BLASTDB

## R Package Dependencies

The `blast_protein()` function requires these R packages:
- **cli** - for nice error messages (in Imports)
- **dplyr** - for data manipulation (in Imports)
- **readr** - for reading BLAST output (in Imports)
- **tibble** - for tidy data frames (in Imports)
- **processx** - for better subprocess handling (in Suggests)

All required packages are automatically installed with micromicon.

## Verifying Setup

Test your setup:

```r
# Check BLASTDB
Sys.getenv("BLASTDB")

# Check blastp binary
Sys.which("blastp")

# List available databases
list.files(Sys.getenv("BLASTDB"), pattern = "\\.psq$")

# Quick test
blast_protein("MTEYKLVVVG", database = "swissprot", max_hits = 5)
```

## Performance Tips

1. **Use more threads**: `threads = 8` or more for large databases
2. **Adjust E-value**: Higher values (1e-5) for distant homologs, lower (1e-10) for close matches
3. **Limit hits**: `max_hits = 10` for quick results
4. **SwissProt first**: Test with SwissProt before using larger databases like nr

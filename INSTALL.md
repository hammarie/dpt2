# DPT2: Multiplex dPCR Primer Design Tool

A Python tool for designing and scoring primer pairs for multiplex digital PCR (dPCR) applications.

## Installation

### Prerequisites

- Python 3.9 or higher
- pip (Python package manager)

### Option 1: Install from source (Recommended for development)

```bash
# Clone or download the repository
cd /path/to/dpt2

# Install in development mode
pip install -e .
```

### Option 2: Standard installation

```bash
cd /path/to/dpt2
pip install .
```

### Option 3: Install with requirements file

```bash
cd /path/to/dpt2
pip install -r requirements.txt
```

## BLAST Database Setup (Important!)

The tool can optionally use local BLAST screening for primer specificity checking. To enable this feature:

### 1. Download RefSeq RNA data

```bash
# Make a directory for BLAST databases
mkdir -p ~/blast_databases
cd ~/blast_databases

# Download human RNA sequences (RefSeq mRNA)
# For humans - replace with other organisms as needed
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Protein/human.rna.fna.gz
gunzip human.rna.fna.gz
```

### 2. Create BLAST database

```bash
# Install BLAST+ if not already installed
# macOS: brew install blast
# Ubuntu/Debian: sudo apt-get install ncbi-blast+
# Windows: Download from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

# Create the database
makeblastdb -in human.rna.fna -dbtype nucl -out human_rna
```

This creates files: `human_rna.ndb`, `human_rna.nhr`, `human_rna.nin`, etc.

## Usage

### Basic usage (without BLAST screening)

```bash
python scripts/design_and_score_primers.py --genes RPP30 MRGPRX1
```

### With BLAST specificity screening

```bash
python scripts/design_and_score_primers.py \
  --genes RPP30 MRGPRX1 \
  --blast-screen \
  --blast-db /path/to/blast/database/human_rna
```

### Custom parameters

```bash
python scripts/design_and_score_primers.py \
  --genes RPP30 MRGPRX1 \
  --num-return 10 \
  --product-min 70 \
  --product-max 400 \
  --out-prefix MyResults \
  --results-dir my_results
```

### Available arguments

```
--genes GENE1 GENE2 ...
    Gene symbols or NCBI gene IDs (required)
    Example: RPP30 MRGPRX1

--num-return N (default: 10)
    Maximum primer pairs per gene to design

--product-min LENGTH (default: 70)
    Minimum PCR product size in bp

--product-max LENGTH (default: 400)
    Maximum PCR product size in bp

--out-prefix NAME (default: PrimerAnalyse_local)
    Prefix for output files

--results-dir PATH (default: results)
    Directory to store numbered result files

--blast-screen
    Enable BLAST specificity screening

--blast-db PATH
    Path to BLAST database (required with --blast-screen)

--blast-max-mismatches N (default: 1)
    Maximum mismatches for target hits

--blast-min-total-mismatches N (default: 2)
    Minimum mismatches for off-target hits

--blast-min-3prime-mismatches N (default: 2)
    Minimum mismatches in 3' end (last 5 bp)
```

## Output

The tool generates two files per run in the `results/` folder:

### 1. Excel file (`1_PrimerAnalyse_local.xlsx`)

Contains three sheets:
- **primer_sequences**: Forward and reverse primer sequences
- **hairpin_self_dimer_primer3**: Hairpin and self-dimer ΔG values
- **heterodimer_dG_primer3**: Cross-dimer ΔG values between primers

### 2. CSV file (`1_PrimerAnalyse_local_primerblast.csv`)

BLAST screening results (when `--blast-screen` is used)

Files are numbered automatically: `1_`, `2_`, `3_`, etc. for each run.

## Troubleshooting

### BLAST database error

If you see: "BLAST Database error: Database memory map file error"

1. Verify the database files exist and are readable
2. Rebuild the database:
   ```bash
   makeblastdb -in source.fna -dbtype nucl -out database_name -overwrite
   ```

### Missing sequences

If genes are not found:
- Verify gene symbols are correct (case-sensitive NCBI symbols)
- Check internet connection (fetches from NCBI)
- Provide alternative gene ID or accession number

### Installation issues on Windows

If you encounter `blastn` not found:
1. Download BLAST+ from: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
2. Add the BLAST `bin` directory to your PATH environment variable
3. Restart terminal/command prompt

## Dependencies

- **primer3-py** (>=2.0.0): Primer design
- **biopython** (>=1.80): Sequence fetching from NCBI
- **pandas** (>=1.3.0): Data processing
- **pyyaml** (>=6.0): Configuration

### Optional

- **blastn** (from NCBI BLAST+ toolkit): For specificity screening

## Development

To contribute or modify the tool:

```bash
# Install in development mode
pip install -e .

# Modifications to source files in src/mpx_dpcr/ take effect immediately
```

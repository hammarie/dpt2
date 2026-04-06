# DPT2: Multiplex dPCR Primer Design Tool

A Python-based tool for automated design and thermodynamic analysis of primer pairs for multiplex digital PCR (dPCR) applications.

## Features

✨ **Primer Design**
- Automatic primer pair design using Primer3
- Customizable product size ranges
- Multiple primer candidates per gene
- Thermodynamic quality filtering

🔬 **Thermodynamic Analysis**
- Hairpin structure free energy (ΔG) calculations
- Self-dimer formation analysis
- Cross-dimer (heterodimer) scoring between primer pairs

🎯 **Specificity Screening** (Optional)
- Local BLAST screening against human RefSeq RNA
- 3' mismatch enforcement for specificity
- Customizable stringency parameters

📊 **Output**
- Excel files with all thermodynamic analysis
- CSV files with BLAST screening results

## Quick Start

### 1. Installation (One-time setup)

**macOS/Linux:**
```bash
cd /path/to/dpt2
chmod +x install.sh
./install.sh
```

**Windows:**
```bash
cd \path\to\dpt2
install.bat
```

Or manually:
```bash
pip install -e .
```

### 2. Run Primer Design

```bash
python scripts/design_and_score_primers.py --genes RPP30 MRGPRX1
```

Results appear in `results/` folder as `1_PrimerAnalyse_local.xlsx`

### 3. (Optional) Set up BLAST Specificity Screening

Follow the detailed instructions in [INSTALL.md](INSTALL.md#blast-database-setup-important)

Then run with BLAST enabled:
```bash
python scripts/design_and_score_primers.py \
  --genes RPP30 MRGPRX1 \
  --blast-screen \
  --blast-db /path/to/human_rna
```

## System Requirements

- **Python:** 3.9+
- **Operating System:** Windows, macOS, or Linux
- **Disk Space:** ~50MB (including BLAST database optional)
- **Internet:** Required to fetch gene sequences from NCBI (one-time per gene)

## Directory Structure

```
dpt2/
├── scripts/
│   └── design_and_score_primers.py    # Main script to run
├── src/mpx_dpcr/
│   ├── fetch.py                       # NCBI sequence fetching
│   ├── design.py                      # Primer3 design
│   ├── evaluate.py                    # Thermodynamic analysis
│   ├── blast_check.py                 # BLAST specificity screening
│   └── cli.py                         # Command-line interface
├── results/                           # Output folder (auto-created)
├── requirements.txt                   # Python dependencies
├── pyproject.toml                     # Project configuration
├── INSTALL.md                         # Detailed installation guide
└── README.md                          # This file
```

## Usage Examples

### Basic Primer Design

```bash
# Design primers for two genes
python scripts/design_and_score_primers.py --genes RPP30 MRGPRX1

# Specify different parameters
python scripts/design_and_score_primers.py \
  --genes RPP30 MRGPRX1 \
  --num-return 15 \
  --product-min 80 \
  --product-max 350
```

### With BLAST Screening

```bash
python scripts/design_and_score_primers.py \
  --genes RPP30 MRGPRX1 \
  --blast-screen \
  --blast-db /path/to/blast/database/human_rna \
  --blast-max-mismatches 3 \
  --blast-min-total-mismatches 2 \
  --blast-min-3prime-mismatches 2
```

### Adjust Output Directory

```bash
python scripts/design_and_score_primers.py \
  --genes RPP30 MRGPRX1 \
  --results-dir my_custom_results_folder \
  --out-prefix MyProject
```

Output files: `my_custom_results_folder/1_MyProject.xlsx`

## Output Explanation

### Excel Spreadsheet

**Sheet 1: primer_sequences**
- Gene ID, forward/reverse primer sequences, role

**Sheet 2: hairpin_self_dimer_primer3**
- Per-oligo thermodynamic metrics
- Hairpin ΔG and self-dimer ΔG values
- Negative values indicate stability

**Sheet 3: heterodimer_dG_primer3**
- Cross-dimer scores between all primer pairs
- Identifies potential unwanted primer-primer interactions


### CSV File (BLAST results)

When using `--blast-screen`:
- Gene ID and pair index
- Forward/reverse primer sequences
- BLAST hit summaries
- Specificity pass/fail status

## Parameters Reference

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--genes` | *required* | Gene symbols or NCBI IDs |
| `--num-return` | 10 | Primer pairs per gene |
| `--product-min` | 70 | Min fragment size (bp) |
| `--product-max` | 400 | Max fragment size (bp) |
| `--out-prefix` | PrimerAnalyse_local | Output filename prefix |
| `--results-dir` | results | Output directory |
| `--blast-screen` | OFF | Enable BLAST specificity check |
| `--blast-db` | none | Path to BLAST database |
| `--blast-max-mismatches` | 1 | Mismatches for target hits |
| `--blast-min-total-mismatches` | 2 | Mismatches for off-targets |
| `--blast-min-3prime-mismatches` | 2 | 3' end mismatches required |
| `--blast-3prime-window` | 5 | 3' end window size (bp) |
| `--blast-ignore-mismatches-ge` | 6 | Ignore targets with ≥N mismatches |

## Troubleshooting

### "BLAST Database error"
- Verify database path is correct
- Rebuild BLAST database (see INSTALL.md)

### "No mRNA found for gene"
- Check gene symbol spelling
- Try NCBI gene ID instead of symbol
- Verify internet connection

### "No primers passed filters"
- Relax constraints: increase `--product-max`, decrease `--max-self-end`
- Try without `--blast-screen` to isolate issue
- Check sequence quality

### Windows: "blastn not found"
- Download BLAST+ from: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
- Add BLAST bin directory to PATH
- Restart terminal

## File Naming Convention

Output files are numbered automatically:
- **First run:** `1_PrimerAnalyse_local.xlsx`, `1_PrimerAnalyse_local_primerblast.csv`
- **Second run:** `2_PrimerAnalyse_local.xlsx`, `2_PrimerAnalyse_local_primerblast.csv`
- And so on...

This prevents overwriting previous results.

## Dependencies

- **primer3-py** – Primer design engine
- **biopython** – NCBI sequence fetching
- **pandas** – Data handling and Excel export
- **pyyaml** – Configuration (optional)
- **BLAST+** – Optional local sequence screening

## Development & Contributing

The source code is organized as a Python package:

```bash
# Installation for development
pip install -e .

# Changes to src/mpx_dpcr/ files take effect immediately
```

Core modules:
- `mpx_dpcr.fetch` – Sequence retrieval from NCBI
- `mpx_dpcr.design` – Primer3 wrapper
- `mpx_dpcr.evaluate` – Thermodynamic calculations
- `mpx_dpcr.blast_check` – BLAST specificity screening

## Performance

- **Typical runtime:** 1-5 minutes per gene (depends on BLAST screening)
- **Memory usage:** ~500 MB (includes BLAST database)
- **Network:** 2-5 MB total (gene sequence downloads)

## Known Limitations

1. NCBI gene fetching requires active internet connection
2. Primer3 single-threaded (not parallelized)
3. Maximum 7 genes typically recommended for multiplex applications

---

**Last Updated:** April 2026
**Version:** 0.1.0

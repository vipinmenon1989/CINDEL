# CINDEL — CRISPR-Cpf1 sgRNA Activity Scoring Tool

**Author:** Vipin Menon, BIG Lab
**Original:** May 2016 | **Modernized:** 2026
**Language:** Python 3.8+

---

## Overview

CINDEL calculates predicted activity scores for **CRISPR-Cpf1 (Cas12a)** guide RNAs. It helps researchers identify the most active guides from a pool of candidates.

The scoring model combines:
- **RNA free energy** of the 23-bp target region (via ViennaRNA)
- **Global nucleotide composition** (A, AC, CG, CC, TA counts)
- **Position-dependent nucleotide features** (50 coefficients for specific positions)

The final score is a logistic-transformed linear combination of these features, returning a value between **0** (inactive) and **1** (highly active).

---

## Sequence Requirements

| Property        | Requirement                          |
|-----------------|--------------------------------------|
| Total length    | **27 bp**                            |
| PAM (pos 0–3)   | **TTTV** where V = A, G, or C        |
| Target (pos 4–27) | 23 bp                              |

Valid PAM sequences: `TTTA`, `TTTG`, `TTTC`

---

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/vipinmenon1989/CINDEL.git
cd CINDEL
```

### 2. Create a virtual environment (recommended)
```bash
python3 -m venv venv
source venv/bin/activate        # macOS / Linux
venv\Scripts\activate           # Windows
```

### 3. Install dependencies
```bash
pip install -r requirements.txt
```

> **Note:** ViennaRNA must be installed separately via conda (recommended):
> ```bash
> conda install -c bioconda viennarna
> ```

---

## Usage

### Mode A — Batch scoring from CSV

Score all guide RNAs in a CSV file. The input must have a header row with columns `sequence_id` and `sequence`.

```bash
python CINDEL.py -a input.csv
python CINDEL.py -a input.csv --output results.csv
```

**Input CSV format:**
```
sequence_id,sequence
guide_001,TTTACGTAGCTAGCTAGCTAGCTAGCT
guide_002,TTTGCGATCGATCGATCGATCGATCGA
```

**Output CSV format:**
```
sequence_id,sequence,score
guide_001,TTTACGTAGCTAGCTAGCTAGCTAGCT,0.734521
guide_002,TTTGCGATCGATCGATCGATCGATCGA,0.412038
```

---

### Mode B — Single sequence scoring

Score one 27-bp guide RNA directly on the command line.

```bash
python CINDEL.py -b TTTACGTAGCTAGCTAGCTAGCTAGCT
```

**Output:**
```
Sequence          : TTTACGTAGCTAGCTAGCTAGCTAGCT
PAM               : TTTA
Target (23 bp)    : CGTAGCTAGCTAGCTAGCTAGCT
CINDEL Score      : 0.734521
```

---

### Mode C — Guide RNA finder

Scan a long genomic sequence for all valid Cpf1 guide RNA candidates (must be > 100 bp).

```bash
python CINDEL.py -c ATCGATCGATCGATCG...LONG_SEQUENCE...
```

**Output:** A ranked table of all candidate guides sorted by score (highest first).

---

## Command-Line Options

| Option              | Description                                        |
|---------------------|----------------------------------------------------|
| `-a`, `--batch`     | Batch mode: path to input CSV file                 |
| `-b`, `--single`    | Single mode: one 27-bp sequence string             |
| `-c`, `--find`      | Finder mode: long genomic sequence to scan         |
| `-o`, `--output`    | Output CSV path for batch mode (default: Score.csv)|
| `-v`, `--verbose`   | Enable verbose/debug logging                       |
| `-h`, `--help`      | Show help message and exit                         |

---

## Dependencies

| Package    | Purpose                              | Install                          |
|------------|--------------------------------------|----------------------------------|
| ViennaRNA  | RNA secondary structure free energy  | `conda install -c bioconda viennarna` |
| Python ≥ 3.8 | Runtime                            | [python.org](https://www.python.org) |

See `requirements.txt` for the full list.

---

---

## Snakemake Workflow

Run the scoring pipeline via Snakemake instead of calling CINDEL.py directly:

```bash
pip install snakemake
snakemake -n                       # dry run - shows planned jobs
snakemake --cores 1                # real run using config/config.yaml
snakemake --cores 1 --use-conda    # real run, auto-creating the conda env in envs/environment.yaml
```

Edit `config/config.yaml` to switch between `mode: batch` (score a CSV) and
`mode: single` (score one sequence), and to point at your own input file.
A small demo dataset is provided at `demo/input.csv`.

---

## CI/CD

Every push and pull request to `main`/`develop` runs three GitHub Actions jobs:
1. **syntax-check** - compiles CINDEL.py and lints it with pyflakes
2. **snakemake-dry-run** - validates the workflow graph (`snakemake -n`)
3. **smoke-test** - runs the full pipeline on the demo data, checks the output, and runs the unit tests in `tests/`

---

## Project Structure

```
CINDEL/
├── CINDEL.py                    # Main scoring tool
├── Snakefile                    # Snakemake workflow (batch / single modes)
├── config/config.yaml           # Workflow configuration
├── envs/environment.yaml        # Conda environment (Python + ViennaRNA)
├── demo/input.csv               # Small demo dataset
├── tests/test_cindel.py         # Unit tests (pytest)
├── .github/workflows/ci.yml     # CI: syntax check, dry-run, smoke test
├── requirements.txt              # Python dependencies
├── README.md                    # This file
└── .gitignore                   # Git ignore rules
```

---


---

## Known Issues & Recent Fixes (2026)

The 2025 "modernization" commit updated this README to describe an
argparse-based CLI, but the underlying `CINDEL.py` code had not actually
been changed and was still Python 2 (unparenthesized `print`, `xrange`,
mixed tabs/spaces indentation, binary-mode CSV writes). It has now been
rewritten to match: real Python 3, `argparse`-based CLI, type hints, and
docstrings.

A real logic bug was also found and fixed: the batch-mode PAM check
```python
if set[j][0:3] == 'TTTA' or 'TTTG' or 'TTTC':
```
always evaluated to `True` (non-empty strings are truthy in a boolean
`or`), so invalid PAMs were never actually rejected, and the slice was
3 bp instead of the required 4. This is now `is_valid_pam()`, which
correctly checks all 4 PAM bases against `TTTA`/`TTTG`/`TTTC`.

## Citation

If you use CINDEL in your research, please cite the original work from the BIG Lab.

---

## License

Copyright © BIG Lab. All rights reserved.

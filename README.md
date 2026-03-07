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

## Project Structure

```
CINDEL/
├── CINDEL.py          # Main scoring tool
├── requirements.txt   # Python dependencies
├── README.md          # This file
└── .gitignore         # Git ignore rules
```

---

## Citation

If you use CINDEL in your research, please cite the original work from the BIG Lab.

---

## License

Copyright © BIG Lab. All rights reserved.

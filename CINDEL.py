#!/usr/bin/env python3
"""
CINDEL - CRISPR-Cpf1 (Cas12a) sgRNA activity scoring tool.

Author: Vipin Menon, BIG Lab
Original: May 2016 | Modernized: 2026

The scoring model combines RNA free energy of the 23-bp target region,
global mono/dinucleotide composition, and position-dependent nucleotide
features into a logistic-transformed activity score in the range [0, 1].

Sequence requirement: total length 27 bp, where positions 0-3 are the
PAM (must be TTTV, V = A, G, or C) and positions 4-27 are the 23-bp
target region.
"""

import argparse
import csv
import logging
import math
import sys
from typing import List, Tuple

import RNA

logger = logging.getLogger("CINDEL")

VALID_PAMS: Tuple[str, ...] = ("TTTA", "TTTG", "TTTC")

# Position-dependent coefficients for single-nucleotide and dinucleotide
# features, as (position, motif, weight). Unchanged from the original 2016
# model - do not edit without re-deriving the fit.
POSITION_PARAMETERS: List[Tuple[int, str, float]] = [
    (9, 'A', 0.037608778), (13, 'AA', 0.167340348), (25, 'AA', 0.364227683),
    (6, 'AA', -0.37092459), (6, 'AC', 0.083159983), (9, 'AC', 0.459109261),
    (16, 'AG', -0.257994467), (22, 'AT', 0.469569195), (26, 'C', 0.166195809),
    (14, 'CA', 0.11154741), (5, 'CA', -0.317283912), (23, 'CC', 0.130179283),
    (4, 'CC', 1.065828073), (24, 'CG', -0.186173348), (6, 'CG', -0.231597251),
    (7, 'CG', 0.275263279), (13, 'CT', 0.025086518), (19, 'CT', 0.164219823),
    (24, 'CT', 0.083624311), (6, 'CT', 0.120799685), (7, 'CT', -0.336424167),
    (9, 'CT', -0.094751841), (19, 'G', -0.076228205), (22, 'G', -0.066751042),
    (23, 'G', -0.102884761), (26, 'G', -0.043253074), (4, 'G', 0.489354734),
    (10, 'GA', 0.059114363), (18, 'GA', 0.077528126), (21, 'GA', 0.216211044),
    (4, 'GA', 0.077064537), (20, 'GC', -0.123286737), (21, 'GG', -0.124324225),
    (22, 'GG', -0.358412408), (5, 'GG', -0.058789262), (8, 'GG', 0.126113307),
    (23, 'GT', -0.15567756), (21, 'T', 0.048550229), (4, 'T', -0.781196013),
    (9, 'T', -0.028800337), (23, 'TA', 0.444556704), (24, 'TA', 0.274164613),
    (7, 'TA', 0.391209275), (15, 'TC', 0.123925716), (25, 'TC', 0.194045073),
    (10, 'TG', -0.190097449), (16, 'TG', 0.036768788), (6, 'TG', 0.12530035),
    (10, 'TT', -0.173588498), (12, 'TT', -0.035141319), (16, 'TT', -0.127599558),
]

INTERCEPT = -1.830804463
FREE_ENERGY_COEF = 0.11018728

# Global mono/dinucleotide composition coefficients, keyed by motif.
NUCLEOTIDE_COMPOSITION_COEFS = {
    'A': 0.101900819,
    'AC': 0.038389043,
    'CG': 0.006835173,
    'CC': -0.24886262,
    'TA': 0.037868724,
}


def calculate_score(seq: str) -> float:
    """Compute the CINDEL activity score for a single 27-bp guide sequence.

    Args:
        seq: A 27-bp sequence: 4-bp PAM (TTTV) + 23-bp target region.

    Returns:
        A logistic-transformed activity score between 0 (inactive) and
        1 (highly active).

    Raises:
        ValueError: If seq is not exactly 27 bp long.
    """
    if len(seq) != 27:
        raise ValueError(f"Sequence must be exactly 27 bp, got {len(seq)} bp: {seq!r}")

    score = INTERCEPT

    # RNA free energy of the 23-bp target region (positions 4-27).
    target_region = seq[4:27]
    free_energy = round(RNA.fold(target_region)[-1], 0)
    score += free_energy * FREE_ENERGY_COEF

    # Global nucleotide composition.
    for motif, coef in NUCLEOTIDE_COMPOSITION_COEFS.items():
        score += coef * seq.count(motif)

    # Position-dependent nucleotide features.
    for position, motif, weight in POSITION_PARAMETERS:
        if seq[position:position + len(motif)] == motif:
            score += weight

    return 1.0 / (1.0 + math.exp(-score))


def is_valid_pam(seq: str) -> bool:
    """Return True if seq starts with a valid Cpf1 PAM (TTTA, TTTG, or TTTC).

    Note: the original 2016 code had a logic bug here
    (`if seq[0:3] == 'TTTA' or 'TTTG' or 'TTTC':`), which evaluated the
    `or` clauses as always-truthy non-empty strings, so PAM validation
    never actually rejected anything. This version fixes that.
    """
    return seq[0:4] in VALID_PAMS


def run_batch(input_csv: str, output_csv: str = "Score.csv") -> None:
    """Score every guide RNA in a CSV file.

    Input CSV must have a header row and columns: sequence_id, sequence.
    Output CSV has columns: sequence_id, sequence, score.
    """
    results = []
    skipped = 0
    with open(input_csv, newline="") as fh:
        reader = csv.reader(fh)
        next(reader, None)  # header
        for row in reader:
            if len(row) < 2:
                continue
            seq_id, seq = row[0].strip(), row[1].strip()
            if not is_valid_pam(seq):
                logger.warning(
                    "Skipping %s: PAM must be TTTV (V=A/G/C), got %r",
                    seq_id, seq[0:4],
                )
                skipped += 1
                continue
            results.append((seq_id, seq, calculate_score(seq)))

    with open(output_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sequence_id", "sequence", "score"])
        for seq_id, seq, score in results:
            writer.writerow([seq_id, seq, f"{score:.6f}"])

    logger.info(
        "Wrote %d scored guides to %s (%d skipped for invalid PAM)",
        len(results), output_csv, skipped,
    )


def run_single(sequence: str) -> None:
    """Score one 27-bp guide RNA and print a human-readable summary."""
    seq = sequence.strip()
    if not is_valid_pam(seq):
        raise ValueError(f"Invalid PAM {seq[0:4]!r}; must be TTTA, TTTG, or TTTC")
    score = calculate_score(seq)
    print(f"Sequence          : {seq}")
    print(f"PAM               : {seq[0:4]}")
    print(f"Target (23 bp)    : {seq[4:27]}")
    print(f"CINDEL Score      : {score:.6f}")


def run_finder(sequence: str) -> None:
    """Scan a long sequence for candidate Cpf1 guides, ranked by score."""
    sequence = sequence.strip()
    if len(sequence) <= 100:
        raise ValueError("Input sequence for finder mode must be > 100 bp")

    candidates = []
    search_from = 0
    while True:
        idx = sequence.find("TTT", search_from)
        if idx == -1:
            break
        if idx + 27 <= len(sequence):
            candidate = sequence[idx:idx + 27]
            if is_valid_pam(candidate):
                score = calculate_score(candidate)
                candidates.append((candidate, candidate[0:4], candidate[4:27], score))
        search_from = idx + 1

    candidates.sort(key=lambda c: c[3], reverse=True)

    print(f"{'sequence':<30}{'PAM':<6}{'target':<25}score")
    for seq, pam, target, score in candidates:
        print(f"{seq:<30}{pam:<6}{target:<25}{score:.6f}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="CINDEL.py",
        description="CINDEL: CRISPR-Cpf1 sgRNA activity scoring tool.",
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("-a", "--batch", metavar="INPUT_CSV",
                       help="Batch mode: path to input CSV file")
    mode.add_argument("-b", "--single", metavar="SEQUENCE",
                       help="Single mode: one 27-bp sequence string")
    mode.add_argument("-c", "--find", metavar="SEQUENCE",
                       help="Finder mode: long genomic sequence to scan")
    parser.add_argument("-o", "--output", default="Score.csv",
                         help="Output CSV path for batch mode (default: Score.csv)")
    parser.add_argument("-v", "--verbose", action="store_true",
                         help="Enable verbose/debug logging")
    return parser


def main(argv: List[str] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    if args.batch:
        run_batch(args.batch, args.output)
    elif args.single:
        run_single(args.single)
    elif args.find:
        run_finder(args.find)


if __name__ == "__main__":
    main(sys.argv[1:])

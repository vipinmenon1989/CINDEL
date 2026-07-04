"""Unit tests for CINDEL.py."""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import CINDEL


def test_is_valid_pam_accepts_valid():
    assert CINDEL.is_valid_pam("TTTACGTAGCTAGCTAGCTAGCTAGCT")
    assert CINDEL.is_valid_pam("TTTGCGATCGATCGATCGATCGATCGA")
    assert CINDEL.is_valid_pam("TTTCGGATCGATCGATCGATCGATCGA")


def test_is_valid_pam_rejects_invalid():
    # Regression test: the original code's PAM check
    # (`x == 'TTTA' or 'TTTG' or 'TTTC'`) always evaluated True.
    assert not CINDEL.is_valid_pam("AAAACGTAGCTAGCTAGCTAGCTAGCT")
    assert not CINDEL.is_valid_pam("GGGGCGATCGATCGATCGATCGATCGA")


def test_calculate_score_rejects_wrong_length():
    try:
        CINDEL.calculate_score("TOOSHORT")
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_calculate_score_returns_valid_probability():
    score = CINDEL.calculate_score("TTTACGTAGCTAGCTAGCTAGCTAGCT")
    assert 0.0 <= score <= 1.0
    assert not math.isnan(score)


def test_run_batch_skips_invalid_pam_and_scores_rest(tmp_path):
    input_csv = tmp_path / "input.csv"
    input_csv.write_text(
        "sequence_id,sequence\n"
        "guide_001,TTTACGTAGCTAGCTAGCTAGCTAGCT\n"
        "guide_002,AAAACGTAGCTAGCTAGCTAGCTAGCT\n"
    )
    output_csv = tmp_path / "out.csv"
    CINDEL.run_batch(str(input_csv), str(output_csv))

    content = output_csv.read_text()
    assert "guide_001" in content
    assert "guide_002" not in content

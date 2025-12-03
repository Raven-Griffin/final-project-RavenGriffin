# tests/test_pd_analyze_and_plot.py
import os
import pandas as pd
import pytest

from ratez.pd import analyze, plot_rates


def test_analyze_basic_rates(simple_alignment_dict):
    """
    Check that analyze() computes correct mismatch/indel rates
    for a tiny toy alignment.
    """
    data = simple_alignment_dict
    ref = "s1"

    df, avg_mismatch, avg_indel = analyze(data, ref)

    # Positions 1..4
    assert list(df["Position"]) == [1, 2, 3, 4]

    # Expected counts:
    # s1: A C G T
    # s2: A C - T  -> 1 indel at pos 3
    # s3: A G G T  -> 1 mismatch at pos 2
    # comparisons = 2 (s2, s3 vs s1)
    expected_mismatch = [0.0, 0.5, 0.0, 0.0]
    expected_indel = [0.0, 0.0, 0.5, 0.0]

    assert pytest_approx_list(df["MismatchRate"], expected_mismatch)
    assert pytest_approx_list(df["IndelRate"], expected_indel)

    # avg_mismatch = 1 mismatch / (4 positions * 2 seqs) = 0.125
    # avg_indel    = 1 indel    / (4 positions * 2 seqs) = 0.125
    assert avg_mismatch == pytest.approx(0.125)
    assert avg_indel == pytest.approx(0.125)


def pytest_approx_list(series, expected, tol=1e-6):
    """
    Helper to compare a pandas Series to a list with approximate floats.
    """
    values = list(series)
    assert len(values) == len(expected)
    for v, e in zip(values, expected):
        assert v == pytest.approx(e, rel=0, abs=tol)
    return True


def test_plot_rates_saves_png(simple_alignment_dict, tmp_path):
    """
    plot_rates() should write a PNG file to the given path.
    """
    data = simple_alignment_dict
    ref = "s1"
    df, _, _ = analyze(data, ref)

    out_path = tmp_path / "mutation_rates.png"
    plot_rates(df, out_path=str(out_path))

    assert out_path.exists()
    assert out_path.stat().st_size > 0

# tests/test_read_aln_and_ratescalc.py
import os
import pandas as pd

from ratez.bio import read_aln
from ratez.pd import analyze
from main import RatesCalc


def test_read_aln_returns_dict(tmp_alignment_fasta):
    """
    read_aln() should return a dict mapping IDs to aligned strings.
    """
    aln_dict = read_aln(str(tmp_alignment_fasta))

    assert isinstance(aln_dict, dict)
    assert set(aln_dict.keys()) == {"s1", "s2", "s3"}
    assert all(isinstance(v, str) for v in aln_dict.values())
    assert len(aln_dict["s1"]) == len(aln_dict["s2"]) == len(aln_dict["s3"])


def test_ratescalc_uses_read_aln(tmp_alignment_fasta):
    """
    RatesCalc.read() + .calc() should populate .data and .rates.
    """
    calc = RatesCalc(aln_file=str(tmp_alignment_fasta), reference="s1")
    calc.read()
    assert isinstance(calc.data, dict)
    assert "s1" in calc.data

    calc.calc()
    # After calc(), analyze() should have been called
    assert hasattr(calc, "rates")
    assert hasattr(calc, "avg_mismatch")
    assert hasattr(calc, "avg_indel")

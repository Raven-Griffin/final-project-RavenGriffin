# tests/conftest.py
import pytest
import sys
from pathlib import Path


# Add project root to PYTHONPATH
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
@pytest.fixture
def simple_alignment_dict():
    """
    Small toy alignment for testing analyze().
    s1 is the reference.
    """
    return {
        "s1": "ACGT",
        "s2": "AC-T",
        "s3": "AGGT",
    }

@pytest.fixture
def tmp_fasta_file(tmp_path: Path):
    """
    Tiny FASTA file with two short DNA sequences that contain ORFs.
    """
    content = """>s1
ATGAAAGGGTAA
>s2
ATGCCCCCCCTAA
"""
    fasta_path = tmp_path / "test.fasta"
    fasta_path.write_text(content)
    return fasta_path

@pytest.fixture
def tmp_alignment_fasta(tmp_path: Path):
    """
    Tiny aligned FASTA file (all sequences same length).
    """
    content = """>s1
ACGT
>s2
AC-T
>s3
AGGT
"""
    fasta_path = tmp_path / "aligned.fasta"
    fasta_path.write_text(content)
    return fasta_path

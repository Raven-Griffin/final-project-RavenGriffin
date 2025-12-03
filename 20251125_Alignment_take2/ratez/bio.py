#!/usr/bin/env python3
"""
bio.py

Utility function to read a multiple sequence alignment in FASTA format
using Biopython’s AlignIO module.

This function returns a dictionary mapping sequence IDs to their
aligned sequences (as uppercase strings with gaps represented by '-').

Example
-------
    from ratez.bio import read_aln

    alignment_dict = read_aln("alignment.fasta")
    for sid, seq in alignment_dict.items():
        print(sid, seq)
"""

from Bio import AlignIO
from typing import Dict


def read_aln(fname: str = "alignment.fasta") -> Dict[str, str]:
    """
    Read a FASTA alignment file and return sequences as a dictionary.

    Parameters
    ----------
    fname : str, optional
        Path to the input FASTA alignment file (default: "alignment.fasta").

    Returns
    -------
    dict[str, str]
        Dictionary where keys are sequence identifiers and values are
        uppercase aligned sequences as strings.

    Notes
    -----
    - Sequences are read using Biopython’s AlignIO interface.
    - All whitespace is stripped and bases are converted to uppercase.
    - The function assumes the alignment format is 'fasta'.
    """
    # Read the alignment using Biopython (returns an Alignment object)
    alignment = AlignIO.read(fname, "fasta")

    # Convert to a simple Python dictionary: {ID: SEQUENCE}
    dic = {
        str(record.id).strip(): str(record.seq).strip().upper()
        for record in alignment
    }

    return dic
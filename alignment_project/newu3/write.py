#!/usr/bin/env python3
"""
Utility to write sequences and their IDs into a FASTA-formatted file.
"""

import sys
from typing import List

def write_fa(ids: List[str], sequences: List[str], filepath: str) -> None:
	"""
	Write sequence IDs and sequences into a file in FASTA format.

	Parameters
	----------
	ids : List[str]
		List of sequence identifiers (e.g., ["Seq1", "Seq2"]).
	sequences : List[str]
		List of corresponding sequences (e.g., ["ATGCC", "GATTACA"]).
	filepath : str
		Path to the output FASTA file.

	Returns
	-------
	None

	Notes
	-----
	- Each sequence is written as:
		  >ID
		  SEQUENCE
	- Sequences are automatically uppercased and stripped of whitespace.
	- The function assumes `len(ids) == len(sequences)`.
	"""
	if len(ids) != len(sequences):
		raise ValueError("IDs and sequences must have the same length.")
	
	sys.stdout.write("Alignment:\n")
	with open(filepath, "w", encoding="utf-8") as f:
		for seq_id, seq in zip(ids, sequences):
			clean_seq = seq.strip().upper()
			f.write(f">{seq_id}\n{clean_seq}\n")
			sys.stdout.write(f">{seq_id}\n{clean_seq}\n")
	sys.stdout.write(f"Alignment saved into {filepath}\n//\n")
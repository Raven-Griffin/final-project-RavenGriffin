#!/usr/bin/env python3
"""
readFastA.py

Utility functions for reading FASTA files using regular expressions only.
This script is part of the newu2 module.
"""

import re
from typing import List


def read_fa(filepath: str) -> str:
	"""
	Read sequences from a FASTA file and return them as a formatted string.

	Parameters
	----------
	filepath : str
		Path to the FASTA file.

	Returns
	-------
	str
		A formatted string where each sequence entry is represented by two lines:
		- Line 1: Sequence ID (no description, only the ID).
		- Line 2: Sequence (single line, uppercase).

	Notes
	-----
	- Uses only regular expressions for parsing.
	- Sequence lines are joined and uppercased.
	- Example output for two sequences:

		>ID1
		ATGCGTAC
		>ID2
		GATTACA
	"""
	with open(filepath, "r", encoding="utf-8") as f:
		content = f.read()

	# Regex to capture FASTA entries: header + sequence
	# Header starts with ">", take first word as ID
	entries = re.findall(r">(.*?)\n([^>]*)", content, flags=re.DOTALL|re.MULTILINE)

	formatted_entries: List[str] = []
	for header, seq in entries:
		# Extract only sequence ID (first token before any space)
		seq_id_match = re.match(r"^(\S+)", header.strip())
		seq_id = seq_id_match.group(1) if seq_id_match else "UNKNOWN"

		# Clean sequence: remove whitespace, make uppercase
		sequence = re.sub(r"\s+", "", seq).upper()

		formatted_entries.append(f">{seq_id}\n{sequence}")

	return "\n".join(formatted_entries)
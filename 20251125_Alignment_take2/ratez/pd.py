#!/usr/bin/env python3
"""
pd.py

Compute per-site and overall mutation (mismatch) and indel rates
from a nucleotide alignment using pandas.

Each sequence in the input dictionary is compared to a reference sequence.
For every position:
  - If both residues are gaps ('-'), the site is skipped.
  - If exactly one residue is a gap, it counts as an indel.
  - If both residues are nucleotides but different, it counts as a mismatch.
The output includes per-position rates and overall averages.

Example
-------
	from ratez.pd import analyze
	data = {
		's1': 'GATACAGATACAG-AGATACA',
		's2': '--GATACAGATTA-AGATACA',
		's3': '--GATACACTCTA-AGATACA'
	}
	result_df, avg_mismatch, avg_indel = analyze(data, 's1')
"""

import sys
import pandas as pd
from typing import Dict, Tuple


def analyze(data: Dict[str, str], reference: str) -> Tuple[pd.DataFrame, float, float]:
	"""
	Analyze an alignment and compute mismatch and indel rates versus a reference sequence.

	Parameters
	----------
	data : dict[str, str]
		Dictionary mapping sequence IDs to aligned nucleotide strings.
	reference : str
		Key of the reference sequence to which all others are compared.

	Returns
	-------
	tuple
		A 3-element tuple containing:
		1. pd.DataFrame
			A table with columns:
			  - 'Position' (1-based index)
			  - 'MismatchRate' (fraction of mismatches per site)
			  - 'IndelRate' (fraction of indels per site)
		2. float
			Average mismatch rate across all positions.
		3. float
			Average indel rate across all positions.

	Notes
	-----
	- Positions with both residues equal to '-' are ignored (not counted).
	- Rates are normalized by the number of sequences compared to the reference.
	"""
	ref_seq = data[reference]
	seq_ids = [sid for sid in data if sid != reference]

	# Convert to a DataFrame for structure inspection (optional)
	df = pd.DataFrame({sid: list(seq) for sid, seq in data.items()})
	n_positions = len(ref_seq)

	# Initialize counters for mismatches and indels
	mismatch_counts = [0] * n_positions
	indel_counts = [0] * n_positions
	comparisons = 0  # number of non-reference sequences

	# Compare each non-reference sequence to the reference
	for sid in seq_ids:
		comparisons += 1
		for i, (ref_base, base) in enumerate(zip(ref_seq, data[sid])):
			# Skip positions where both are gaps
			if ref_base == '-' and base == '-':
				continue
			# Count indels (one gap, one nucleotide)
			elif ref_base == '-' or base == '-':
				indel_counts[i] += 1
			# Count mismatches (two nucleotides, different)
			elif ref_base != base:
				mismatch_counts[i] += 1

	# Compute per-position rates (normalized by number of comparisons)
	mismatch_rates = [m / comparisons for m in mismatch_counts]
	indel_rates = [g / comparisons for g in indel_counts]

	# Combine into a summary DataFrame
	result_df = pd.DataFrame({
		'Position': range(1, n_positions + 1),
		'MismatchRate': mismatch_rates,
		'IndelRate': indel_rates
	})

	# Compute overall averages (across all sites and sequences)
	avg_mismatch = sum(mismatch_counts) / (n_positions * comparisons)
	avg_indel = sum(indel_counts) / (n_positions * comparisons)

	# Print results to stdout for convenience
	sys.stdout.write(f"""Mutation rates per position:
{result_df.to_string(index=False)}
//
Average mismatch rate per position: {avg_mismatch:.4f}
Average indel rate per position: {avg_indel:.4f}
""")

	return result_df, avg_mismatch, avg_indel
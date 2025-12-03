#!/usr/bin/env python3
"""
main.py

Driver script for demonstrating the Needleman–Wunsch global alignment algorithm.
It relies on helper modules in the `newu3` and `ratez` packages to handle
sequence I/O, scoring, matrix filling, traceback, and mutation-rate analysis.

Example
-------
    python3 main.py -i sequences.fasta -o alignment.fasta -r s1
"""

##
# MODULES
##

import argparse
import sys
import re
from typing import Any, Tuple, Dict, List, Iterable
from collections import defaultdict

import numpy as np
from Bio import SeqIO

from newu3.read import read_fa
from newu3.write import write_fa
from newu3.num import parse_fa, init_mat, fill_matrix, trace_matrix
from ratez.bio import read_aln
from ratez.pd import analyze
from kmers.kmers import extract_kmers
from kmers.kmers import build_tree_dict_from_kmers


##
# AUXILIARY FUNCTIONS
##

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the program."""
    parser = argparse.ArgumentParser(
        description="Needleman–Wunsch global sequence alignment demonstration."
    )
    parser.add_argument(
        "-i", "--input",
        default="sequences.fasta",
        help="Path to input FASTA file."
    )
    parser.add_argument(
        "-o", "--output",
        default="alignment.fasta",
        help="Path to output FASTA alignment file."
    )
    parser.add_argument(
        "-r", "--reference",
        default="s1",
        help="Sequence ID to use as reference for mutation-rate analysis."
    )
    parser.add_argument(
        "-ma", "--match",
        default=1,
        type=int,
        help="Score for a match."
    )
    parser.add_argument(
        "-mi", "--mismatch",
        default=-1,
        type=int,
        help="Penalty for a mismatch."
    )
    parser.add_argument(
        "-in", "--indel",
        default=-1,
        type=int,
        help="Penalty for an insertion/deletion."
    )
    parser.add_argument(
        "-t", "--table",
        default=1,
        type=int,
        help="Translation table."
    )
    parser.add_argument(
        "-k", "--kmer",
        default=3,
        type=int,
        help="K-mer size."
    )
    return parser.parse_args()


def read_nucs(fasta_file: str) -> Dict:
    """
    Read a FASTA file and return a dictionary mapping sequence IDs to sequences.

    Parameters
    ----------
    fasta_file : str
        Path to a FASTA file.

    Returns
    -------
    dict[str, Bio.Seq.Seq]
        A dictionary where keys are FASTA record IDs (strings) and values
        are Biopython Seq objects representing nucleotide sequences.
    """
    dic: Dict[str, Any] = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seqid: str = str(record.id)
        seq = record.seq
        dic[seqid] = seq
    return dic


def find_lorfs(fasta_dic: Dict, ttable: int) -> Dict[str, np.ndarray]:
    """
    For each sequence in a dictionary, identify all ORFs (open reading frames)
    in forward, reverse, and reverse-complement orientations across all three
    reading frames, then return the longest ORF as an array of sets.

    Parameters
    ----------
    fasta_dic : dict[str, Bio.Seq.Seq]
        Dictionary of sequence IDs to nucleotide sequences.
    ttable : int
        NCBI translation table number to use when translating sequences.

    Returns
    -------
    dict[str, np.ndarray]
        Dictionary mapping sequence IDs to their longest ORF represented as an
        array of sets (one set of residues per column).
    """
    lorfs: Dict[str, np.ndarray] = {}

    for seqid, raw_seq in fasta_dic.items():
        orfs_found: List[str] = []

        # Three orientations: forward, reverse, reverse complement
        fwd = raw_seq
        rev = raw_seq[::-1]
        rco = raw_seq.reverse_complement()

        for my_seq in (fwd, rev, rco):
            for frame in (0, 1, 2):
                this_frame = my_seq[frame:]
                remainder = len(this_frame) % 3
                if remainder != 0:
                    this_frame = this_frame[:-remainder]

                aa = str(
                    this_frame.translate(
                        table=ttable,
                        to_stop=False,
                        cds=False,
                    )
                )

                orfs = [
                    orf for orf in re.compile(r"(M.+?\*)").findall(aa)
                    if len(orf) > 3
                ]
                orfs_found.extend(orfs)

        longest: str = max(orfs_found, key=len) if orfs_found else ""
        if longest:
            col_sets = [set(ch) for ch in longest]
            lorfs[seqid] = np.array(col_sets, dtype=object)

    return lorfs


def build_tree_dict(terminals: List[str]) -> Tuple[Dict[int, str], Dict[int, List[Any]]]:
    """
    Build a simple right-branching binary tree as a dictionary.

    Each new node joins the previous internal node with the next terminal,
    forming a progressively nested structure.

    Example
    -------
    >>> build_tree_dict(["t1", "t2", "t3", "t4"])
    (
        {1: 't1', 2: 't2', 3: 't3', 4: 't4'},
        {5: ['t1', 't2'], 6: [5, 't3'], 7: [6, 't4']}
    )
    """
    terminals_dic: Dict[int, str] = {}
    n = len(terminals)

    # Assign IDs to terminal nodes
    for i, name in enumerate(terminals, start=1):
        terminals_dic[i] = name

    # Create internal nodes by sequentially joining terminals
    internal_nodes: Dict[int, List[Any]] = {}
    current_node = n + 1
    left, right = 1, 2
    internal_nodes[current_node] = [terminals_dic[left], terminals_dic[right]]

    for j in range(3, n + 1):
        current_node += 1
        label1 = terminals_dic.get(current_node - 1, current_node - 1)
        label2 = terminals_dic.get(j, j)
        internal_nodes[current_node] = [label1, label2]

    return terminals_dic, internal_nodes


##
# CLASS
##

class AlignmentRunner:
    """
    Manages sequence alignment tasks using the Needleman–Wunsch algorithm.

    This class loads sequences from a FASTA file, builds a simple progressive
    tree, performs pairwise alignments following that topology, and writes
    the resulting multiple alignment to disk.
    """

    def __init__(
        self,
        fasta_file: str = "sequences.fasta",
        aln_file: str = "alignment.fasta",
        match: int = 1,
        mismatch: int = -1,
        indel: int = -1,
        ttable: int = 1,
        ksize: int = 3,
    ) -> None:
        """
        Initialize the AlignmentRunner.

        Parameters
        ----------
        fasta_file : str
            Path to the input FASTA file containing sequences to align.
        aln_file : str
            Path to the output FASTA file where the final alignment is written.
        match : int
            Score (or reward) for a match between residues.
        mismatch : int
            Penalty for a mismatch between residues.
        indel : int
            Penalty for an insertion or deletion (gap event).
        ttable : int
            Number of the translation table.
        ksize : int
            K-mer size.
        """
        # Load and parse sequences
        nuc_dic = read_nucs(fasta_file)
        lorfs = find_lorfs(nuc_dic, ttable)
        if len(lorfs) <= 1:
            sys.stderr.write("Cannot continue: don't have enough LORFs.\n")
            sys.exit(1)

        self.sequences: Dict[Any, np.ndarray] = lorfs
        self.kmers = extract_kmers(self.sequences, ksize)
        self.terminals, self.nodes = build_tree_dict_from_kmers(self.kmers)

        # Initialize empty placeholders for internal nodes
        for node in self.nodes:
            self.sequences[node] = np.array([], dtype=object)

        # Scoring scheme
        self.match = match
        self.mismatch = mismatch
        self.indel = indel

        # Output path
        self.aln_file = aln_file

    def matrix_initialization(self, seq1: np.ndarray, seq2: np.ndarray) -> None:
        """Initialize the dynamic programming matrix for a pair of sequences."""
        self.seq1 = seq1
        self.seq2 = seq2
        self.mat = init_mat(seq1, seq2)

    def fill(self) -> None:
        """Fill the dynamic programming matrix with scores."""
        self.mat = fill_matrix(
            matrix=self.mat,
            seq1_array=self.seq1,
            seq2_array=self.seq2,
            match=self.match,
            mismatch=self.mismatch,
            indel=self.indel,
        )

    def trace(self) -> None:
        """Perform traceback to reconstruct the optimal alignment."""
        self.aln = trace_matrix(
            matrix=self.mat,
            seq1_array=self.seq1,
            seq2_array=self.seq2,
        )

    def update_sequences(self, id1: Any, id2: Any) -> None:
        """Store the aligned sequences for two terminals after traceback."""
        self.sequences[id1] = self.aln[0]
        self.sequences[id2] = self.aln[1]

    def create_consensus(self, id1: Any, id2: Any, node: int) -> None:
        """
        Build a consensus sequence (set of residues per site)
        from two aligned sequences and assign it to an internal node.
        """
        seq1, seq2 = self.aln
        consensus = [seq1[i].union(seq2[i]) for i in range(len(seq1))]
        self.consensus = np.array(consensus, dtype=object)
        self.sequences[node] = self.consensus

    def patch(self, history: List[Any]) -> None:
        """
        Propagate newly introduced gaps to all previously aligned sequences.

        Parameters
        ----------
        history : list
            Ordered list of sequence IDs that have already been aligned.
        """
        # Positions of gaps in the first sequence of the new alignment
        gap_positions = [i for i, col in enumerate(self.aln[0]) if col == {'-'}]

        # Insert corresponding gaps into all earlier sequences
        for header in history:
            seq_arr = self.sequences[header]
            seq_list = list(seq_arr)

            # Insert from right to left so indices don't shift
            for pos in sorted(gap_positions, reverse=True):
                if pos <= len(seq_list):
                    seq_list.insert(pos, {'-'})
                else:
                    # If somehow position is past end, pad with gaps
                    seq_list.extend([{'-'}] * (pos - len(seq_list) + 1))

            self.sequences[header] = np.array(seq_list, dtype=object)

	def write(self) -> None:
    """Write the final multiple alignment to a FASTA file."""

    # --- FIX: Normalize sequence lengths ---
    max_len = max(len(seq) for seq in self.sequences.values())
    for key, seq in self.sequences.items():
        if len(seq) < max_len:
            pad = max_len - len(seq)
            seq = list(seq) + [{' - '}] * pad
            self.sequences[key] = np.array(seq, dtype=object)

    # Now write FASTA
    seqids = [self.terminals[key] for key in self.terminals]
    seqs = [
        ''.join(next(iter(s)) for s in self.sequences[seqid])
        for seqid in seqids
    ]
    write_fa(ids=seqids, sequences=seqs, filepath=self.aln_file)



class RatesCalc:
    """
    Wrapper for computing mutation and indel rates from an alignment file.
    """

    def __init__(self, aln_file: str = "alignment.fasta", reference: str = "s1") -> None:
        self.fasta_file = aln_file
        self.reference = reference

    def read(self) -> None:
        """Read a FASTA alignment into a dictionary using ratez.bio.read_aln."""
        self.data = read_aln(self.fasta_file)

    def calc(self) -> None:
        """
        Compute per-site and average mismatch/indel rates
        using ratez.pd.analyze().
        """
        self.rates, self.avg_mismatch, self.avg_indel = analyze(
            self.data,
            self.reference,
        )


def main() -> None:
    """
    Execute the complete alignment workflow.

    Steps
    -----
    1. Parse command-line arguments.
    2. Build progressive tree and perform successive pairwise alignments.
    3. Write the final multiple alignment to file.
    4. Compute and print mutation/indel rates versus a reference sequence.
    """
    args = parse_args()

    runner = AlignmentRunner(
        fasta_file=args.input,
        aln_file=args.output,
        match=args.match,
        mismatch=args.mismatch,
        indel=args.indel,
        ttable=args.table,
        ksize=args.kmer,
    )

    history: List[Any] = []

    # Sequentially align pairs according to the progressive tree
    for internal_node in runner.nodes:
        head1, head2 = runner.nodes[internal_node]
        seq1, seq2 = runner.sequences[head1], runner.sequences[head2]

        runner.matrix_initialization(seq1, seq2)
        runner.fill()
        runner.trace()
        runner.update_sequences(head1, head2)
        runner.create_consensus(head1, head2, internal_node)
        runner.patch(history)

        # Track all sequences that have been processed
        for name in (head1, head2):
            if name not in history:
                history.append(name)

    # Write final alignment to disk
    runner.write()

    # Calculate mutation and indel rates relative to the reference
    calc = RatesCalc(aln_file=args.output, reference=args.reference)
    calc.read()
    calc.calc()


if __name__ == "__main__":
    main()

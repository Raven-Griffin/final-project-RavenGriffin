# tests/test_newu3_num.py
import numpy as np
from Bio.Seq import Seq

from newu3.num import init_mat, fill_matrix, trace_matrix


def test_init_mat_shape():
    """
    init_mat() should create a matrix of correct shape for two sequences.
    """
    seq1 = np.array([{"A"}, {"C"}, {"G"}])
    seq2 = np.array([{"A"}, {"G"}])

    mat = init_mat(seq1, seq2)
    # +1 for the initial gap row/column
    assert mat.shape[0] == len(seq1) 
    assert mat.shape[1] == len(seq2) + 2
    print("Matrix shape is:", mat.shape)



def test_fill_and_trace_basic_alignment():
    """
    Simple sanity check for fill_matrix + trace_matrix:
    Align 'ACG' vs 'AG' with match/mismatch/indel scoring.
    """
    seq1 = np.array([{"A"}, {"C"}, {"G"}])
    seq2 = np.array([{"A"}, {"G"}])

    mat = init_mat(seq1, seq2)
    mat = fill_matrix(matrix=mat,
                      seq1_array=seq1,
                      seq2_array=seq2,
                      match=1,
                      mismatch=-1,
                      indel=-1)

    aln1, aln2 = trace_matrix(matrix=mat,
                              seq1_array=seq1,
                              seq2_array=seq2)

    # Check same length and composed of sets
    assert len(aln1) == len(aln2)
    assert all(isinstance(x, set) for x in aln1)
    assert all(isinstance(x, set) for x in aln2)

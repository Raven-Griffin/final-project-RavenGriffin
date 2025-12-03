import numpy as np
from Needlemen import init_mat, fill_matrix, trace_matrix

def test_matrix_initialization():
    s1 = np.array([{"A"}, {"C"}])
    s2 = np.array([{"A"}, {"T"}])
    mat = init_mat(s1, s2)
    assert mat.shape == (3, 3)

def test_fill_and_trace():
    s1 = np.array([{"A"}, {"C"}])
    s2 = np.array([{"A"}, {"C"}])

    mat = init_mat(s1, s2)
    mat = fill_matrix(mat, s1, s2, match=1, mismatch=-1, indel=-1)

    aln1, aln2 = trace_matrix(mat, s1, s2)

    assert "".join(next(iter(x)) for x in aln1) == "AC"
    assert "".join(next(iter(x)) for x in aln2) == "AC"

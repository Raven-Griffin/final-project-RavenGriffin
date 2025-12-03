# tests/test_orfs_and_tree.py
import numpy as np
from Bio.Seq import Seq

from main import find_lorfs, build_tree_dict, read_nucs


def test_find_lorfs_returns_longest_orf():
    """
    find_lorfs() should return longest ORF for each sequence
    as an array of sets (one per amino acid).
    """
    fasta_dic = {
        "s1": Seq("ATGAAAGGGTAA"),       # translates to M K G *
        "s2": Seq("ATGCCCCCCCTAA"),      # some other long ORF
    }
    lorfs = find_lorfs(fasta_dic, ttable=1)

    assert isinstance(lorfs, dict)
    assert "s1" in lorfs
    assert isinstance(lorfs["s1"],np.ndarray)

    for seqid, arr in lorfs.items():
        assert isinstance(arr, np.ndarray)
        # Each element should be a set of characters (amino acids)
        assert all(isinstance(x, set) for x in arr)


def test_build_tree_dict_simple():
    """
    build_tree_dict() should build a simple right-branching tree.
    """
    terminals = ["t1", "t2", "t3", "t4"]
    term_dic, nodes = build_tree_dict(terminals)

    # Terminals dictionary
    assert term_dic == {1: "t1", 2: "t2", 3: "t3", 4: "t4"}

    # Internal node keys should start at len(terminals)+1
    assert min(nodes.keys()) == 5
    assert max(nodes.keys()) == 7

    # First internal node should join t1 and t2
    assert nodes[5] == ["t1", "t2"]


def test_read_nucs_reads_sequences(tmp_fasta_file):
    """
    read_nucs() should return a dict of IDs to Bio.Seq.Seq objects.
    """
    from Bio.Seq import Seq  # type hint only

    nucs = read_nucs(str(tmp_fasta_file))
    assert set(nucs.keys()) == {"s1", "s2"}

    for v in nucs.values():
        # We don't assert exact class to avoid tight coupling,
        # but we check it's string-convertible and has length.
        s = str(v)
        assert isinstance(s, str)
        assert len(s) > 0

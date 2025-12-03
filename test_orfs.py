from main import find_lorfs
from Bio.Seq import Seq

def test_find_lorfs_basic():
    fasta_dic = {"s1": Seq("ATGGGGTAA")}
    lorfs = find_lorfs(fasta_dic, 1)
    assert "s1" in lorfs
    # ORF should translate to: M G G *
    assert len(lorfs["s1"]) == 3

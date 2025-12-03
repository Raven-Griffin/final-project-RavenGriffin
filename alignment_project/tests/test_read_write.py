import os
from newu3.read import read_fa
from newu3.write import write_fa

def test_write_and_read_fasta(tmp_path):
    out = tmp_path / "test.fa"
    ids = ["s1", "s2"]
    seqs = ["ACGT", "TTTT"]

    write_fa(ids, seqs, out)

    data = read_fa(out)

    assert data["s1"] == "ACGT"
    assert data["s2"] == "TTTT"

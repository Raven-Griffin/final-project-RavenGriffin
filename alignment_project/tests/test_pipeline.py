import os
from main import AlignmentRunner
from ratez.pd import mutation_rates
from ratez.bio import read_aln

def test_full_pipeline(tmp_path):
    # Create minimal FASTA
    fasta = tmp_path / "mini.fa"
    fasta.write_text(">s1\nATGGGGTAA\n>s2\nATGAGGTAA\n")

    out = tmp_path / "align.fa"

    runner = AlignmentRunner(
        fasta_file=str(fasta),
        aln_file=str(out),
        match=1, mismatch=-1, indel=-1,
        ttable=1, ksize=3
    )

    history = []
    for node in runner.nodes:
        h1, h2 = runner.nodes[node]
        seq1, seq2 = runner.sequences[h1], runner.sequences[h2]

        runner.matrix_initialization(seq1, seq2)
        runner.fill()
        runner.trace()
        runner.update_sequences(h1, h2)
        runner.create_consensus(h1, h2, node)
        runner.patch(history)

        for h in (h1, h2):
            if h not in history:
                history.append(h)

    runner.write()

    # Now check mutation rates
    aln = read_aln(str(out))
    rates = mutation_rates(aln, "s1")
    assert len(rates) >= 1

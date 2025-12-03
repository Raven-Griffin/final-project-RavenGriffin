from ratez.pd import mutation_rates

def test_mut_rates():
    aln = {
        "s1": "ACGT",
        "s2": "AC-T",
        "s3": "TCGT",
    }

    results = mutation_rates(aln, "s1")

    r2 = next(r for r in results if r["seq_id"] == "s2")
    assert r2["mismatches"] == 0
    assert r2["gaps"] == 1
    assert r2["mutation_rate"] == 1 / 4

    r3 = next(r for r in results if r["seq_id"] == "s3")
    assert r3["mismatches"] == 1

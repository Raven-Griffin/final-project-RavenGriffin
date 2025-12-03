from kmers.kmers import extract_kmers, build_tree_dict_from_kmers

def test_kmer_extraction():
    sequences = {
        "s1": [{"A"}, {"C"}, {"G"}, {"T"}],
        "s2": [{"A"}, {"A"}, {"G"}, {"T"}],
    }
    km = extract_kmers(sequences, 2)
    assert "s1" in km
    assert "AC" in km["s1"]

def test_tree_building():
    kmers = {"s1": {"AA": 2}, "s2": {"AC": 1}}
    terminals, nodes = build_tree_dict_from_kmers(kmers)
    assert len(terminals) == 2
    assert len(nodes) == 1

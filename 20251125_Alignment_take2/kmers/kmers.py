from typing import Any, Tuple, Dict, List, Iterable
from collections import defaultdict
import numpy as np

def extract_kmers(
    sequences: Dict[str, np.ndarray], 
    k: int = 3
) -> Dict[str, List[str]]:
    """
    Extract k-mers from a dictionary of amino acid sequences.

    Parameters
    ----------
    sequences : Dict[str, np.ndarray]
        A dictionary where each key is a sequence name (e.g., 's1')
        and each value is a NumPy array of single-character sets,
        e.g., array([{'M'}, {'N'}, {'S'}, ...], dtype=object).

    k : int, optional
        Length of the k-mers to extract (default is 3).

    Returns
    -------
    Dict[str, List[str]]
        A dictionary with the same keys, where each value is a list
        of k-mer strings extracted from the sequence.

    Notes
    -----
    - Residues are extracted by converting each set (e.g., {'M'}) to
      its single element ('M'), preserving order.
    - The last incomplete k-mer is discarded.
    """
    kmers_dict: Dict[str, List[str]] = {}

    for key, arr in sequences.items():
        # Convert array of sets into a list of characters
        seq_list = [next(iter(x)) for x in arr]
        seq_str = "".join(seq_list)

        # Extract sliding-window k-mers
        kmers = [seq_str[i:i+k] for i in range(len(seq_str) - k + 1)]
        kmers_dict[key] = kmers

    return kmers_dict

def jaccard_similarity(a: Iterable[str], b: Iterable[str]) -> float:
    """
    Compute the Jaccard similarity between two collections of k-mers.

    Parameters
    ----------
    a, b : Iterable[str]
        Collections (e.g., list or set) of k-mer strings.

    Returns
    -------
    float
        Jaccard similarity: |A ∩ B| / |A ∪ B|. Returns 0.0 if both are empty.
    """
    set_a = set(a)
    set_b = set(b)
    if not set_a and not set_b:
        return 0.0
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    return intersection / union

def order_by_kmer_similarity(kmer_dict: Dict[str, List[str]]) -> List[str]:
    """
    Order sequence IDs so that the most similar pairs (by k-mer Jaccard)
    are grouped first, and additional sequences are added in order of
    highest similarity to any already-chosen sequence.

    Parameters
    ----------
    kmer_dict : Dict[str, List[str]]
        Maps sequence IDs (e.g., 's1', 's2', ...) to lists of k-mers.

    Returns
    -------
    List[str]
        A list of sequence IDs in an order suitable for building a
        right-branching guide tree where more similar sequences
        tend to be joined earlier.
    """
    names = list(kmer_dict.keys())
    n = len(names)
    
    if n <= 1:
        return names
    
    # Precompute similarities between all pairs
    sim: Dict[Tuple[str, str], float] = {}
    for i in range(n):
        for j in range(i + 1, n):
            a, b = names[i], names[j]
            sim[(a, b)] = jaccard_similarity(kmer_dict[a], kmer_dict[b])
            
    # 1) Start with the pair of sequences that has the highest similarity
    best_pair = None
    best_sim = -1.0
    for (a, b), value in sim.items():
        if value > best_sim:
            best_sim = value
            best_pair = (a, b)
            
    if best_pair is None:
        # All similarities are zero or we only have one sequence
        return names
    
    ordered: List[str] = [best_pair[0], best_pair[1]]
    remaining = set(names) - set(ordered)
    
    # 2) Iteratively add the sequence that is most similar to any sequence
    #    already in 'ordered'.
    while remaining:
        best_candidate = None
        best_candidate_sim = -1.0
        
        for candidate in remaining:
            # Similarity of this candidate to the current cluster:
            # take the maximum similarity to any member of 'ordered'.
            max_sim_to_cluster = 0.0
            for existing in ordered:
                pair = (existing, candidate) if (existing, candidate) in sim else (candidate, existing)
                max_sim_to_cluster = max(max_sim_to_cluster, sim.get(pair, 0.0))
                
            if max_sim_to_cluster > best_candidate_sim:
                best_candidate_sim = max_sim_to_cluster
                best_candidate = candidate
                
        ordered.append(best_candidate)  # type: ignore[arg-type]
        remaining.remove(best_candidate)  # type: ignore[arg-type]
        
    return ordered

def build_tree_dict_from_kmers(
        kmer_dict: Dict[str, List[str]]
) -> Tuple[Dict[int, str], Dict[int, List[Any]]]:
        """
        Build a right-branching binary tree (as dictionaries) from a dictionary of
        k-mers, joining sequences in an order that favors higher k-mer similarity.

        Parameters
        ----------
        kmer_dict : Dict[str, List[str]]
                A dictionary mapping sequence IDs (e.g., 's1', 's2', ...) to lists of
                k-mer strings.

        Returns
        -------
        terminals_dic : Dict[int, str]
                Maps integer node IDs for terminal nodes to sequence labels.
                Example: {1: 's1', 2: 's2', 3: 's3'}

        internal_nodes : Dict[int, List[Any]]
                Maps integer node IDs for internal nodes to their children,
                which can be either sequence labels or integer IDs of previous
                internal nodes (right-branching pattern).
                Example (for 3 sequences):
                {4: ['s1', 's2'], 5: [4, 's3']}

        Notes
        -----
        - The order in which terminals appear is determined by k-mer similarity:
            first we pick the most similar pair, then iteratively add the sequence
            that is most similar to any already-added sequence.
        - The resulting structure is a simple right-branching guide tree, suitable
            for determining the order of progressive multiple sequence alignment.
        """
        # Determine an order of sequence IDs based on k-mer similarity
        ordered_terminals: List[str] = order_by_kmer_similarity(kmer_dict)
    
        # Build terminal node dictionary
        terminals_dic: Dict[int, str] = {}
        for i, name in enumerate(ordered_terminals, start=1):
                terminals_dic[i] = name
            
        n = len(ordered_terminals)
        internal_nodes: Dict[int, List[Any]] = {}
    
        if n == 0:
                return terminals_dic, internal_nodes
        if n == 1:
                # Only one sequence, no internal node
                return terminals_dic, internal_nodes
    
        # First internal node joins the first two terminals
        current_node = n + 1
        left, right = 1, 2
        internal_nodes[current_node] = [terminals_dic[left], terminals_dic[right]]
    
        # Subsequent internal nodes join the previous internal node with the next terminal
        for j in range(3, n + 1):
                current_node += 1
                # Left child is the previous internal node (by ID),
                # right child is the next terminal label.
                label1 = current_node - 1
                label2 = terminals_dic[j]
                internal_nodes[current_node] = [label1, label2]
            
        return terminals_dic, internal_nodes
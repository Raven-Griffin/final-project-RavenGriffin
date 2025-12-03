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
    Build a *more balanced* binary tree (as dictionaries) from a dictionary of
    k-mers, joining sequences in an order that favors higher k-mer similarity.

    Parameters
    ----------
    kmer_dict : Dict[str, List[str]]
        Maps sequence IDs (e.g., 's1', 's2', ...) to lists of k-mers.

    Returns
    -------
    terminals_dic : Dict[int, str]
        Maps integer terminal IDs (1..n) to sequence names.

    internal_nodes : Dict[int, List[Any]]
        Maps integer internal node IDs (n+1 .. n+(n-1)) to a list of two children.
        Each child is either:
          - a sequence name (string) for a leaf, or
          - an internal node ID (int) for a subtree.

    Notes
    -----
    - Uses `order_by_kmer_similarity` to get a similarity-aware ordering.
    - Builds a more balanced tree by pairing leaves/subtrees and then combining
      those pairs, instead of chaining everything onto the previous node.
    """
    # 1) Get a similarity-informed ordering of sequence IDs
    ordered_terminals: List[str] = order_by_kmer_similarity(kmer_dict)

    # 2) Build the terminal node dictionary: 1..n -> sequence name
    terminals_dic: Dict[int, str] = {
        i: name for i, name in enumerate(ordered_terminals, start=1)
    }

    n = len(ordered_terminals)
    internal_nodes: Dict[int, List[Any]] = {}

    if n == 0:
        return terminals_dic, internal_nodes
    if n == 1:
        # Only one sequence, nothing to join
        return terminals_dic, internal_nodes

    # 3) Create a queue of IDs representing current "nodes" in the tree.
    #    Start with terminal IDs 1..n.
    queue: List[int] = list(terminals_dic.keys())

    # Internal node IDs start after the terminal IDs: n+1, n+2, ...
    current_node_id = n

    # 4) Repeatedly take pairs from the queue and combine them into new internal nodes.
    #    Append each new internal node ID back to the queue.
    #
    #    This tends to produce a more balanced tree:
    #      - First combine (1,2), (3,4), ...
    #      - Then combine their parents, etc.
    #
    #    When there's an odd node out, it gets carried forward.
    while len(queue) > 1:
        # Take the first node
        left_id = queue.pop(0)

        # Take the second if present; if not, it just bubbles up
        if queue:
            right_id = queue.pop(0)
        else:
            # Nothing to pair with; push it back and break
            queue.insert(0, left_id)
            break

        current_node_id += 1

        # Convert IDs to labels:
        # - If it's a terminal ID, look up the sequence name.
        # - If it's an internal node ID, keep the int (subtree reference).
        def node_label(node_id: int) -> Any:
            return terminals_dic.get(node_id, node_id)

        left_label = node_label(left_id)
        right_label = node_label(right_id)

        # Record this internal node
        internal_nodes[current_node_id] = [left_label, right_label]

        # The new internal node participates in the next round
        queue.append(current_node_id)

    return terminals_dic, internal_nodes

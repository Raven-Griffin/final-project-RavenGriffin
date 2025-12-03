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
        
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap



class KmerAnalyzer:
    def __init__(self, kmer_dict):
        self.dict = kmer_dict
        self.all_kmers = None
        self.row_names = None
        self.col_names = None

    def count_kmers(self):
        overall = {}
        global_counts = {}

        for seq_id, kmers in self.dict.items():
            per_seq = {}
            for k in kmers:
                per_seq[k] = per_seq.get(k, 0) + 1
                global_counts[k] = global_counts.get(k, 0) + 1

            sorted_seq = sorted(per_seq.items(), key=lambda x: x[1], reverse=True)
            overall[seq_id] = sorted_seq

        self.all_kmers = {
            "per_sequence": overall,
            "global": sorted(global_counts.items(), key=lambda x: x[1], reverse=True)
        }

    def build_kmer_matrix(self, top_k=8):
        per_seq = self.all_kmers["per_sequence"]
        global_sorted = self.all_kmers["global"]

        top_kmers = [k for (k, _) in global_sorted[:top_k]]

        matrix = []
        row_names = []

        for seq_id, kmer_list in per_seq.items():
            d = dict(kmer_list)
            row = [d.get(k, 0) for k in top_kmers]
            matrix.append(row)
            row_names.append(seq_id)

        return row_names, top_kmers, np.array(matrix)

    def plot_kmer_heatmap(self, top_k=8, out_path="results/kmer_heatmap.png"):
        self.row_names, self.col_names, M = self.build_kmer_matrix(top_k)

        plt.figure(figsize=(12, 6))
        cmap = LinearSegmentedColormap.from_list(
            "purple_indigo", ["#C8A2C8", "#7a5195", "#4B0082"]
        )
        norm = plt.Normalize(vmin=0, vmax=M.max())
        plt.imshow(M, aspect='auto', cmap=cmap, norm=norm)

        cbar = plt.colorbar(label="k-mer frequency")
        ticks = list(range(0, int(M.max()) + 1))
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([str(t) for t in ticks])

        plt.xticks(range(len(self.col_names)), self.col_names, rotation=90)
        plt.yticks(range(len(self.row_names)), self.row_names)
        plt.title(f"Top {top_k} k-mer Frequency Heatmap")
        plt.tight_layout()
        
        plt.savefig(out_path, dpi=300)
        print(f"[Info]Saved heatmap PNG → {out_path}")
        
        plt.close()
       

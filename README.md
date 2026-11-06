# final-project-RavenGriffin

Name: Raven Griffin




## This project implements a complete bioinformatics workflow:

- Read nucleotide sequences from a FASTA file
- Detect longest ORFs across frames and orientations
- Build a progressive alignment guide tree using k-mer similarity
- Perform Needleman–Wunsch global alignments
- Produce a multiple sequence alignment (MSA)
- Compute mismatch and indel mutation rates per position
- Generate and save a k-mer heatmap as a PNG image

Everything is driven from **`main.py`**.


---

## 1. Installation

You’ll need **Python 3.10+**.

Install dependencies:


pip install biopython numpy pandas matplotlib pytest
## 2. Running Code 
There is a sample file called sequences.fasta

input:  python3 main.py -i sequences.fasta -o out.fasta -r s1

## 3. WorkFlow
### 3.1 Load Sequences
read_nucs() reads the input FASTA and stores sequences in a dictionary keyed by ID

### 3.2 ORF Detection
find_lorfs():
- Scans forward and reverse-complement strands
- Checks frames 0, 1, 2
- Translates using the chosen NCBI translation table
- Finds ORFs matching M ... *
- Keeps the longest ORF per sequence
- Stores ORFs as arrays of amino acid sets

### 3.3 Guide Tree via Kmers
From the ORFs:
- extract_kmers() builds kmer dictionarier per sequence
- build_tree_dict_from_kmers() constructs a simple progressive guide tree
- This determines the order in which pairwise alignments are performed

### 3.4 Needleman-Wunsch Alignment
Using functions from newu3.num
- init_mat() - innitalizes scoring matrix
- fill_matrix() - dynamic programming score fill
- trace_matrix() - traceback to get optimal alignment
  - pairwise alignments along the guide tree
  - gap propagation across previously aligned sequences (patch)
  - writing the final MSA to out.fasta

### 3.5 Mutation Rate Calculation
Using ratez.pd.analyze():
- Compares all seuences to a reference ID
- Computes per-position:
  - Mismatch rate
  - Indel rate
- Compute average mismatch and average indel rates acros positions
Results are saved to mutation_rates.csv

### 3.6 Kmer Heatmap
Using KmerAnalyzer in kmers/kmers/py
1. count_kmers()
   - Counts k-mer frequencies per sequence
   - Builds a global ranking of most frequent k-mers
2. build_kmer_matrix(top_k)
  - Builds a matrix of kmer counts:
      - Rows = sequences
      - Columns = top k kmers
3. plot_kmer_heatmap(top_k=8, out_path="results/kmer_heatmap.png")
  - Creates a heatmap using matplotlib
  - Automatically saves to results/kmer_heatmap.png
  - Does not pop up a window (if implemented as save-only)

## 5. Output Files
After running: python3 main.py -i sequences.fasta -o out.fasta -r s1
You will have:
- out.fasta - final multiple sequence alignment
- mutation_rates.csv
- results/kmer_heatmap.png
- results/mutation_rates_png

## 6. Running Tests
In bash run: pytest

## 7. Project Structure
```text
20251125_Alignment_take2/
├── main.py
├── newu3/
│   ├── __init__.py
│   ├── read.py
│   ├── write.py
│   ├── num.py
├── ratez/
│   ├── __init__.py
│   ├── bio.py
│   ├── pd.py
├── kmers/
│   ├── __init__.py
│   ├── kmers.py
├── tests/
│   ├── test_newu3_num.py
│   ├── test_orfs_and_tree.py
│   ├── test_pd_analyze_and_plot.py
│   ├── test_read_aln_and_ratescalc.py
├── results/
│   ├── kmer_heatmap.png (generated)
│   └── mutation_rates.csv (generated)
└── README.md
```
## 8. AI Usage
- Testing
- Debugging
- Error Handling




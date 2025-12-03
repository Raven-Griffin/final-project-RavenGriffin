# final-project-RavenGriffin

Name: Raven Griffin
Email: rgriff62@charlotte.edu
UNCC ID: 801282357



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
- 




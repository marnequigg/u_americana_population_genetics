from Bio import SeqIO
from Bio.Seq import Seq
import random

# Define IUPAC ambiguity codes and their possible resolutions
iupac_map = {
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G']
}

# Choose resolution strategy: "first", "random", or a fixed base like "A"
resolution_strategy = "random"  # options: "first", "random", or "A"

def resolve_base(base):
    base = base.upper()
    if base in iupac_map:
        if resolution_strategy == "first":
            return iupac_map[base][0]
        elif resolution_strategy == "random":
            return random.choice(iupac_map[base])
        else:
            return resolution_strategy  # e.g., always "A"
    return base

# Input and output FASTA paths
input_fasta = "/data/labs/Fant/Quigg/03e.hybpiper/04.hybpiper_pt2/glabra_sub/reference/CBG01_b2_consensus.fasta"
output_fasta = "/data/labs/Fant/Quigg/03e.hybpiper/04.hybpiper_pt2/glabra_sub/reference/ulmus_glabra_sub_ref.fasta"

# Convert and write new FASTA
with open(output_fasta, "w") as out_handle:
    for record in SeqIO.parse(input_fasta, "fasta"):
        resolved_seq = ''.join(resolve_base(b) for b in str(record.seq))
        record.seq = Seq(resolved_seq)
        SeqIO.write(record, out_handle, "fasta")

print(f"✅ Ambiguity codes resolved using strategy '{resolution_strategy}' → saved to {output_fasta}")

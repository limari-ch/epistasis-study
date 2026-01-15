#!/usr/bin/env python3
"""
trim_fasta.py

This script trims a multiple sequence alignment (MSA) to a specific region
based on nucleotide positions. It uses Biopython's AlignIO module.

Example usage:
    python trim_fasta.py input_msa.fasta output_msa_trimmed.fasta 106 281487
"""

from Bio import AlignIO
import sys

if len(sys.argv) != 5:
    print("Usage: python trim_fasta.py <input_fasta> <output_fasta> <start> <end>")
    sys.exit(1)

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])

# Read alignment
alignment = AlignIO.read(input_fasta, "fasta")

# Trim alignment
trimmed_alignment = alignment[:, start:end]

# Write output
AlignIO.write(trimmed_alignment, output_fasta, "fasta")

print(f"Trimmed alignment saved to: {output_fasta}")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MI_parallel.py

Compute pairwise Mutual Information (MI) across all sites in a multiple
sequence alignment (MSA) using parallel processing.

"""

import numpy as np
from Bio import AlignIO
from math import log
import csv
import multiprocessing as mp

# Global variables initialized in each child process
msa = None
freqs = None


def init_globals(msa_arg, freqs_arg):
    """
    Initialize global variables within each subprocess.
    This avoids sending large objects repeatedly to worker processes.
    """
    global msa, freqs
    msa = msa_arg
    freqs = freqs_arg


def compute_joint_frequency_pair(pos_pair):
    """
    Compute the Mutual Information (MI) for a given pair of positions (pos_i, pos_j)
    in the multiple sequence alignment (MSA).

    Parameters
    ----------
    pos_pair : tuple
        A pair of column indices (pos_i, pos_j) corresponding to two sites.

    Returns
    -------
    tuple
        (pos_i, pos_j, mi) where mi is the Mutual Information value.
    """
    pos_i, pos_j = pos_pair
    num_seqs = len(msa)
    joint_counts = {}

    # Count all observed symbol pairs across sequences
    for record in msa:
        a = record.seq[pos_i]
        b = record.seq[pos_j]
        key = (a, b)
        joint_counts[key] = joint_counts.get(key, 0) + 1

    # Compute joint frequency distribution
    joint_freq = {k: v / num_seqs for k, v in joint_counts.items()}

    # Calculate Mutual Information (base 2)
    mi = 0.0
    for (a, b), p_ab in joint_freq.items():
        p_a = freqs[pos_i].get(a, 1e-10)
        p_b = freqs[pos_j].get(b, 1e-10)
        if p_ab > 0:
            mi += p_ab * log(p_ab / (p_a * p_b), 2)

    return (pos_i, pos_j, mi)


def compute_frequencies(msa_local):
    """
    Compute marginal frequencies for each alignment position.

    Parameters
    ----------
    msa_local : Bio.Align.MultipleSeqAlignment
        The multiple sequence alignment loaded via Biopython.

    Returns
    -------
    tuple
        freqs_local : list of dicts
            Each dictionary stores the frequency of symbols for one alignment column.
        alignment_length : int
            Total number of positions in the alignment.
    """
    num_seqs = len(msa_local)
    alignment_length = msa_local.get_alignment_length()
    freqs_local = []

    for i in range(alignment_length):
        column = [record.seq[i] for record in msa_local]
        symbols, counts = np.unique(column, return_counts=True)
        freq = {s: count / num_seqs for s, count in zip(symbols, counts)}
        freqs_local.append(freq)

    return freqs_local, alignment_length


def main():
    """
    Main function to execute the MI computation workflow.

    Notes
    -----
    For demonstration purposes, this script processes only the first 100 positions.
    For full-scale analyses, generate all unique position pairs:
        L = alignment length
        Total pairs = L * (L - 1) / 2
    Be aware that a full MSA of ~10,000 sites produces ~50 million pairs.
    """
    # Path to the multiple sequence alignment (FASTA format)
    msa_file = "/Users/datasets/msa.fasta"
    msa_local = AlignIO.read(msa_file, "fasta")

    # Compute marginal frequencies
    freqs_local, L = compute_frequencies(msa_local)
    print("Alignment length:", L)

    # Demonstration: use only the first 100 positions
    L_demo = 100
    pos_pairs = [(i, j) for i in range(L_demo) for j in range(i + 1, L_demo)]
    print("Number of position pairs to process (demo):", len(pos_pairs))

    # Initialize multiprocessing pool
    pool = mp.Pool(processes=mp.cpu_count(),
                   initializer=init_globals,
                   initargs=(msa_local, freqs_local))

    # Parallel computation
    results = pool.map(compute_joint_frequency_pair, pos_pairs)
    pool.close()
    pool.join()

    # Save results to CSV
    output_file = "/Users/datasets/mi_results.csv"
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Position_i", "Position_j", "MutualInformation"])
        for pos_i, pos_j, mi_val in results:
            writer.writerow([pos_i, pos_j, mi_val])

    print("Results saved to:", output_file)


if __name__ == "__main__":
    main()

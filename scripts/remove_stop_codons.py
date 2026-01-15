#!/usr/bin/env python3
"""
Removes stop codons and ambiguous codons from a multiple sequence alignment (MSA).
This script ensures that all sequences retain a valid codon length (multiple of 3).

Codons containing gaps ('-'), dots ('.'), or ambiguous bases ('N') are ignored.
Sequences with invalid final lengths are reported and discarded.

Example usage:
    python remove_stop_codons.py input_alignment.fasta cleaned_alignment.fasta
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def clean_msa(input_fasta, output_fasta):
    """
    Removes stop codons and ambiguous codons from a multiple sequence alignment.

    Parameters
    ----------
    input_fasta : str
        Path to the input FASTA alignment file.
    output_fasta : str
        Path to the output FASTA file to save cleaned sequences.
    """

    # Set of stop codons and degenerate stop codon representations
    stop_codons = {
        "TAA", "TAG", "TGA",
        "TAR", "TGR", "TAN", "TGN",
        "TAATAG", "TGATAATAA", "TAATGA", "TGATAA", "TAATAA"
    }

    cleaned_records = []
    discarded = {}

    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = str(record.seq).upper()
        cleaned_codons = []
        removed_codons = 0

        # Process in codon triplets
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if "-" in codon or "." in codon or "N" in codon:
                continue  # skip codons with gaps or ambiguous bases
            if codon in stop_codons:
                removed_codons += 1
                continue
            cleaned_codons.append(codon)

        cleaned_seq = "".join(cleaned_codons)

        # Save only sequences with valid codon lengths
        if len(cleaned_seq) % 3 == 0 and cleaned_seq:
            cleaned_record = SeqRecord(
                Seq(cleaned_seq),
                id=record.id,
                description=""
            )
            cleaned_records.append(cleaned_record)
        else:
            discarded[record.id] = len(cleaned_seq)

    # Write the cleaned alignment
    SeqIO.write(cleaned_records, output_fasta, "fasta")

    # Summary report
    print(f"\n[✓] Cleaned sequences saved to: {output_fasta}")
    print(f"    Total cleaned: {len(cleaned_records)}")
    print(f"    Total discarded: {len(discarded)}")

    if discarded:
        print("\n[!] Discarded sequences (invalid final length):")
        for sid, slen in discarded.items():
            print(f" - {sid}: {slen} nt")


# === Example usage ===
if __name__ == "__main__":
    input_file = "/Users/path/msa.fasta"         # Replace with your input file path
    output_file = "/Users/path/cleaned_msa.fasta"  # Replace with your output file path
    clean_msa(input_file, output_file)

from Bio import SeqIO

def filter_sequences(input_file, output_file, threshold_percentage):
    with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            sequence = str(record.seq)
            n_count = sequence.count('N') + sequence.count('n')
            percentage_n = (n_count / len(sequence)) * 100

            if percentage_n < threshold_percentage:
                SeqIO.write(record, output_handle, "fasta")

# Example usage
input_file = "/Users/path/name_msa_trimed.fasta"
output_file = "/Users/path/output_name_msa_trimed.fasta"
threshold_percentage = 1.0  # Set your desired threshold percentage here

filter_sequences(input_file, output_file, threshold_percentage)

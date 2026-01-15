import pandas as pd
import numpy as np

def sample_sequences(file_path, output_path):
    # Read the TSV file
    df = pd.read_csv(file_path, delimiter="\t")
    
    # Check for required columns
    date_column = "date"  # Column containing the collection dates
    variant_column = "Lineage"  # Column containing the viral variant information

    if date_column not in df.columns or variant_column not in df.columns:
        raise KeyError(f"Required columns '{date_column}' or '{variant_column}' are missing. Available columns: {df.columns}")
    
    # Extract collection month from the date column
    df['Collection_Month'] = pd.to_datetime(df[date_column], errors='coerce').dt.to_period('M')
    
    # Drop rows with invalid dates
    df = df.dropna(subset=['Collection_Month'])
    
    # Group by month
    sampled_data = []
    for month, month_group in df.groupby('Collection_Month'):
        print(f"Processing month: {month}, Total sequences: {len(month_group)}")
        total_sequences = len(month_group)
        
        if total_sequences <= 10:
            # If 10 or fewer sequences, take all
            sampled_data.append(month_group)
        else:
            # Sample proportionally by variant
            month_group_counts = month_group[variant_column].value_counts()
            month_group_proportions = month_group_counts / total_sequences

            samples_per_variant = (month_group_proportions * 10).round().astype(int)
            
            # Fail-safe to prevent infinite loop
            max_attempts = 1000
            attempts = 0

            while samples_per_variant.sum() != 10:
                if samples_per_variant.sum() < 10:
                    difference = (month_group_proportions * 10 - samples_per_variant).sort_values(ascending=False)
                    samples_per_variant[difference.idxmax()] += 1
                elif samples_per_variant.sum() > 10:
                    difference = (samples_per_variant - month_group_proportions * 10).sort_values(ascending=True)
                    candidate = difference.idxmin()
                    if samples_per_variant[candidate] > 1:
                        samples_per_variant[candidate] -= 1

                # Check if stuck
                attempts += 1
                if attempts >= max_attempts:
                    print(f"Warning: Sampling adjustment for month {month} reached max attempts. Proceeding with current allocation.")
                    break

            # Sample the required number of sequences per variant
            for variant, n_samples in samples_per_variant.items():
                sampled_data.append(month_group[month_group[variant_column] == variant].sample(n=n_samples, random_state=42))
    
    # Combine all sampled data
    final_sampled_data = pd.concat(sampled_data, ignore_index=True)
    
    # Save the sampled data to a new TSV file
    final_sampled_data.to_csv(output_path, sep="\t", index=False)
    print("Sampling complete. Output saved to:", output_path)

# File paths
file_path = "/Users/path/name_metadata.tsv"  # Replace with your input TSV file path
output_path = "/Users/path/random_name_metadata.tsv"  # Replace with your output TSV file path

# Run the function
sample_sequences(file_path, output_path)

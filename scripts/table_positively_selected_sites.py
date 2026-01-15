import pandas as pd

# Load each CSV file
df_california = pd.read_csv("/Users/path/California_datamonkey_table_hyphy.csv")
df_florida = pd.read_csv("/Users/path/Florida_datamonkey_table_hyphy.csv")
df_puertorico = pd.read_csv("/Users/path/PuertoRico_datamonkey_table_hyphy.csv")
df_newyork = pd.read_csv("/Users/path/NewYork_datamonkey_table_hyphy.csv")
df_hawaii = pd.read_csv("/Users/path/Hawaii_datamonkey_table_hyphy.csv")

# Add a column to identify the source
df_california["Location"] = "California"
df_florida["Location"] = "Florida"
df_puertorico["Location"] = "Puerto Rico"
df_newyork["Location"] = "New York"
df_hawaii["Location"] = "Hawaii"

# Combine all into one DataFrame
combined_df = pd.concat([df_california, df_florida, df_puertorico, df_newyork, df_hawaii], ignore_index=True)

# Filter significant sites (p-value < 0.05)
significant_sites = combined_df[combined_df["p-value"] < 0.05]

# Select key columns to display
columns_to_show = ["Site", "Location", "LRT", "p-value", "&beta;<sup>+</sup>", "p<sup>+</sup>"]
significant_table = significant_sites[columns_to_show]

# Save or print the table
print(significant_table)
# To save:
significant_table.to_csv("/Users/path/significant_sites_table.csv", index=False)

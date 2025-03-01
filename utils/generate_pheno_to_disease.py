import pandas as pd

# Read the input file into a DataFrame
df = pd.read_csv('phenotype_to_genes.txt', sep='\t')

# Create a new DataFrame with the required format
output_df = df[['hpo_id', 'disease_id']]

# Save the output DataFrame to a new file
output_df.to_csv('pheno_to_diseases.txt', sep='\t', index=False, header=False)

print("File transformation complete.")

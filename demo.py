import json
import pandas as pd

# Open the JSONL file and read it line by line
with open('simulated_patients.json', 'r') as file:
    data = []
    for line in file:
        # Parse each line as a JSON object and append it to the data list
        data.append(json.loads(line.strip()))

# Flatten the JSON data
flattened_data = []
for record in data:
    disease_id = record.get('disease_id', '')
    true_genes = ','.join(record.get('true_genes', []))
    age = record.get('age', '')
    positive_phenotypes = ','.join(record.get('positive_phenotypes', {}).keys())
    negative_phenotypes = ','.join(record.get('negative_phenotypes', {}).keys())
    n_distractor_genes = record.get('n_distractor_genes', '')
    distractor_genes = ','.join(record.get('distractor_genes', {}).keys())
    dropout_phenotypes = ','.join(record.get('dropout_phenotypes', {}).keys())
    corruption_phenotypes = ','.join(record.get('corruption_phenotypes', {}).keys())
    id_val = record.get('id', '')

    # Append each record as a list
    flattened_data.append([disease_id, true_genes, age, positive_phenotypes, negative_phenotypes, n_distractor_genes, distractor_genes, dropout_phenotypes, corruption_phenotypes, id_val])

# Create DataFrame
df = pd.DataFrame(flattened_data, columns=['disease_id', 'true_genes', 'age', 'positive_phenotypes', 'negative_phenotypes', 'n_distractor_genes', 'distractor_genes', 'dropout_phenotypes', 'corruption_phenotypes', 'id'])

# Save to CSV
df.to_csv('output.csv', index=False)
print('CSV file saved as output.csv')

import csv

# Input and output file names
input_file = "genes_to_disease.txt"
output_file = "mapped_genes.txt"

# Open the input and output files
with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    reader = csv.reader(infile, delimiter="\t")  # Tab-separated file
    next(reader)  # Skip header
    
    outfile.write("gene_symbol\tDisease_id\n")  # Write new header
    
    for row in reader:
        if len(row) >= 4:  # Ensure row has necessary data
            gene_symbol = row[1]  # Second column
            disease_id = row[3]   # Fourth column
            outfile.write(f"{gene_symbol}\t{disease_id}\n")

print(f"Output written to {output_file}")
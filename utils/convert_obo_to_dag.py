import re

def parse_obo_to_dag(obo_file, output_file):
    # Open the OBO file for reading
    with open(obo_file, 'r') as f:
        obo_data = f.read()

    # Split the data into individual terms
    terms = obo_data.split('[Term]')
    
    # Initialize a dictionary to hold 'id' -> 'is_a' relationships
    dag = []

    for term in terms:
        # Extract term ID
        term_id = re.search(r'id:\s*(HP:\d+)', term)
        if term_id:
            term_id = term_id.group(1)
        
            # Find all 'is_a' relationships in the term
            is_a_matches = re.findall(r'is_a:\s*(HP:\d+)', term)

            # For each 'is_a' relationship, create a (child, parent) pair in DAG
            for is_a in is_a_matches:
                dag.append((term_id, is_a))

    # Write the DAG to the output file
    with open(output_file, 'w') as out_file:
        for relation in dag:
            out_file.write(f"{relation[0]}\t{relation[1]}\n")
    
    print(f"DAG successfully written to {output_file}")

# Usage
obo_file = "./hp.obo"  # Update with your actual OBO file path
output_file = "output_dag.txt"  # The file to save the DAG relations
parse_obo_to_dag(obo_file, output_file)

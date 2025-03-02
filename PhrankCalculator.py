import pronto
import math
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
import argparse
import sys
import os
import json

class PhrankCalculator:
    def __init__(self, obo_file, hpoa_file, disease_gene_file=None):
        self.ontology = self.load_ontology(obo_file)
        self.child_to_parent, self.parent_to_children = self.build_dag_maps()
        self.disease_to_hpo, self.disease_to_name = self.read_disease_annotations(hpoa_file)
        self.disease_to_gene = self.load_disease_gene(disease_gene_file) if disease_gene_file else {}
        self.gene_pheno_map = self.compute_gene_pheno_map()
        self.hpo_prob = self.hpo_term_probability()
        self.IC, self.marginal_IC = self.compute_information_content()
        # Save the DAG to a file
        self.save_dag_to_file("output_dag.txt")
        self.save_disease_annotations("output.txt")

    def load_ontology(self, path_to_obo):
        """Load the HPO ontology from an OBO file."""
        print(f"Loading ontology from {path_to_obo}...")
        return pronto.Ontology(path_to_obo)

    def build_dag_maps(self):
        """Build directed acyclic graph maps for parent-child relationships."""
        print("Building DAG maps...")
        child_to_parent = defaultdict(list)
        parent_to_children = defaultdict(list)
        root_term = "HP:0000118"
        
        # Iterate through all terms using ontology.terms() for explicit access
        for term in self.ontology.terms():
            if term.id == root_term:
                continue
            # Get direct parents (not all superclasses)
            for parent in term.superclasses(distance=1):
                if parent.id != term.id:  # Avoid self-loops
                    child_to_parent[term.id].append(parent.id)
                    parent_to_children[parent.id].append(term.id)
        return child_to_parent, parent_to_children

    def read_disease_annotations(self, hpo_disease_annotations):
        """Read disease to HPO annotations from a file."""
        print(f"Reading disease annotations from {hpo_disease_annotations}...")
        disease_to_hpo = defaultdict(set)
        disease_to_name = {}
        
        with open(hpo_disease_annotations, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                disease_id, disease_name, _, hpo_id = parts[:4]
                disease_to_hpo[disease_id].add(hpo_id)
                disease_to_name[disease_id] = disease_name
        return disease_to_hpo, disease_to_name

    def load_disease_gene(self, disease_gene_file):
        """Load disease to gene mappings from a file."""
        print(f"Loading disease-gene mappings from {disease_gene_file}...")
        disease_to_gene = defaultdict(set)
        with open(disease_gene_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                gene, disease = parts[:2]
                disease_to_gene[disease].add(gene)
        return disease_to_gene

    def compute_gene_pheno_map(self):
        """Compute gene to phenotype mappings based on disease annotations."""
        print("Computing gene-phenotype mappings...")
        gene_pheno_map = defaultdict(set)
        for disease, genes in self.disease_to_gene.items():
            for gene in genes:
                gene_pheno_map[gene].update(self.disease_to_hpo.get(disease, set()))
        return gene_pheno_map

    def hpo_term_probability(self):
        """Calculate probability of each HPO term based on disease annotations."""
        print("Calculating HPO term probabilities...")
        hpo_counts = defaultdict(int)
        total_diseases = len(self.disease_to_hpo)
        
        for hpo_set in self.disease_to_hpo.values():
            for hpo in hpo_set:
                hpo_counts[hpo] += 1
                
        return {hpo: count/total_diseases for hpo, count in hpo_counts.items()}

    def compute_information_content(self):
        """Compute information content for each HPO term."""
        print("Computing information content...")
        IC = {}
        marginal_IC = {}
        
        # Calculate standard IC
        for hpo, prob in self.hpo_prob.items():
            IC[hpo] = -math.log2(prob) if prob > 0 else 0
        
        # Calculate marginal IC
        for hpo in IC:
            parents = self.child_to_parent.get(hpo, [])
            if not parents:
                marginal_IC[hpo] = IC[hpo]
            else:
                parent_ICs = [IC.get(parent, 0) for parent in parents]
                max_parent_IC = max(parent_ICs) if parent_ICs else 0
                marginal_IC[hpo] = IC[hpo] - max_parent_IC
                
        return IC, marginal_IC

    def get_closure(self, phenotypes):
        """Get ancestor closure for a set of phenotypes."""
        closure = set(phenotypes)
        stack = list(phenotypes)
        
        while stack:
            current = stack.pop()
            parents = self.child_to_parent.get(current, [])
            for parent in parents:
                if parent not in closure:
                    closure.add(parent)
                    stack.append(parent)
        return closure

    def compute_phenotype_similarity(self, query_phenotypes, target_phenotypes):
        """Compute similarity between two sets of phenotypes."""
        query_closure = self.get_closure(query_phenotypes)
        target_closure = self.get_closure(target_phenotypes)
        overlap = query_closure & target_closure
        return sum(self.marginal_IC.get(p, 0) for p in overlap)

    def rank_diseases(self, patient_phenotypes, top_n=10):
        """Rank diseases based on similarity to patient phenotypes."""
        print(f"Ranking diseases for {len(patient_phenotypes)} patient phenotypes...")
        scores = []
        for disease, disease_phenos in self.disease_to_hpo.items():
            score = self.compute_phenotype_similarity(patient_phenotypes, disease_phenos)
            scores.append((score, disease, self.disease_to_name[disease]))
        return sorted(scores, reverse=True)[:top_n]

    def rank_genes(self, patient_phenotypes, top_n=10):
        """Rank genes based on similarity to patient phenotypes."""
        print(f"Ranking genes for {len(patient_phenotypes)} patient phenotypes...")
        gene_scores = defaultdict(float)
        for gene, gene_phenos in self.gene_pheno_map.items():
            score = self.compute_phenotype_similarity(patient_phenotypes, gene_phenos)
            gene_scores[gene] = max(gene_scores[gene], score)
        sorted_genes = sorted(gene_scores.items(), key=lambda x: x[1], reverse=True)
        return sorted_genes[:top_n]

    def suggest_terms(self, top_diseases, patient_phenotypes, n=3):
        """Suggest additional phenotype terms to consider."""
        print("Suggesting additional terms...")
        suggestions = []
        for _, disease_id, _ in top_diseases[:3]:
            disease_terms = self.disease_to_hpo[disease_id]
            for term in disease_terms:
                if term not in patient_phenotypes:
                    suggestions.append((term, self.marginal_IC.get(term, 0)))
        return sorted(suggestions, key=lambda x: x[1], reverse=True)[:n]
    
    def save_dag_to_file(self, output_file):
        """Save the DAG structure to a file for compatibility with original Phrank."""
        print(f"Saving DAG to {output_file}...")
        with open(output_file, 'w') as f:
            # Format: term<tab>parentTerm1,parentTerm2,...
            for child, parents in self.child_to_parent.items():
                f.write(f"{child}\t{','.join(parents)}\n")
            
            # Add terms with no parents
            for term in self.ontology.terms():
                if term.id not in self.child_to_parent and term.id != "HP:0000001":  # Skip root
                    f.write(f"{term.id}\t\n")
    
    def save_disease_annotations(self, output_file):
        """Save disease annotations to a file for compatibility with original Phrank."""
        print(f"Saving disease annotations to {output_file}...")
        with open(output_file, 'w') as f:
            # Format: disease<tab>phenotype1,phenotype2,...
            for disease, phenotypes in self.disease_to_hpo.items():
                f.write(f"{disease}\t{','.join(phenotypes)}\n")

    def visualize_dag(self, focus_term=None, depth=2, max_nodes=50):
        """Visualize the HPO DAG structure."""
        print("Visualizing DAG...")
        G = nx.DiGraph()
        
        nodes_to_include = set()
        term_labels = {}
        
        if focus_term:
            # Start with focus term and explore to limited depth
            current_level = {focus_term}
            nodes_to_include.update(current_level)
            term_labels[focus_term] = self.get_term_label(focus_term)
            
            for _ in range(depth):
                next_level = set()
                for term in current_level:
                    # Add children
                    children = self.parent_to_children.get(term, [])
                    for child in children[:5]:  # Limit to avoid overcrowding
                        next_level.add(child)
                        nodes_to_include.add(child)
                        term_labels[child] = self.get_term_label(child)
                    
                    # Add parents
                    parents = self.child_to_parent.get(term, [])
                    for parent in parents:
                        next_level.add(parent)
                        nodes_to_include.add(parent)
                        term_labels[parent] = self.get_term_label(parent)
                
                current_level = next_level
                if len(nodes_to_include) > max_nodes:
                    break
        else:
            # Include the root and some important categories
            root = "HP:0000118"  # Phenotypic abnormality
            nodes_to_include.add(root)
            term_labels[root] = self.get_term_label(root)
            
            # Add main categories
            for child in self.parent_to_children.get(root, []):
                nodes_to_include.add(child)
                term_labels[child] = self.get_term_label(child)
                
                # Add some grandchildren
                for grandchild in self.parent_to_children.get(child, [])[:3]:
                    nodes_to_include.add(grandchild)
                    term_labels[grandchild] = self.get_term_label(grandchild)
        
        # Create the graph
        for node in nodes_to_include:
            G.add_node(node)
            for parent in self.child_to_parent.get(node, []):
                if parent in nodes_to_include:
                    G.add_edge(node, parent)  # Child -> Parent direction
        
        # Visualize
        plt.figure(figsize=(12, 10))
        pos = nx.spring_layout(G, k=0.8, iterations=50, seed=42)
        
        nx.draw_networkx_nodes(G, pos, node_color='lightblue', node_size=300, alpha=0.8)
        nx.draw_networkx_edges(G, pos, edge_color='gray', arrows=True, 
                              arrowstyle='->', arrowsize=10, width=1.0)
        
        # Adjust label positions slightly
        label_pos = {node: (coords[0], coords[1] - 0.02) for node, coords in pos.items()}
        nx.draw_networkx_labels(G, label_pos, labels=term_labels, font_size=6)
        
        plt.title(f"HPO Ontology DAG{' (focused on '+focus_term+')' if focus_term else ''}")
        plt.axis('off')
        plt.tight_layout()
        
        # Save the visualization
        output_file = f"hpo_dag{'_'+focus_term.replace(':','_') if focus_term else ''}.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved visualization to {output_file}")
        plt.close()
        
        return output_file
    
    def get_term_label(self, term_id):
        """Get a shortened label for a term."""
        if term_id in self.ontology:
            name = self.ontology[term_id].name
            return f"{term_id}\n{name[:15]}..." if len(name) > 15 else f"{term_id}\n{name}"
        return term_id
    
    def visualize_disease_similarity_network(self, patient_phenotypes, top_n=5):
        """Visualize the disease similarity network for a patient."""
        print("Visualizing disease similarity network...")
        G = nx.Graph()
        
        # Get top diseases
        top_diseases = self.rank_diseases(patient_phenotypes, top_n)
        
        # Add patient node
        G.add_node("Patient", type="patient")
        
        # Add disease nodes and edges
        for score, disease_id, name in top_diseases:
            G.add_node(disease_id, type="disease", name=name, score=score)
            G.add_edge("Patient", disease_id, weight=score)
            
            # Add a few key phenotypes for this disease
            for hpo in list(self.disease_to_hpo[disease_id])[:3]:
                G.add_node(hpo, type="phenotype", 
                          name=self.ontology[hpo].name if hpo in self.ontology else hpo)
                G.add_edge(disease_id, hpo)
        
        # Add patient phenotypes
        for pheno in patient_phenotypes:
            G.add_node(pheno, type="patient_phenotype",
                      name=self.ontology[pheno].name if pheno in self.ontology else pheno)
            G.add_edge("Patient", pheno)
        
        # Visualize
        plt.figure(figsize=(12, 10))
        
        # Node colors and sizes
        node_colors = []
        node_sizes = []
        for node in G.nodes():
            node_type = G.nodes[node].get('type', '')
            if node_type == 'patient':
                node_colors.append('red')
                node_sizes.append(1000)
            elif node_type == 'disease':
                # Color by score (blue gradient)
                score = G.nodes[node].get('score', 0)
                intensity = min(255, int(score * 30))
                node_colors.append(f'#{0:02x}{0:02x}{intensity:02x}')
                node_sizes.append(800)
            elif node_type == 'patient_phenotype':
                node_colors.append('green')
                node_sizes.append(500)
            else:  # Regular phenotype
                node_colors.append('yellow')
                node_sizes.append(400)
        
        # Edge weights
        edge_weights = [G.edges[e].get('weight', 1) * 2 for e in G.edges()]
        
        # Layout
        pos = nx.spring_layout(G, k=0.3, iterations=50, seed=42)
        
        # Draw
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.8)
        nx.draw_networkx_edges(G, pos, width=edge_weights, alpha=0.5, edge_color='gray')
        
        # Labels
        labels = {}
        for node in G.nodes():
            if G.nodes[node].get('type') == 'disease':
                labels[node] = G.nodes[node].get('name', node)
            elif G.nodes[node].get('type') in ['phenotype', 'patient_phenotype']:
                term_name = G.nodes[node].get('name', '')
                if len(term_name) > 15:
                    term_name = term_name[:15] + "..."
                labels[node] = f"{node}\n{term_name}"
            else:
                labels[node] = node
        
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=8)
        
        # Legend
        legend_elements = [
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=15, label='Patient'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=15, label='Disease'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=15, label='Patient Phenotype'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='yellow', markersize=15, label='Disease Phenotype')
        ]
        plt.legend(handles=legend_elements, loc='upper right')
        
        plt.title("Patient-Disease-Phenotype Similarity Network")
        plt.axis('off')
        plt.tight_layout()
        
        # Save
        output_file = "disease_similarity_network.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved visualization to {output_file}")
        plt.close()
        
        return output_file


class Phrank:
    """
    Implementation compatible with the original Phrank from the Bejerano lab.
    This is based on the code at https://bitbucket.org/bejerano/phrank/src/master/
    """
    def __init__(self, dagfile, diseaseannotationsfile=None, weightsfile=None):
        self.root = "HP:0000118"  # Phenotypic abnormality
        self.dag = self.read_dag(dagfile)
        self.disease_to_pheno = self.read_disease_pheno(diseaseannotationsfile) if diseaseannotationsfile else {}
        self.weights = self.read_weights(weightsfile) if weightsfile else {}
        
        # Compute closures if disease annotations are provided
        if self.disease_to_pheno:
            self.disease_to_closure = self.compute_closures(self.disease_to_pheno)
        
        # Compute information content
        self.hpo_probabilities = self.compute_probabilities()
        self.term_ic = self.compute_information_content()
    
    def read_dag(self, dagfile):
        """Read the DAG structure from a file."""
        dag = defaultdict(list)
        with open(dagfile, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                term = parts[0]
                parents = parts[1].split(',') if len(parts) > 1 and parts[1] else []
                dag[term] = parents
        return dag
    
    def read_disease_pheno(self, file_path):
        """Read disease to phenotype mappings."""
        disease_to_pheno = defaultdict(set)
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                disease = parts[0]
                phenotypes = parts[1].split(',') if len(parts) > 1 and parts[1] else []
                disease_to_pheno[disease].update(phenotypes)
        return disease_to_pheno
    
    def read_weights(self, file_path):
        """Read term weights if provided."""
        weights = {}
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    term, weight = parts
                    weights[term] = float(weight)
        return weights
    
    def get_closure(self, terms):
        """Compute ancestor closure for a set of terms."""
        closure = set(terms)
        for term in terms:
            stack = [term]
            while stack:
                current = stack.pop()
                parents = self.dag.get(current, [])
                for parent in parents:
                    if parent and parent not in closure:
                        closure.add(parent)
                        stack.append(parent)
        return closure
    
    def compute_closures(self, entity_to_terms):
        """Compute closures for all entities."""
        entity_to_closure = {}
        for entity, terms in entity_to_terms.items():
            entity_to_closure[entity] = self.get_closure(terms)
        return entity_to_closure
    
    def compute_probabilities(self):
        """Compute term probabilities for information content."""
        if not self.disease_to_pheno:
            return {}  # No disease annotations available
            
        term_counts = defaultdict(int)
        total_diseases = len(self.disease_to_pheno)
        
        # Count term occurrences in closures
        for disease, terms in self.disease_to_closure.items():
            for term in terms:
                term_counts[term] += 1
        
        # Convert to probabilities
        return {term: count/total_diseases for term, count in term_counts.items()}
    
    def compute_information_content(self):
        """Compute information content for each term."""
        return {term: -math.log2(prob) if prob > 0 else 0 
                for term, prob in self.hpo_probabilities.items()}
    
    def compute_similarity(self, query_terms, target_terms):
        """Compute semantic similarity using information content."""
        query_closure = self.get_closure(query_terms)
        target_closure = self.get_closure(target_terms)
        
        # Find overlapping terms
        overlap = query_closure & target_closure
        
        # Sum the information content
        similarity = sum(self.term_ic.get(term, 0) for term in overlap)
        
        # Apply weights if available
        if self.weights:
            weight_factor = sum(self.weights.get(term, 1.0) for term in query_terms)
            similarity *= weight_factor
            
        return similarity
    
    def rank_diseases(self, query_genes, query_phenotypes, top_n=10):
        """Rank diseases by similarity to query phenotypes."""
        if not self.disease_to_pheno:
            return []  # No disease data available
            
        results = []
        for disease, phenotypes in self.disease_to_pheno.items():
            similarity = self.compute_similarity(query_phenotypes, phenotypes)
            results.append((disease, similarity))
        
        # Sort by similarity, descending
        sorted_results = sorted(results, key=lambda x: x[1], reverse=True)
        return sorted_results[:top_n]
    
    def export_results(self, results, output_file):
        """Export ranking results to a file."""
        with open(output_file, 'w') as f:
            for disease, score in results:
                f.write(f"{disease}\t{score}\n")


def main():
    parser = argparse.ArgumentParser(description='Phrank: Phenotype ranking tool')
    parser.add_argument('--obo', required=True, help='Path to HPO .obo file')
    parser.add_argument('--hpoa', required=True, help='Path to phenotype.hpoa file')
    parser.add_argument('--genes', help='Path to gene-disease mapping file')
    parser.add_argument('--patient', required=True, help='File with patient phenotypes, one per line')
    parser.add_argument('--visualize', action='store_true', help='Generate visualizations')
    parser.add_argument('--output', default='results', help='Output directory')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    # Load patient phenotypes
    patient_phenotypes = []
    with open(args.patient, 'r') as f:
        for line in f:
            hpo_id = line.strip()
            if hpo_id and not hpo_id.startswith('#'):
                patient_phenotypes.append(hpo_id)
    
    print(f"Loaded {len(patient_phenotypes)} patient phenotypes")
    
    # Initialize and run Phrank
    calculator = PhrankCalculator(args.obo, args.hpoa, args.genes)
    
    # Generate results
    disease_rankings = calculator.rank_diseases(patient_phenotypes)
    gene_rankings = calculator.rank_genes(patient_phenotypes) if args.genes else []
    term_suggestions = calculator.suggest_terms(disease_rankings, patient_phenotypes)
    
    # Save results
    with open(os.path.join(args.output, 'disease_rankings.txt'), 'w') as f:
        for score, disease_id, name in disease_rankings:
            f.write(f"{disease_id}\t{name}\t{score}\n")
    
    if gene_rankings:
        with open(os.path.join(args.output, 'gene_rankings.txt'), 'w') as f:
            for gene, score in gene_rankings:
                f.write(f"{gene}\t{score}\n")
    
    with open(os.path.join(args.output, 'term_suggestions.txt'), 'w') as f:
        for term, ic in term_suggestions:
            name = calculator.ontology[term].name if term in calculator.ontology else ""
            f.write(f"{term}\t{name}\t{ic}\n")
    
    # Generate visualizations if requested
    if args.visualize:
        print("Generating visualizations...")
        calculator.visualize_dag()
        calculator.visualize_dag(focus_term=patient_phenotypes[0] if patient_phenotypes else None)
        calculator.visualize_disease_similarity_network(patient_phenotypes)
    
    print(f"Results saved to {args.output} directory")


if __name__ == "__main__":
    # Example usage
    if len(sys.argv) > 1:
        main()
    else:
        # Demo with sample data
        print("Running demo with sample data...")
        phrank = PhrankCalculator(
            obo_file='data/hp.obo',
            hpoa_file='data/phenotype.hpoa',
            disease_gene_file='data/genes_to_disease.txt'
        )
        
        # Example patient phenotypes
        patient_phenotypes = {
            'HP:0000824', 'HP:0000831', 'HP:0000938', 'HP:0001260',
            'HP:0001266', 'HP:0001268', 'HP:0001332', 'HP:0001510',
            'HP:0001596', 'HP:0003077', 'HP:0005135', 'HP:0040171',
            'HP:0008619', 'HP:0008669', 'HP:0008734', 'HP:0100840',
            'HP:0001627', 'HP:0000062', 'HP:0010646'
        }
        
        # Generate visualizations
        phrank.visualize_dag()
        phrank.visualize_dag(focus_term='HP:0000707')  # Abnormality of the nervous system
        phrank.visualize_disease_similarity_network(patient_phenotypes)
        
        # Get disease rankings
        disease_ranking = phrank.rank_diseases(patient_phenotypes)
        print("\nTop Disease Matches:")
        for score, disease_id, name in disease_ranking:
            print(f"{name} (Score: {score:.2f})")
        
        # # Get gene rankings
        # gene_ranking = phrank.rank_genes(patient_phenotypes)
        # print("\nTop Gene Matches:")
        # for gene, score in gene_ranking:
        #     print(f"{gene} (Score: {score:.2f})")
#!/usr/bin/env python3
import pronto
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
import random
import math

def visualize_hpo_dag(obo_file, max_nodes=100, focus_term=None, depth_limit=3):
    """
    Visualize the HPO ontology as a DAG.
    
    Parameters:
    - obo_file: Path to the HPO .obo file
    - max_nodes: Maximum number of nodes to display (to avoid overcrowding)
    - focus_term: Optional HPO term ID to focus the visualization around (e.g., 'HP:0000118')
    - depth_limit: When focus_term is provided, how many levels to display
    """
    # Load ontology
    print("Loading HPO ontology...")
    ontology = pronto.Ontology(obo_file)
    
    # Build DAG maps
    child_to_parent = defaultdict(list)
    parent_to_children = defaultdict(list)
    
    print("Building DAG structure...")
    # Use HP:0000118 (Phenotypic abnormality) as the root
    root_term = "HP:0000118" if focus_term is None else focus_term
    
    # Create graph
    G = nx.DiGraph()
    
    # Add nodes and edges
    terms_to_process = set()
    term_labels = {}
    
    if focus_term:
        # Start with the focus term and explore limited depth
        current_terms = {focus_term}
        seen_terms = set(current_terms)
        
        for depth in range(depth_limit + 1):
            next_terms = set()
            
            for term_id in current_terms:
                if term_id in ontology:
                    term = ontology[term_id]
                    term_labels[term_id] = f"{term_id}\n{term.name[:20]+'...' if len(term.name) > 20 else term.name}"
                    terms_to_process.add(term_id)
                    
                    # Add children
                    for child in term.subclasses():
                        if child.id != term_id:  # Avoid self-loops
                            child_to_parent[child.id].append(term_id)
                            parent_to_children[term_id].append(child.id)
                            if child.id not in seen_terms:
                                next_terms.add(child.id)
                                seen_terms.add(child.id)
                    
                    # Add parents
                    for parent in term.superclasses():
                        if parent.id != term_id:  # Avoid self-loops
                            child_to_parent[term_id].append(parent.id)
                            parent_to_children[parent.id].append(term_id)
                            if parent.id not in seen_terms:
                                next_terms.add(parent.id)
                                seen_terms.add(parent.id)
            
            current_terms = next_terms
            if len(terms_to_process) > max_nodes:
                break
    else:
        # Add all terms up to max_nodes
        for term in list(ontology.terms())[:max_nodes]:
            term_id = term.id
            term_labels[term_id] = f"{term_id}\n{term.name[:20]+'...' if len(term.name) > 20 else term.name}"
            terms_to_process.add(term_id)
            
            # Add parent relationships
            for parent in term.superclasses():
                if parent.id != term_id:  # Avoid self-loops
                    child_to_parent[term_id].append(parent.id)
                    parent_to_children[parent.id].append(term_id)
    
    # Limit to max_nodes if we have too many
    if len(terms_to_process) > max_nodes:
        terms_to_process = set(list(terms_to_process)[:max_nodes])
    
    # Add nodes and edges to graph
    for term_id in terms_to_process:
        G.add_node(term_id)
        for parent_id in child_to_parent[term_id]:
            if parent_id in terms_to_process:
                G.add_edge(term_id, parent_id)  # Direction is from child to parent
    
    # Visualize
    plt.figure(figsize=(16, 12))
    pos = nx.spring_layout(G, seed=42, k=0.8)  # k controls spacing
    
    # Draw nodes
    nx.draw_networkx_nodes(G, pos, 
                          node_color='lightblue', 
                          node_size=500, 
                          alpha=0.8)
    
    # Draw edges
    nx.draw_networkx_edges(G, pos, 
                          edge_color='gray', 
                          arrows=True, 
                          arrowstyle='->', 
                          arrowsize=15, 
                          width=1.0,
                          alpha=0.6)
    
    # Draw labels with smaller font and better positioning
    label_pos = {node: (coords[0], coords[1] - 0.02) for node, coords in pos.items()}
    nx.draw_networkx_labels(G, label_pos, 
                           labels=term_labels, 
                           font_size=8, 
                           font_family='sans-serif',
                           verticalalignment='top')
    
    plt.title(f"HPO Ontology DAG Visualization{' (focused on '+focus_term+')' if focus_term else ''}")
    plt.axis('off')
    plt.tight_layout()
    
    # Save and show
    output_file = f"hpo_dag{'_'+focus_term.replace(':','_') if focus_term else ''}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Visualization saved to {output_file}")
    plt.show()

def visualize_gene_disease_network(phrank_calculator, patient_phenotypes=None, top_n=5):
    """
    Visualize the gene-disease-phenotype network based on the PhrankCalculator
    and optionally highlight connections to a patient's phenotypes.
    
    Parameters:
    - phrank_calculator: Initialized PhrankCalculator object
    - patient_phenotypes: Optional set of patient HPO terms
    - top_n: Number of top diseases/genes to include 
    """
    G = nx.Graph()
    
    # If patient phenotypes are provided, get rankings
    if patient_phenotypes:
        disease_ranking = phrank_calculator.rank_diseases(patient_phenotypes, top_n)
        gene_ranking = phrank_calculator.rank_genes(patient_phenotypes, top_n)
        
        # Add patient node
        G.add_node("Patient", type="patient")
        
        # Add disease nodes from ranking
        for _, disease_id, name in disease_ranking:
            G.add_node(disease_id, type="disease", name=name)
            G.add_edge("Patient", disease_id, type="patient_disease")
            
            # Add HPO terms from this disease
            for hpo in phrank_calculator.disease_to_hpo[disease_id]:
                if hpo in patient_phenotypes:
                    G.add_node(hpo, type="shared_phenotype", 
                              name=phrank_calculator.ontology[hpo].name if hpo in phrank_calculator.ontology else hpo)
                    G.add_edge("Patient", hpo, type="patient_phenotype")
                    G.add_edge(disease_id, hpo, type="disease_phenotype")
            
            # Add associated genes
            for gene in phrank_calculator.disease_to_gene.get(disease_id, []):
                if len([g for g, _ in gene_ranking if g == gene]) > 0:  # If gene is in top genes
                    G.add_node(gene, type="gene")
                    G.add_edge(disease_id, gene, type="disease_gene")
    else:
        # Sample visualization of the overall network structure
        # Take a sample of diseases
        sample_diseases = list(phrank_calculator.disease_to_hpo.keys())[:top_n]
        
        for disease_id in sample_diseases:
            G.add_node(disease_id, type="disease", name=phrank_calculator.disease_to_name.get(disease_id, disease_id))
            
            # Add sample HPO terms
            sample_hpos = list(phrank_calculator.disease_to_hpo[disease_id])[:3]
            for hpo in sample_hpos:
                G.add_node(hpo, type="phenotype", 
                          name=phrank_calculator.ontology[hpo].name if hpo in phrank_calculator.ontology else hpo)
                G.add_edge(disease_id, hpo, type="disease_phenotype")
            
            # Add associated genes
            for gene in list(phrank_calculator.disease_to_gene.get(disease_id, []))[:2]:
                G.add_node(gene, type="gene")
                G.add_edge(disease_id, gene, type="disease_gene")
    
    # Visualize
    plt.figure(figsize=(16, 12))
    
    # Define node colors and sizes by type
    node_colors = []
    node_sizes = []
    for node in G.nodes():
        node_type = G.nodes[node].get('type', '')
        if node_type == 'patient':
            node_colors.append('red')
            node_sizes.append(1000)
        elif node_type == 'disease':
            node_colors.append('skyblue')
            node_sizes.append(800)
        elif node_type == 'gene':
            node_colors.append('lightgreen')
            node_sizes.append(600)
        elif node_type == 'shared_phenotype':
            node_colors.append('orange')
            node_sizes.append(500)
        else:  # regular phenotype
            node_colors.append('yellow')
            node_sizes.append(400)
    
    # Define edge colors by type
    edge_colors = []
    for u, v in G.edges():
        edge_type = G.edges[u, v].get('type', '')
        if edge_type == 'patient_disease':
            edge_colors.append('red')
        elif edge_type == 'patient_phenotype':
            edge_colors.append('orange')
        elif edge_type == 'disease_gene':
            edge_colors.append('green')
        else:  # disease_phenotype
            edge_colors.append('blue')
    
    # Position nodes
    pos = nx.spring_layout(G, seed=42, k=0.5)
    
    # Draw the network
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.8)
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=1.5, alpha=0.6)
    
    # Create labels with names when available
    labels = {}
    for node in G.nodes():
        if G.nodes[node].get('type') == 'disease':
            labels[node] = G.nodes[node].get('name', node)
        elif G.nodes[node].get('type') in ['phenotype', 'shared_phenotype']:
            hpo_name = G.nodes[node].get('name', node)
            labels[node] = f"{node}\n{hpo_name[:15]}..." if len(hpo_name) > 15 else f"{node}\n{hpo_name}"
        else:
            labels[node] = node
    
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, font_family='sans-serif')
    
    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=15, label='Patient'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='skyblue', markersize=15, label='Disease'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='lightgreen', markersize=15, label='Gene'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', markersize=15, label='Shared Phenotype'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='yellow', markersize=15, label='Phenotype')
    ]
    plt.legend(handles=legend_elements, loc='upper right')
    
    plt.title("Gene-Disease-Phenotype Network" + (" (Patient Analysis)" if patient_phenotypes else ""))
    plt.axis('off')
    plt.tight_layout()
    
    # Save and show
    output_file = "gene_disease_network.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Visualization saved to {output_file}")
    plt.show()

# Example usage
if __name__ == "__main__":
    # 1. Visualize the HPO ontology DAG
    visualize_hpo_dag('hp.obo', max_nodes=50)
    
    # 2. Visualize a specific section of the HPO ontology
    visualize_hpo_dag('hp.obo', focus_term='HP:0000707', depth_limit=2)  # Abnormality of the nervous system
    
    # 3. Visualize gene-disease network using a PhrankCalculator
    from paste import PhrankCalculator  # Import from your original file
    
    phrank = PhrankCalculator(
        obo_file='hp.obo',
        hpoa_file='phenotype.hpoa',
        disease_gene_file='genes_to_disease.txt'
    )
    
    # Example with patient phenotypes
    patient_phenotypes = {
        'HP:0000059', 'HP:0000164', 'HP:0000248', 'HP:0000252',
        'HP:0000276', 'HP:0000308', 'HP:0000411', 'HP:0000448'
    }
    
    visualize_gene_disease_network(phrank, patient_phenotypes)
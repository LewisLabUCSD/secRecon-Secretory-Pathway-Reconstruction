import requests
import json
import pickle
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from tqdm import tqdm


############################################
## Color mapping dictionary for plotting ###
############################################

system_colors = {
    'Protein conformation': (1.0, 0.6, 0.0, 1.0),
    'Post-translational modifications': (0.1, 0.9, 0.1, 1.0), 
    'Proteostasis': (0.3, 0.3, 0.9, 1.0),
    'Translocation': (0.4, 0.7, 0.9, 1.0),
    'Vesicle trafficking': (0.7, 0.7, 0.3, 1.0),
    'Secreted': (0.9, 0.1, 0.1, 1.0)
}

localization_colors = {
    'Actin Cytoskeleton': (0.7, 0.4, 0.2, 1.0),
    'Centrosome': (0.0, 0.5, 0.5, 1.0),
    'ERGIC': (1.0, 0.6, 0.0, 1.0),
    'Endosome': (0.0, 0.7, 0.3, 1.0), 
    'Recycling Endosome': (0.3, 0.0, 0.7, 1.0),
    'Late Endosome': (0.8, 0.0, 0.5, 1.0),
    'Early Endosome': (0.7, 0.0, 0.4, 1.0), 
    'Plasma Membrane': (0.5, 0.5, 0.0, 1.0), 
    'Golgi': (0.1, 0.9, 0.1, 1.0), 
    'cis-Golgi': (0.2, 0.8, 0.6, 1.0),
    'trans-Golgi': (0.4, 0.8, 0.4, 1.0), 
    'Nucleus': (0.9, 0.1, 0.1, 1.0),
    'Nucleolus': (0.5, 0.1, 0.1, 1.0),
    'Phagosome': (0.6, 0.3, 0.0, 1.0),
    'Proteasome': (0.3, 0.3, 0.9, 1.0),
    'Mitochondria': (0.9, 0.5, 0.5, 1.0),
    'Lysosome': (0.5, 0.9, 0.8, 1.0),  
    'Cytosol': (0.3, 0.8, 0.95, 1.0),
    'Cytoplasm': (0.4, 0.7, 0.9, 1.0),
    'Endoplasmic Reticulum': (0.7, 0.7, 0.3, 1.0),
    'Secreted': (0.8, 0.4, 0.0, 1.0),
    'Stress Granule': (0.9, 0.3, 0.05, 1.0),
    'Translation': (0.8, 0.2, 0.8, 1.0),
    'Unknown': (0.0, 0.0, 0.0, 1.0)
}

def get_radius(value, categories):
    """
    Determine the radius of the node based on the normalized_log value.
    
    Parameters:
    value (float): The normalized_log value of the node.
    
    Returns:
    float: The radius corresponding to the given value based on predefined categories.
    """
    # Iterate over the categories to find the appropriate radius for the given value
    for (low, high), radius in categories.items():
        if low <= value < high:
            return radius
    
    # Return a default small value if the value does not fit into any category
    return 0.001

def fetch_string_interactions(main_genes, extra_genes=[]):
    """
    Fetch protein-protein interactions from STRING database for a list of genes.
    
    Parameters:
    - main_genes (list): List of gene names to fetch interactions for.
    - extra_genes (list): Additional genes to consider in interactions.
    
    Returns:
    - G (networkx.Graph): A graph where nodes are proteins and edges represent interactions.
    """
    # Initialize an empty graph
    G = nx.Graph()

    # Combine main_genes and extra_genes for convenience
    all_genes = set(main_genes + extra_genes)
    
    # Initialize set to track seen interactions
    seen_interactions = set() 
    
    # Fetch interactions for each gene in the list
    for gene in tqdm(main_genes, desc="Fetching Interactions"):
        url = f"https://string-db.org/api/json/network?identifiers={gene}"
        response = requests.get(url)
        
        if response.status_code == 200:
            interactions = json.loads(response.text)
            
            for interaction in interactions:
                protein1 = interaction['preferredName_A']
                protein2 = interaction['preferredName_B']
                interaction_tuple = (protein1, protein2)

                # Check if interaction is already seen, considering both directions
                if interaction_tuple not in seen_interactions and (protein2, protein1) not in seen_interactions:
                    if protein1 in all_genes and protein2 in all_genes:
                        score = interaction['score']
                        # Add nodes and edges to the graph if interaction is new
                        G.add_node(protein1)
                        G.add_node(protein2)
                        G.add_edge(protein1, protein2, weight=score)
                        # Add to seen interactions in both orders
                        seen_interactions.add(interaction_tuple)
                        seen_interactions.add((protein2, protein1))
        else:
            print(f'No interaction found for {gene}')
                
    return G


def visualize_network(G, gene_dict, pos, node_size=0.010, filename=None, color_by='systems', legends=True, categories=None):
    """
    Visualize a protein-protein interaction network using matplotlib.

    Parameters:
    - G (networkx.Graph): The graph to visualize.
    - pos (dict): Positions of nodes.
    - gene_dict (dict): Dictionary containing gene information.
    - node_size (float, optional): Size of the nodes.
    - filename (str, optional): If provided, save the plot to this filename.
    - color_by (str, optional): Attribute to color nodes by ('systems' or 'localization').
    - legends (bool, optional): Whether to display legends.
    """
    # Create a figure and axes
    fig, ax = plt.subplots(figsize=(40, 40))

    # Scale edge widths
    edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
    min_width = 0.1
    max_width = 4.0
    epsilon = 1e-10
    if edge_weights:
        edge_weights = [min_width + (w - min(edge_weights)) * (max_width - min_width) / (max(edge_weights) - min(edge_weights) + epsilon) for w in edge_weights]

    # Draw edges with scaled widths
    nx.draw_networkx_edges(G, pos, width=edge_weights, edge_color='lightgray', ax=ax)

    # Draw nodes as pie charts
    if color_by == 'systems':
        for node, (x, y) in pos.items():
            systems = gene_dict[node]['systems']
            colors = [system_colors[sys] for sys in systems]

            if node_size == 'exp':
                # Draw pie chart at node position with edgecolor and linewidth
                raw_radius = gene_dict[node]['Expression']
                radius = get_radius(raw_radius, categories)
                ax.pie([1]*len(systems), colors=colors, radius=radius, center=(x, y), wedgeprops=dict(edgecolor='black', linewidth=0.5))
            else:
                ax.pie([1]*len(systems), colors=colors, radius=node_size, center=(x, y), wedgeprops=dict(edgecolor='black', linewidth=0.5))

    elif color_by == 'localization':
        for node, (x, y) in pos.items():
            localizations = gene_dict[node]['subcellular_localization']
            colors = [localization_colors[loc] for loc in localizations]

            if node_size == 'exp':
                # Draw pie chart at node position with edgecolor and linewidth
                raw_radius = gene_dict[node]['Expression']
                radius = get_radius(raw_radius, categories)
                ax.pie([1]*len(localizations), colors=colors, radius=radius, center=(x, y), wedgeprops=dict(edgecolor='black', linewidth=0.5))
            else:
                ax.pie([1]*len(localizations), colors=colors, radius=0.012, center=(x, y), wedgeprops=dict(edgecolor='black', linewidth=0.5))

    # Get the current axis limits
    x_values, y_values = zip(*pos.values())
    min_x, max_x = min(x_values), max(x_values)
    min_y, max_y = min(y_values), max(y_values)

    # Set new axis limits
    ax.set_xlim(min_x - 0.1, max_x + 0.1)
    ax.set_ylim(min_y - 0.1, max_y + 0.1)

    # Legend
    if legends:
        if color_by == 'systems':
            legend_patches = [mpatches.Patch(color=color, label=category) for category, color in system_colors.items()]
            color_legend = plt.legend(handles=legend_patches, prop={'size': 30}, loc='lower left', bbox_to_anchor=(0.6, 0.8), title="System" , title_fontsize=35)
        elif color_by == 'localization':
            legend_patches = [mpatches.Patch(color=color, label=category) for category, color in localization_colors.items()]
            color_legend = plt.legend(handles=legend_patches, prop={'size': 25}, loc='lower left', bbox_to_anchor=(0.7, 0.4), title="Subcellular Localization" , title_fontsize=35)

        # Create custom legend handles for node sizes
        size_legend_patches = []
        size_labels = ['0 to 2', '2 to 4', '> 4']
        scaling_factor = 1800

        for (low, high), radius, label in zip(categories.keys(), categories.values(), size_labels):
            handle = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=radius * scaling_factor, label=label)
            size_legend_patches.append(handle)

        # Add the second legend for node sizes
        size_legend = plt.legend(handles=size_legend_patches, prop={'size': 40}, loc='lower left', bbox_to_anchor=(0, 0.3), title="Log2 Fold Change", title_fontsize=35)

        # Add both legends to the plot
        ax.add_artist(color_legend)
        ax.add_artist(size_legend)

        plt.subplots_adjust(right=0.75)

    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')

    plt.show()
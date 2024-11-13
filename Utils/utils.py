import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from collections import Counter

def integrate_dicts(gene_dict, process_dict):
    """
    Integrate gene_dict with additional system, subsystem, and process information from process_dict
    based on matching processes, organizing them under separate subkeys for each gene.

    Args:
    gene_dict (dict): Dictionary of genes with their processes and subcellular localization.
    process_dict (dict): Dictionary categorizing biological processes into systems and subsystems.

    Returns:
    dict: Updated gene_dict with separate subkeys for system, subsystem, and process information.
    """
    
    # Function to find and update system, subsystem, and process information for a given process name.
    def update_system_subsystem_process(gene_details, process_name):
        for system, process_list in process_dict.items():
            for process in process_list:
                # Check if the process is in the category of the systems
                if system == process_name:
                    gene_details['system'].add(system)
                # Check if the process is in the category of the subsystems
                elif process['Subsystem'] == process_name:
                    gene_details['system'].add(system)
                    gene_details['subsystem'].add(process['Subsystem'])
                # Check if the process matches the given process_name and update accordingly.
                elif 'Process' in process and process['Process'] == process_name:
                    gene_details['system'].add(system)
                    gene_details['subsystem'].add(process['Subsystem'])
                    gene_details['process'].add(process_name)
                elif 'Subprocess' in process and process['Subprocess'] == process_name:
                    gene_details['system'].add(system)
                    gene_details['subsystem'].add(process['Subsystem'])
                    gene_details['process'].add(process['Process'])
    
    # Iterate through each gene in gene_dict and update with system, subsystem, and process info.
    for gene, details in gene_dict.items():
        # Initialize sets to ensure uniqueness.
        details['system'] = set()
        details['subsystem'] = set()
        details['process'] = set()
        
        for process in details['processes']:
            update_system_subsystem_process(details, process)  # Update details with the found information.
        
        # Convert sets back to lists for JSON compatibility.
        details['system'] = list(details['system'])
        details['subsystem'] = list(details['subsystem'])
        details['process'] = list(details['process'])
    
    return gene_dict

    

def get_gene_color(gene, gene_dict, process_dict, category_colors):

    """
    Function to determine the color of a gene based on its associated processes.

    It takes into consideration all processes associated with the gene, and assigns 
    a color based on the parent process that has the majority of processes associated 
    with the gene.

    Args:
    gene (str): The gene for which the color needs to be determined.
    gene_dict (dict): A dictionary with genes as keys and associated processes as values.
    process_dict (dict): A dictionary with categories (parent processes) as keys and associated processes as values.
    category_colors (dict): A dictionary with categories as keys and associated colors as values.

    Returns:
    str: The color (as a string) associated with the category that has the majority of processes associated with the gene.
    """

    # Initialize a dict to count category occurrences
    category_counts = {}
    
    # Loop over categories and processes
    for category, processes in process_dict.items():
        # Check if any process is in the gene's processes
        for process in flatten_processes(category, processes):
            if process in gene_dict.get(gene, [])['processes']:
                # If process is found, increment the count for the category
                category_counts[category] = category_counts.get(category, 0) + 1

    # Get the category with maximum count
    max_category = max(category_counts, key=category_counts.get, default=None)

    # If no category is found (i.e., the gene has no processes), return 'black'
    if max_category is None:
        return 'black'
    else:
        # Return the color for the category with maximum count
        return category_colors.get(max_category, 'black')


def flatten_processes(category, process_list):
    flat_list = []
    for process in process_list:
        if isinstance(process, dict):
            flat_list.extend(flatten_processes(category, list(process.values())))
        elif isinstance(process, list):
            flat_list.extend(flatten_processes(category, process))
        else:
            flat_list.append(process)
            
    flat_list.append(category)
    flat_list = list(set(flat_list))
    
    return flat_list

def adjust_color_alpha(hex_color, alpha):
    # Convert the hex color to RGB
    rgb_color = mcolors.hex2color(hex_color)
    # Add the alpha component
    rgba_color = (*rgb_color, alpha)
    return rgba_color


def categorize_location(text):
    """
    Categorizes the cellular location of a protein based on the text description.
    
    Parameters:
        text (str): The text description containing information about the cellular location.
        
    Returns:
        list: A list of unique categories representing the cellular locations.
    """

    # Check for missing values
    if text is None or (isinstance(text, float) and np.isnan(text)):
        return []
    
    categories = []  # List to store the categories of cellular locations
    text = text.lower()  # Convert the text to lowercase for uniformity
    
    # Remove text in curly braces and split by punctuation marks
    clean_text = ''.join(c for c in text if c not in '{}[]()')
    locations = clean_text.split(';')[0]
    locations = locations.split('.')
    locations = [location.strip() for loc in locations for location in loc.split(',')]
    
    # Remove text from notes
    text = text.split('note=')[0]
    
    # Check for the presence of the keywords in each parsed location
    for location in locations:
        if 'golgi' in location:
            if ('cis-golgi') in location:
                categories.append('cis-Golgi')
            elif ('trans-golgi') in location:
                categories.append('trans-Golgi')
            else:
                categories.append('Golgi')
        elif 'cytoplasm' in location:
            categories.append('Cytoplasm')
        elif 'nucleus' in location:
            categories.append('Nucleus')
        elif ('mitochondrion' or 'mitochondria') in location:
            categories.append('Mitochondria')
        elif ('endosome') in location:
            if ('early endosome') in location:
                categories.append('Early Endosome')
            elif ('late endosome') in location:
                categories.append('Late Endosome')
            elif ('recycling endosome') in location:
                categories.append('Recycling Endosome')
            else:
                categories.append('Endosome')
        elif ('lysosome') in location:
            categories.append('Lysosome')
        elif ('endoplasmic reticulum') in location:
            categories.append('Endoplasmic Reticulum')
        elif ('membrane') in location:
            categories.append('Plasma Membrane')
        elif ('phagosome') in location:
            categories.append('Phagosome')
            
    return list(set(categories))  # Using set to remove duplicate categories, if any




def identify_enriched_processes(clusters_dict, gene_dict, process_key='processes'):
    """
    Identifies and normalizes the enrichment of processes across gene clusters.
    
    This function calculates the count of each process within each cluster,
    normalizes these counts by the total occurrences of each process in the
    original dataset, and returns a DataFrame with these normalized counts
    alongside their raw counts for comparison.

    Parameters:
    - clusters_dict (dict): A dictionary where each key is a cluster identifier and
      each value is a dict of genes, with their associated data including processes.
    - gene_dict (dict): A dictionary representing the original dataset, where each key
      is a gene identifier and each value is a dict of gene data, including processes.
    - process_key (str): The key used in gene data dicts to access the list of processes.
      Defaults to 'processes'.

    Returns:
    - pd.DataFrame: A DataFrame containing the cluster identifier, process names,
      raw counts of processes in each cluster, total counts of processes in the
      original dataset, and the normalized count (raw count / total count).
    """
    cluster_data = []
    original_process_counts = Counter()

    # Count occurrences of each process in the original dataset to establish a baseline for normalization
    for gene, data in gene_dict.items():
        if data.get(process_key):
            original_process_counts.update(data[process_key])

    # Analyze each cluster to count occurrences of processes and prepare data for normalization
    for cluster_id, genes in clusters_dict.items():
        all_processes = [data[process_key] for gene, data in genes.items() if data.get(process_key)]
        process_counts = Counter([process for sublist in all_processes for process in sublist])
        
        # Collect data for each process in the current cluster, including total counts from the original dataset
        for process, count in process_counts.items():
            cluster_data.append({
                'Cluster': cluster_id,
                'Process': process,
                'Count': count,
                'Total_Count': original_process_counts.get(process, 0)  # Default to 0 if process not found
            })

    df = pd.DataFrame(cluster_data)

    # Normalize the counts (raw count in cluster / total count in original dataset)
    # Replace 0 with pd.NA to avoid division by zero and maintain data integrity
    df['Normalized_Count'] = df['Count'] / df['Total_Count'].replace(0, pd.NA)

    return df


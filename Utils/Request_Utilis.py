from Bio import Entrez
import time
import requests
from bs4 import BeautifulSoup
import json

def get_entrez_id(gene_symbol, database="gene", organism="Homo sapiens"):
    search_term = f"{gene_symbol}[Gene Name] AND {organism}[Organism]"
    
    try:
        handle = Entrez.esearch(db=database, term=search_term)
        record = Entrez.read(handle)
        handle.close()

        # Check if any IDs were found
        if record["IdList"]:
            return record["IdList"][0]  # Returns the first ID found
        else:
            return "No ID found for given gene symbol"
    except Exception as e:
        return f"An error occurred: {str(e)}"


def read_gene_ID(gene_EntrezID):
    gene_fetch_tries = 0
    gene_fetch_found = True
    record = []
    while gene_fetch_tries < 3 and gene_fetch_found:
        try:
            search_results = Entrez.esearch(db="gene", term = gene_EntrezID)            
            record = Entrez.read(search_results)
            gene_fetch_found = False
        except Exception as e:
            time.sleep(3)
            gene_fetch_tries += 1
    if gene_fetch_tries > 2:        
        print(f"EntrezID of {gene_EntrezID} could not be read; An error occurred")

    if int(record["Count"]) > 0:

        # Get the Gene ID of the first result
        gene_id = record["IdList"][0]

        # Fetch the gene information, maximum 3 tries
        fetch_tries = 0
        fetch_found = True
        while fetch_tries < 3 and fetch_found:
            gene_fetch_tries = 0
            gene_fetch_found = True
            gene_info = []
            while gene_fetch_tries < 3 and gene_fetch_found:
                try:
                    gene_info = Entrez.efetch(db="gene", id=gene_id, retmode="xml")    
                    gene_fetch_found = False
                except Exception as e:
                    time.sleep(3)
                    gene_fetch_tries += 1
                    print(f"Gene ID {gene_id}; An error occurred: {e}")
            if gene_fetch_tries > 2:
                print(f"{gene_EntrezID}. Gene ID {gene_id}; An error occurred")

            if gene_info == []:
                KeyError("sda")
            
            try:
                gene_record = Entrez.read(gene_info)
                fetch_found = False
            except Exception as e:
                time.sleep(3)
                fetch_tries += 1                
    return gene_record

def description_fetch_and_try(human_ortholog_EntrezID):
    human_fetch_tries = 0
    human_fetch_found = True
    gene_description = []

    while human_fetch_tries < 3 and human_fetch_found:
        try:
            search_results_human = Entrez.esearch(db="gene", term=human_ortholog_EntrezID)     
            record_human = Entrez.read(search_results_human)
            gene_id_human = record_human["IdList"][0]
            gene_info_human = Entrez.efetch(db="gene", id=gene_id_human, retmode="xml")
            gene_record_human = Entrez.read(gene_info_human)
            gene_description = gene_record_human[0]['Entrezgene_summary']
            human_fetch_found = False
        except Exception as e:
            time.sleep(3)
            human_fetch_tries += 1
            # print(f"Human_Description {human_ortholog_EntrezID}; An error occurred: {e}")     
    if human_fetch_tries > 2:
        print(f"Human_Description {human_ortholog_EntrezID}; An error occurred")


    return gene_description

def get_gene_ids(gene_entrezID, target_tax_id):
    url = f'https://www.ncbi.nlm.nih.gov/gene/{gene_entrezID}/ortholog/'
    response = requests.get(url)
    if response.status_code == 200:
        soup = BeautifulSoup(response.content, 'html.parser')
        script_tag = soup.find('script', string=lambda text: 'var appData =' in text)
        script_content = script_tag.text if script_tag else ""
        gene_data_start = script_content.find('appData.genes =') + len('appData.genes =')
        gene_data_end = script_content.find(';', gene_data_start)
        gene_data_json = script_content[gene_data_start:gene_data_end].strip()

        # Ensure gene_data_json is not empty
        if gene_data_json:
            try:
                genes_data = json.loads(gene_data_json)
                for gene_info in genes_data:
                    if gene_info.get('tax_id') == int(target_tax_id):
                        return gene_info.get('gene_id')
            except json.JSONDecodeError as e:
                print(f"Error decoding JSON: {e}")
                return None
        else:
            print("No gene data found in the page.")
            return None
    else:
        print(f'No accession for gene {gene_entrezID}')
        return None

    return None
    
def Gene_Info_from_EntrezID(EntrezID):
    protein_sequence_ID_list = []
    transcript_sequence_ID_list = []
    gene_ensemble = ''
    gene_record = read_gene_ID(int(EntrezID))
    gene_products = []
    
    org = gene_record[0]['Entrezgene_source']['BioSource']['BioSource_org']['Org-ref']['Org-ref_taxname']
    if 'Gene-ref_syn' in gene_record[0]['Entrezgene_gene']['Gene-ref']:
        gene_synonyms = gene_record[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_syn']
    else:
        gene_synonyms = []

    # Get the gene name and gene NCBI ID
    if 'Gene-ref_locus' in gene_record[0]['Entrezgene_gene']['Gene-ref']:
        gene_symbol = gene_record[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']                        
    else:
        gene_symbol = gene_record[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_locus-tag']

    # Get the Gene description
    if 'Gene-ref_desc' in gene_record[0]['Entrezgene_gene']['Gene-ref']:
        gene_name = gene_record[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_desc']
    else:
        try:
            gene_name = gene_record[0]['Entrezgene_locus'][0]['Gene-commentary_products'][0]['Gene-commentary_label']
        except KeyError:
            print(f'No description anotated for gene {gene_symbol}')


    # Get Ensemble ID
    try:
        for gene_db_refs in gene_record[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_db']:
            if gene_db_refs['Dbtag_db'] == 'Ensembl':
                gene_ensemble = (gene_db_refs['Dbtag_tag']['Object-id']['Object-id_str'])
    except KeyError:
        print(f"No ENSEMBL ID for gene {gene_symbol}")

    for assembly_specific_info in gene_record[0]['Entrezgene_locus']:  
        try:
            if 'Gene-commentary_heading' in assembly_specific_info:
                for assembly_specific_transcript in assembly_specific_info['Gene-commentary_products']:
                    if 'Gene-commentary_products' in assembly_specific_transcript:
                        transcript_sequence_ID = assembly_specific_transcript['Gene-commentary_accession']

                        for Entrez_Comments in gene_record[0]['Entrezgene_comments']:
                            if 'Gene-commentary_comment' in Entrez_Comments:
                                for comments in Entrez_Comments['Gene-commentary_comment']:
                                    if 'Gene-commentary_products' in comments:
                                        for product_per_assembly in comments['Gene-commentary_products']:                                
                                            if product_per_assembly['Gene-commentary_heading'] == 'mRNA Sequence':
                                                mRNA = product_per_assembly['Gene-commentary_accession']
                                                if mRNA == transcript_sequence_ID:  
                                                    for protein_per_transcript_per_assemlby in product_per_assembly['Gene-commentary_products']:
                                                        protein_sequence_ID = protein_per_transcript_per_assemlby['Gene-commentary_accession']

                                                        uniprotID_list = []
                                                        if 'Gene-commentary_comment' in protein_per_transcript_per_assemlby:
                                                            for protein_comments in protein_per_transcript_per_assemlby['Gene-commentary_comment']:
                                                                if protein_comments['Gene-commentary_heading'] == 'UniProtKB':
                                                                    for uniprotID_protein_per_assembly in protein_comments['Gene-commentary_comment'][0]['Gene-commentary_source']:
                                                                        protein_uniprot_id = uniprotID_protein_per_assembly['Other-source_src']['Dbtag']['Dbtag_tag']['Object-id']['Object-id_str']
                                                                        uniprotID_list.append(protein_uniprot_id)
                                                        gene_products.append(((transcript_sequence_ID, protein_sequence_ID, uniprotID_list)))                                       
                                            else:   
                                                if 'Gene-commentary_products' in product_per_assembly:
                                                    for transcript_per_assembly in product_per_assembly['Gene-commentary_products']:      
                                                        mRNA = transcript_per_assembly['Gene-commentary_accession']          
                                                        if mRNA == transcript_sequence_ID:                                                  
                                                            for protein_per_transcript_per_assemlby in transcript_per_assembly['Gene-commentary_products']:                                                                
                                                                protein_sequence_ID = protein_per_transcript_per_assemlby['Gene-commentary_accession']
                                                                uniprotID_list = []
                                                                if 'Gene-commentary_comment' in protein_per_transcript_per_assemlby:
                                                                    for protein_comments in protein_per_transcript_per_assemlby['Gene-commentary_comment']:                                                         
                                                                        if protein_comments['Gene-commentary_heading'] == 'UniProtKB':
                                                                            for uniprotID_protein_per_assembly in protein_comments['Gene-commentary_comment'][0]['Gene-commentary_source']:
                                                                                protein_uniprot_id = uniprotID_protein_per_assembly['Other-source_src']['Dbtag']['Dbtag_tag']['Object-id']['Object-id_str']
                                                                                uniprotID_list.append(protein_uniprot_id)
                                                                else:
                                                                    uniprotID_list = []
                                                                gene_products.append(((transcript_sequence_ID, protein_sequence_ID, uniprotID_list)))
        except KeyError:
            print(f'Gene {gene_symbol} has no products')

    return org, gene_symbol, gene_name, gene_synonyms, gene_ensemble, gene_products



def get_subcellular_localization(uniprot_id):
    """
    Corrected function to fetch the subcellular location of a protein using the UniProt REST API,
    based on the updated understanding of the JSON structure.
    
    Parameters:
    - uniprot_id (str): The UniProt ID of the protein.
    
    Returns:
    - str: The subcellular location of the protein or an error message.
    """
    api_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    
    try:
        response = requests.get(api_url)
        if response.status_code == 200:
            data = response.json()
            for comment in data.get("comments", []):
                if comment.get("commentType") == "SUBCELLULAR LOCATION":
                    locations = comment.get("subcellularLocations", [])
                    if locations:
                        location_descriptions = [loc["location"]["value"] for loc in locations if "location" in loc]
                        return location_descriptions
            return None
        elif response.status_code in [400, 404]:
            return response.json().get("messages", ["Error occurred"])[0]
        else:
            return "Failed to fetch data due to an error on the server."
    except requests.RequestException as e:
        return f"Error fetching data for UniProt ID {uniprot_id}: {e}"
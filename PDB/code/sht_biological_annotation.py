#!/usr/bin/env python3
"""
sht_biological_annotation.py

This script enriches the SHT outlier report with biological annotations by
querying the RCSB PDB GraphQL API. It aims to answer three key questions:
1. Are the supercoiling patterns sequence-dependent?
2. Are the patterns conserved across species?
3. Do the patterns correlate with biological function?

Workflow:
1.  Load the outlier report CSV.
2.  For each PDB ID, construct and execute a GraphQL query to fetch:
    - DNA sequence
    - Source organism
    - Protein function (from keywords, descriptions, GO terms)
    - Chemical modifications
3.  Define functions to parse this data and categorize it.
4.  Merge the annotations back into the outlier dataframe.
5.  Save the final, enriched report to a new CSV file.
"""

import argparse
from pathlib import Path
import pandas as pd
import requests
import time
from tqdm import tqdm

def get_pdb_graphql_query(pdb_id):
    """
    Returns the GraphQL query string for a given PDB ID.
    """
    return f"""
    {{
      entry(entry_id: "{pdb_id}") {{
        struct_keywords {{
          pdbx_keywords
        }}
        polymer_entities {{
          entity_poly {{
            pdbx_seq_one_letter_code_can
            type
          }}
          rcsb_polymer_entity_container_identifiers {{
            source_organism_names
          }}
          rcsb_entity_source_organism {{
            ncbi_scientific_name
          }}
          rcsb_polymer_entity {{
            pdbx_description
          }}
          rcsb_go_terms {{
            name
          }}
        }}
      }}
    }}
    """

def analyze_sequence(sequence):
    """
    Calculates GC content and CpG count for a DNA sequence.
    """
    if not sequence or not isinstance(sequence, str):
        return 0.0, 0
    
    seq_upper = sequence.upper()
    g_count = seq_upper.count('G')
    c_count = seq_upper.count('C')
    total_len = len(seq_upper)
    
    gc_content = (g_count + c_count) / total_len if total_len > 0 else 0.0
    cpg_count = seq_upper.count('CG')
    
    return gc_content, cpg_count

def categorize_function(entity_info, keywords):
    """
    Categorizes the biological function based on keywords and GO terms.
    """
    # Get the nested description and GO terms correctly
    description = entity_info.get('rcsb_polymer_entity', {}).get('pdbx_description', '')
    go_terms = entity_info.get('rcsb_go_terms', [])

    # Combine all text sources for searching
    search_text = " ".join([
        str(description),
        str(keywords),
        " ".join([go['name'] for go in go_terms if go])
    ]).lower()

    if any(k in search_text for k in ['polymerase', 'replication']):
        return "Replication"
    if any(k in search_text for k in ['repair', 'glycosylase', 'ligase', 'nuclease']):
        return "Repair"
    if any(k in search_text for k in ['transcription', 'factor', 'repressor']):
        return "Transcription"
    if any(k in search_text for k in ['nucleosome', 'histone', 'chromatin', 'remodeling']):
        return "Chromatin/Compaction"
    if any(k in search_text for k in ['drug', 'intercalator', 'binding']):
        return "Drug/Ligand Binding"
    
    return "Unknown"

def fetch_and_process_annotations(pdb_id):
    """
    Fetches data from the RCSB API for a single PDB and processes it.
    """
    query = get_pdb_graphql_query(pdb_id)
    response = requests.post("https://data.rcsb.org/graphql", json={'query': query})
    
    if response.status_code != 200:
        return {"error": f"API query failed with status {response.status_code}"}

    data = response.json()
    if not data.get('data') or not data['data'].get('entry'):
        return {"error": "No data in API response"}

    entry = data['data']['entry']
    
    # --- Initialize results ---
    annotations = {
        'Organism': set(),
        'Function': "Unknown",
        'GC_Content': None,
        'CpG_Count': None,
        'DNA_Sequence': None
    }

    # --- Extract keywords ---
    keywords = entry.get('struct_keywords', {}).get('pdbx_keywords', '')

    # --- Process polymer entities ---
    dna_entity_found = False
    for entity in entry.get('polymer_entities', []):
        poly_type = entity.get('entity_poly', {}).get('type', '')
        
        # --- Q1: Sequence Analysis (for DNA) ---
        if 'DNA' in poly_type and not dna_entity_found:
            sequence = entity.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can', '')
            if sequence:
                gc_content, cpg_count = analyze_sequence(sequence)
                annotations['GC_Content'] = round(gc_content, 3)
                annotations['CpG_Count'] = cpg_count
                annotations['DNA_Sequence'] = sequence.replace("\n", "")
                dna_entity_found = True # Process only the first DNA entity for simplicity

        # --- Q2: Species Conservation ---
        # Try to get specific NCBI name first
        source_org_data = entity.get('rcsb_entity_source_organism')
        if source_org_data:
            if not isinstance(source_org_data, list):
                source_org_data = [source_org_data]
            for org in source_org_data:
                if org and org.get('ncbi_scientific_name'):
                    annotations['Organism'].add(org['ncbi_scientific_name'])
        
        # Fallback to general source organism names from container identifiers
        container_ids = entity.get('rcsb_polymer_entity_container_identifiers', {})
        source_names = container_ids.get('source_organism_names')
        if source_names:
            for name in source_names:
                if name:
                    annotations['Organism'].add(name)
        
        # --- Q3: Functional Correlation (from proteins) ---
        if 'protein' in poly_type and annotations['Function'] == 'Unknown':
            annotations['Function'] = categorize_function(entity, keywords)

    # Finalize organism list
    annotations['Organism'] = ", ".join(sorted(list(annotations['Organism']))) if annotations['Organism'] else "N/A"
    
    return annotations

def main():
    """
    Main execution function.
    """
    parser = argparse.ArgumentParser(
        description="Enrich SHT outlier report with biological annotations from RCSB PDB."
    )
    parser.add_argument(
        '--outlier_csv', type=Path, required=True,
        help="Path to the outlier report CSV file (e.g., 'sht_outlier_analysis/outlier_report.csv')."
    )
    parser.add_argument(
        '--output_dir', type=Path, default=Path('sht_annotation_analysis'),
        help="Directory to save the enriched analysis report."
    )
    args = parser.parse_args()

    # --- Setup ---
    args.output_dir.mkdir(exist_ok=True)
    if not args.outlier_csv.exists():
        print(f"Error: Outlier report not found at '{args.outlier_csv}'")
        return

    print(f"--- Starting Biological Annotation of {args.outlier_csv.name} ---")
    outlier_df = pd.read_csv(args.outlier_csv)
    
    # --- API Fetching Loop ---
    all_annotations = []
    print(f"Fetching annotations for {len(outlier_df)} outlier PDBs...")
    
    # Use tqdm for a progress bar
    for index, row in tqdm(outlier_df.iterrows(), total=outlier_df.shape[0]):
        pdb_id = row['PDB ID']
        try:
            annotations = fetch_and_process_annotations(pdb_id)
            annotations['PDB ID'] = pdb_id
            all_annotations.append(annotations)
        except Exception as e:
            print(f"Failed to process {pdb_id}: {e}")
            all_annotations.append({'PDB ID': pdb_id, 'error': str(e)})
        
        # Be polite to the API server
        time.sleep(0.1)

    # --- Merge and Save Results ---
    print("\n--- Merging annotations and saving final report ---")
    
    # Convert list of dicts to DataFrame
    annotation_df = pd.DataFrame(all_annotations)
    
    # Merge with the original outlier data
    enriched_df = pd.merge(outlier_df, annotation_df, on='PDB ID', how='left')
    
    # Ensure essential columns for summary exist and fill NaNs from failed API calls
    for col in ['Function', 'Organism']:
        if col not in enriched_df.columns:
            enriched_df[col] = "N/A"
        else:
            enriched_df[col].fillna('N/A', inplace=True)

    # Reorder columns for better readability
    column_order = [
        'PDB ID', 'Cluster', 'DNA Type', 'Characterization', 'Function',
        'Organism', 'GC_Content', 'CpG_Count', 'Title', 'DNA_Sequence'
    ]
    # Ensure all columns exist before reordering
    final_columns = [col for col in column_order if col in enriched_df.columns]
    enriched_df = enriched_df[final_columns]

    # Save the final report
    output_path = args.output_dir / 'sht_enriched_outlier_report.csv'
    enriched_df.to_csv(output_path, index=False)

    print(f"\nAnalysis complete. Enriched report saved to '{output_path}'")
    
    # --- Display a summary of the findings ---
    print("\n--- Summary of Functional Correlations ---")
    # Filter by the three main discovery categories
    replication_signatures = enriched_df[enriched_df['Characterization'].str.contains("oscillation", na=False)]
    underwound_signatures = enriched_df[enriched_df['Characterization'].str.contains("Hyper-underwound", na=False)]
    overwound_signatures = enriched_df[enriched_df['Characterization'].str.contains("Hyper-overwound", na=False)]

    print("\n[Discovery 1: High-Frequency Oscillation]")
    if not replication_signatures.empty:
        print(replication_signatures[['Function', 'Organism']].value_counts().to_frame('Count').head())

    print("\n[Discovery 2: Hyper-Underwound States]")
    if not underwound_signatures.empty:
        print(underwound_signatures[['Function', 'Organism']].value_counts().to_frame('Count').head())
    
    print("\n[Discovery 3: Hyper-Overwound States]")
    if not overwound_signatures.empty:
        print(overwound_signatures[['Function', 'PDB ID', 'Title']].head())

if __name__ == "__main__":
    main()
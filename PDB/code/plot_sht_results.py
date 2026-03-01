#!/usr/bin/env python3
"""
plot_sht_results.py

Generate a variety of plots to visualize the results from the batch SHT analysis.
This script creates:
1. A 2D "shape space" scatter plot of Twist vs. Writhe.
2. A histogram of the global supercoiling density (sigma).
3. Box plots comparing sigma across different DNA categories.
"""

import argparse
import json
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# --- Plotting Style ---
sns.set_theme(style="whitegrid")

def classify_dna_type(row):
    """
    Classify a PDB entry into a category based on its metadata.
    This uses simple keyword matching on the PDB title.
    """
    title = str(row.get('title', '')).lower()
    pdb_id = str(row.get('pdb_id', '')).lower()

    # Handle canonical forms first
    if pdb_id == '1bna':
        return 'B-DNA (canonical)'
    if pdb_id == '1ana':
        return 'A-DNA (canonical)'
    if 'z-dna' in title or pdb_id in ['1dcg', '131d']:
        return 'Z-DNA'
    if 'a-dna' in title:
        return 'A-DNA'

    # Identify major complex types
    if 'nucleosome' in title or pdb_id == '1aoi':
        return 'Nucleosome'
    protein_keywords = [
        'protein', 'complex', 'factor', 'polymerase', 'repressor',
        'nuclease', 'domain', 'transferase', 'recombinase', 'helicase'
    ]
    if any(keyword in title for keyword in protein_keywords):
        return 'Protein-Complex'

    # Identify drug complexes
    drug_keywords = ['drug', 'netropsin', 'hoechst', 'cisplatin', 'daunomycin']
    if any(keyword in title for keyword in drug_keywords):
        return 'Drug-Complex'

    # Fallback categories
    if 'x-ray' in str(row.get('method', '')).lower():
        return 'B-DNA (other)'

    return 'Other'


def plot_twist_vs_writhe(df, output_dir):
    """
    Generate a 2D scatter plot of SHT Twist vs. SHT Writhe, color-coded by DNA type.
    This visualizes the "shape space" of the analyzed DNA structures.
    """
    print("Generating Twist vs. Writhe shape space plot...")
    plt.figure(figsize=(14, 10))

    # Define a specific color palette for clarity
    unique_types = df['dna_type'].unique()
    palette = sns.color_palette("husl", len(unique_types))
    type_palette = dict(zip(unique_types, palette))
    # Assign specific, distinct colors to key types
    if 'B-DNA (canonical)' in type_palette: type_palette['B-DNA (canonical)'] = 'blue'
    if 'A-DNA (canonical)' in type_palette: type_palette['A-DNA (canonical)'] = 'green'
    if 'Z-DNA' in type_palette: type_palette['Z-DNA'] = 'red'
    if 'Protein-Complex' in type_palette: type_palette['Protein-Complex'] = 'purple'
    if 'Nucleosome' in type_palette: type_palette['Nucleosome'] = 'black'


    ax = sns.scatterplot(
        data=df,
        x='Wr_SHT',
        y='Tw_SHT',
        hue='dna_type',
        palette=type_palette,
        s=50,
        alpha=0.7,
        edgecolor='k',
        linewidth=0.5
    )

    plt.title(r'DNA Shape Space: SHT Twist vs. Writhe', fontsize=18, weight='bold')
    plt.xlabel(r'SHT Writhe ($Wr_{SHT}$)', fontsize=14)
    plt.ylabel(r'SHT Twist ($Tw_{SHT}$)', fontsize=14)
    plt.axhline(0, color='grey', linestyle='--', linewidth=0.8)
    plt.axvline(0, color='grey', linestyle='--', linewidth=0.8)
    plt.legend(title='DNA Type', bbox_to_anchor=(1.02, 1), loc='upper left')

    # Annotate key structures for context
    annotations = {
        '1BNA': 'B-DNA',
        '1ANA': 'A-DNA',
        '1DCG': 'Z-DNA',
        '1AOI': 'Nucleosome',
        '1A2E': 'Cisplatin Adduct'
    }
    for pdb_id, label in annotations.items():
        entry = df[df['pdb_id'] == pdb_id]
        if not entry.empty:
            row = entry.iloc[0]
            ax.text(row['Wr_SHT'] + 0.01, row['Tw_SHT'], label,
                    fontsize=9, ha='left', va='center', weight='bold',
                    bbox=dict(boxstyle="round,pad=0.3", fc="yellow", ec="k", lw=1, alpha=0.8))

    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout(rect=[0, 0, 0.85, 1]) # Adjust layout to make space for legend
    output_path = output_dir / "sht_shape_space_Tw_vs_Wr.png"
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"  -> Saved to {output_path}")


def plot_sigma_distribution(df, output_dir):
    """
    Generate a histogram and a boxplot for the global supercoiling density (sigma_SHT_global).
    """
    print("Generating sigma distribution plots...")
    # 1. Histogram
    palette = sns.color_palette("husl", len(df['dna_type'].unique()))
    plt.figure(figsize=(12, 7))
    sns.histplot(data=df, x='sigma_SHT_global', hue='dna_type', multiple="stack", bins=50, palette=palette)
    plt.title(r'Distribution of Global Supercoiling Density ($\sigma_{SHT}$)', fontsize=16, weight='bold')
    plt.xlabel(r'Global Supercoiling Density ($\sigma_{SHT}$)', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.xlim(-2, 1) # Zoom in on the most populated region
    output_path_hist = output_dir / "sht_sigma_histogram.png"
    plt.savefig(output_path_hist, dpi=300)
    plt.close()
    print(f"  -> Saved histogram to {output_path_hist}")

    # 2. Box Plot
    # Order categories for better visualization
    order = sorted(df['dna_type'].unique())
    plt.figure(figsize=(14, 8))
    sns.boxplot(data=df, x='sigma_SHT_global', y='dna_type', order=order, palette=palette)
    plt.title(r'Global Supercoiling Density by DNA Type', fontsize=16, weight='bold')
    plt.xlabel(r'Global Supercoiling Density ($\sigma_{SHT}$)', fontsize=12)
    plt.ylabel('DNA Type', fontsize=12)
    plt.xlim(-2, 0.5) # Zoom in on the most populated region
    plt.tight_layout()
    output_path_box = output_dir / "sht_sigma_boxplot.png"
    plt.savefig(output_path_box, dpi=300)
    plt.close()
    print(f"  -> Saved boxplot to {output_path_box}")


def main():
    parser = argparse.ArgumentParser(
        description="Visualize results from SHT batch analysis."
    )
    parser.add_argument(
        '--results_csv',
        type=Path,
        # Default path points to the output location of batch_analyze.py
        default=Path(__file__).parent.parent / 'sht_analysis_output' / 'sht_analysis_results.csv',
        help="Path to the aggregated results CSV file."
    )
    parser.add_argument(
        '--metadata_file',
        type=Path,
        # Default path points to the output location of filter_and_prepare.py
        default=Path(__file__).parent.parent / 'dataset' / 'hard_dataset' / 'curated_nonredundant_duplex_metadata.json',
        help="Path to the curated metadata JSON file from filter_and_prepare.py"
    )
    parser.add_argument(
        '--output_dir',
        type=Path,
        # Save plots to a dedicated folder within the main analysis output directory
        default=Path(__file__).parent.parent / 'sht_analysis_output' / 'summary_plots',
        help="Directory to save the generated plots."
    )
    args = parser.parse_args()

    # --- Validation and Setup ---
    if not args.results_csv.exists():
        print(f"Error: Results CSV file not found at {args.results_csv}")
        return
    if not args.metadata_file.exists():
        print(f"Error: Metadata JSON file not found at {args.metadata_file}")
        return

    args.output_dir.mkdir(exist_ok=True)

    # --- Data Loading and Merging ---
    print("Loading and processing data...")
    results_df = pd.read_csv(args.results_csv)
    with open(args.metadata_file, 'r') as f:
        metadata = json.load(f)

    meta_df = pd.DataFrame.from_dict(metadata, orient='index')
    meta_df = meta_df.reset_index().rename(columns={'index': 'pdb_id'})

    # Merge results with metadata
    full_df = pd.merge(results_df, meta_df, on='pdb_id', how='left')

    # Classify each structure
    full_df['dna_type'] = full_df.apply(classify_dna_type, axis=1)

    # Remove rows with NaN values in key columns to avoid plotting errors
    key_cols = ['Tw_SHT', 'Wr_SHT', 'sigma_SHT_global']
    original_rows = len(full_df)
    full_df.dropna(subset=key_cols, inplace=True)
    if len(full_df) < original_rows:
        print(f"  Dropped {original_rows - len(full_df)} rows with missing values.")

    print(f"Data processed. Found {len(full_df)} valid entries.")
    print("DNA types identified:\n", full_df['dna_type'].value_counts())

    # --- Generate Plots ---
    plot_twist_vs_writhe(full_df, args.output_dir)
    plot_sigma_distribution(full_df, args.output_dir)

    print("\nPlotting complete.")


if __name__ == "__main__":
    main()

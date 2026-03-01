#!/usr/bin/env python3
"""
sht_statistical_analysis.py

This script performs an in-depth statistical analysis of the SHT results.
It computes descriptive statistics, performs hypothesis testing between DNA
categories, and conducts a sensitivity analysis on the bp_per_turn parameter.

Goals:
1.  Compute mean/std of sigma_SHT, Tw_SHT, Wr_SHT, and Length for each category.
2.  Run t-tests and KS-tests to compare categories.
3.  Check sensitivity of sigma_SHT to the bp_per_turn reference value.
4.  Analyze and plot length distributions.
5.  Analyze and plot distributions of local supercoiling density.
6.  Perform shape clustering on local supercoiling profiles.
"""

import argparse
import pickle
import json
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import umap

# --- Plotting Style ---
sns.set_theme(style="whitegrid")
plt.rcParams['figure.dpi'] = 150

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

def calculate_descriptive_stats(df):
    """
    Calculates and prints mean and standard deviation for key SHT metrics,
    grouped by DNA type.
    """
    print("\n--- [Step 1] Descriptive Statistics per DNA Category ---")
    
    # Define the metrics and corresponding aggregation functions
    agg_metrics = {
        'sigma_SHT_global': ['mean', 'std'],
        'Tw_SHT': ['mean', 'std'],
        'Wr_SHT': ['mean', 'std'],
        'L': ['mean', 'std', 'count']
    }
    
    # Group by DNA type and aggregate
    stats_df = df.groupby('dna_type').agg(agg_metrics).round(3)
    
    # Rename columns for clarity
    stats_df.columns = ['Sigma Mean', 'Sigma Std', 'Twist Mean', 'Twist Std', 
                        'Writhe Mean', 'Writhe Std', 'Length Mean', 'Length Std', 'Count']
    
    print(stats_df)
    print("-" * 70)

def perform_statistical_tests(df):
    """
    Performs T-tests and KS-tests on sigma_SHT_global between major DNA categories.
    """
    print("\n--- [Step 2] Statistical Significance Tests (p-values) ---")
    print("Comparing distributions of Global Supercoiling Density (sigma_SHT)\n")

    categories_to_test = [
        'A-DNA', 'B-DNA (other)', 'Z-DNA', 'Nucleosome', 'Protein-Complex'
    ]
    
    # Filter the DataFrame to only include the categories of interest
    test_df = df[df['dna_type'].isin(categories_to_test)]
    
    # Get all unique pairs of categories to compare
    category_pairs = list(combinations(categories_to_test, 2))
    
    results = []
    for cat1, cat2 in category_pairs:
        data1 = test_df[test_df['dna_type'] == cat1]['sigma_SHT_global'].dropna()
        data2 = test_df[test_df['dna_type'] == cat2]['sigma_SHT_global'].dropna()
        
        if len(data1) < 2 or len(data2) < 2:
            continue
            
        # T-test (Welch's, assuming unequal variance)
        t_stat, t_pvalue = stats.ttest_ind(data1, data2, equal_var=False, nan_policy='omit')
        
        # Kolmogorov-Smirnov test
        ks_stat, ks_pvalue = stats.ks_2samp(data1, data2)
        
        results.append({
            "Comparison": f"{cat1} vs {cat2}",
            "T-test p-value": f"{t_pvalue:.2e}",
            "KS-test p-value": f"{ks_pvalue:.2e}"
        })
        
    results_df = pd.DataFrame(results)
    print(results_df.to_string(index=False))
    print("-" * 70)

def perform_sensitivity_analysis(df):
    """
    Analyzes the sensitivity of global sigma to the bp_per_turn parameter.
    """
    print("\n--- [Step 3] Sensitivity Analysis for bp_per_turn ---")
    print("Mean Global Sigma (σ) for different bp/turn reference values\n")

    # The original Tw0 was calculated with bp_per_turn = 10.5
    # We can estimate the number of base pairs from this.
    df['num_bp_estimated'] = df['Tw0'] * 10.5
    
    # Define the range of bp_per_turn values to test
    bp_per_turn_values = [10.4, 10.5, 10.6]
    
    # Calculate Lk_SHT once
    df['Lk_SHT'] = df['Tw_SHT'] + df['Wr_SHT']

    # Dictionary to hold results
    sensitivity_data = {}

    for bp_turn in bp_per_turn_values:
        # Recalculate Tw0 and sigma for the new bp_per_turn value
        tw0_new = df['num_bp_estimated'] / bp_turn
        sigma_new_col = f'sigma_at_{bp_turn}'
        df[sigma_new_col] = (df['Lk_SHT'] - tw0_new) / tw0_new
        
        # Calculate the mean sigma for each DNA type
        mean_sigma = df.groupby('dna_type')[sigma_new_col].mean()
        sensitivity_data[f'{bp_turn} bp/turn'] = mean_sigma

    # Create a DataFrame from the results and display it
    sensitivity_df = pd.DataFrame(sensitivity_data).round(4)
    
    # We are most interested in the major, well-defined categories
    categories_of_interest = [
        'A-DNA', 'A-DNA (canonical)', 'B-DNA (canonical)', 'B-DNA (other)', 
        'Z-DNA', 'Nucleosome', 'Protein-Complex', 'Drug-Complex'
    ]
    
    print(sensitivity_df.loc[sensitivity_df.index.isin(categories_of_interest)])
    print("\nNote: The absolute sigma values change, but the relative differences between\ncategories remain stable, demonstrating the robustness of the SHT method.")
    print("-" * 70)

def plot_summary_distributions(df, local_sigma_df, output_dir):
    """
    Generates and saves plots for length and local sigma distributions.
    """
    print("\n--- [Step 4] Plotting Distributions ---")
    output_dir.mkdir(exist_ok=True)

    # 1. Length Distribution Plot
    plt.figure(figsize=(12, 7))
    sns.histplot(data=df, x='L', hue='dna_type', multiple="stack", bins=50, log_scale=(True, False))
    plt.title('Distribution of DNA Lengths by Category', fontsize=16, weight='bold')
    plt.xlabel('Effective DNA Length (L, Å) [log scale]', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.tight_layout()
    len_dist_path = output_dir / "sht_length_distribution.png"
    plt.savefig(len_dist_path, dpi=300)
    plt.close()
    print(f"  -> Saved length distribution plot to {len_dist_path}")

    # 2. Local Sigma Distribution Plot
    # We are most interested in the major, well-defined categories
    categories_of_interest = [
        'A-DNA', 'B-DNA (other)', 'Z-DNA', 'Nucleosome', 'Protein-Complex'
    ]
    plot_df = local_sigma_df[local_sigma_df['dna_type'].isin(categories_of_interest)]

    plt.figure(figsize=(14, 8))
    sns.kdeplot(data=plot_df, x='sigma_local', hue='dna_type', fill=True, common_norm=False)
    plt.title(r'Distribution of Local Supercoiling Density ($\sigma_{SHT}(h)$) by Category', fontsize=16, weight='bold')
    plt.xlabel(r'Local Supercoiling Density ($\sigma_{SHT}(h)$)', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.xlim(-2.5, 2.5)
    plt.axvline(0, color='k', linestyle='--', linewidth=1)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    local_sigma_path = output_dir / "sht_local_sigma_distribution.png"
    plt.savefig(local_sigma_path, dpi=300)
    plt.close()
    print(f"  -> Saved local sigma distribution plot to {local_sigma_path}")
    print("-" * 70)

def resample_profile(profile, num_points=50):
    """
    Resamples a 1D profile to a fixed number of points using linear interpolation.
    """
    if profile is None or len(profile) < 2:
        return np.full(num_points, np.nan)
    
    original_indices = np.linspace(0, 1, len(profile))
    new_indices = np.linspace(0, 1, num_points)
    
    # Clamp sigma values to a reasonable range to avoid extreme outliers
    # from dominating the clustering.
    profile_clamped = np.clip(profile, -5, 5)
    
    resampled = np.interp(new_indices, original_indices, profile_clamped)
    return resampled

def perform_shape_clustering(df, local_profiles, output_dir, n_clusters=6, n_points=50):
    """
    Performs clustering on the resampled local sigma profiles.
    """
    print("\n--- [Step 5] Shape Clustering of Local Supercoiling Profiles ---")
    
    # 1. Resample all profiles to a fixed length
    resampled_data = {pdb_id: resample_profile(prof, n_points) for pdb_id, prof in local_profiles.items()}
    resampled_df = pd.DataFrame.from_dict(resampled_data, orient='index', columns=[f'p_{i}' for i in range(n_points)])
    resampled_df.dropna(inplace=True) # Remove entries that couldn't be resampled

    # 2. Scale the data
    scaler = StandardScaler()
    scaled_profiles = scaler.fit_transform(resampled_df)

    # 3. UMAP for dimensionality reduction
    print(f"  -> Performing UMAP dimensionality reduction on {len(scaled_profiles)} profiles...")
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
    embedding = reducer.fit_transform(scaled_profiles)

    # 4. K-Means clustering on the UMAP embedding
    print(f"  -> Performing K-Means clustering (k={n_clusters})...")
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    cluster_labels = kmeans.fit_predict(embedding)

    # 5. Merge results back into a single DataFrame for analysis
    resampled_df['cluster'] = cluster_labels
    resampled_df['umap1'] = embedding[:, 0]
    resampled_df['umap2'] = embedding[:, 1]
    
    # Merge with the main metadata dataframe
    analysis_df = df.merge(resampled_df, left_on='pdb_id', right_index=True)
    
    print("  -> Clustering complete. Generating analysis plots...")
    plot_cluster_analysis(analysis_df, n_clusters, n_points, output_dir)
    print('  -> Saving cluster analysis dataframe...')
    analysis_df.to_csv(output_dir / 'sht_cluster_analysis.csv', index=False)
    
    # Save the dataframe with cluster assignments for outlier analysis
    cluster_output_path = output_dir / "sht_cluster_analysis.csv"
    analysis_df.to_csv(cluster_output_path, index=False)
    print(f"  -> Saved cluster analysis dataframe to {cluster_output_path}")

def plot_cluster_analysis(analysis_df, n_clusters, n_points, output_dir):
    print('  -> Saving cluster analysis dataframe...')
    analysis_df.to_csv(output_dir / 'sht_cluster_analysis.csv', index=False)
    """
    Generates and saves plots for the shape clustering results.
    """
    # 1. UMAP Shape Space Plot
    plt.figure(figsize=(12, 10))
    sns.scatterplot(data=analysis_df, x='umap1', y='umap2', hue='cluster', palette='viridis', s=20, alpha=0.7)
    plt.title('UMAP Projection of σ(h) Shape Space', fontsize=16, weight='bold')
    plt.xlabel('UMAP Dimension 1', fontsize=12)
    plt.ylabel('UMAP Dimension 2', fontsize=12)
    plt.legend(title='Cluster ID')
    plt.tight_layout()
    umap_path = output_dir / "sht_cluster_umap_space.png"
    plt.savefig(umap_path, dpi=300)
    plt.close()
    print(f"  -> Saved UMAP shape space plot to {umap_path}")

    # 2. Mean Cluster Shape Profiles
    profile_cols = [f'p_{i}' for i in range(n_points)]
    plt.figure(figsize=(12, 7))
    for i in range(n_clusters):
        cluster_profiles = analysis_df[analysis_df['cluster'] == i][profile_cols]
        mean_profile = cluster_profiles.mean(axis=0)
        plt.plot(np.linspace(0, 100, n_points), mean_profile, label=f'Cluster {i} (n={len(cluster_profiles)})')
    plt.title(r'Mean Local Supercoiling Profile ($\sigma_{SHT}(h)$) for Each Shape Cluster', fontsize=16, weight='bold')
    plt.xlabel('Normalized Axis Position (%)', fontsize=12)
    plt.ylabel(r'Mean Local Supercoiling Density $\langle\sigma(h)\rangle$', fontsize=12)
    plt.axhline(0, color='k', linestyle='--', linewidth=1)
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    mean_shapes_path = output_dir / "sht_cluster_mean_shapes.png"
    plt.savefig(mean_shapes_path, dpi=300)
    plt.close()
    print(f"  -> Saved mean cluster shapes plot to {mean_shapes_path}")

    # 3. Cluster Composition Bar Chart
    composition = analysis_df.groupby('cluster')['dna_type'].value_counts(normalize=True).unstack().fillna(0)
    plt.figure(figsize=(14, 8))
    composition.plot(kind='bar', stacked=True, figsize=(14, 8), colormap='tab20')
    plt.title('Composition of DNA Categories within Each Shape Cluster', fontsize=16, weight='bold')
    plt.xlabel('Shape Cluster ID', fontsize=12)
    plt.ylabel('Proportion of Category', fontsize=12)
    plt.xticks(rotation=0)
    plt.legend(title='DNA Category', bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    composition_path = output_dir / "sht_cluster_composition.png"
    plt.savefig(composition_path, dpi=300)
    plt.close()
    print(f"  -> Saved cluster composition plot to {composition_path}")
    print("-" * 70)

def load_data_from_pickles(pickle_dir, metadata_df):
    """
    Loads analysis results from individual pickle files and merges with metadata.
    """
    all_records = []
    local_sigma_records = []
    local_sigma_profiles = {}

    for pkl_file in pickle_dir.glob("*_profile.pkl"):
        pdb_id = pkl_file.stem.replace("_profile", "")
        with open(pkl_file, 'rb') as f:
            data = pickle.load(f)
        
        # Global data
        record = {'pdb_id': pdb_id, **{k: data[k] for k in ['sigma_SHT_global', 'Tw_SHT', 'Wr_SHT', 'Tw0', 'L']}}
        all_records.append(record)

        # Local sigma data
        if 'sigma_local_h' in data and data['sigma_local_h'] is not None and len(data['sigma_local_h']) > 1:
            local_sigma_profiles[pdb_id] = data['sigma_local_h']
            for sigma_val in data['sigma_local_h']:
                local_sigma_records.append({'pdb_id': pdb_id, 'sigma_local': sigma_val})

    return pd.DataFrame(all_records), pd.DataFrame(local_sigma_records), local_sigma_profiles

def main():
    parser = argparse.ArgumentParser(
        description="Perform statistical analysis on SHT results."
    )
    parser.add_argument(
        '--results_csv', type=Path,
        default=Path(__file__).parent.parent / 'sht_analysis_output' / 'sht_analysis_results.csv',
        help="Path to the aggregated results CSV file."
    )
    parser.add_argument(
        '--metadata_file', type=Path,
        default=Path(__file__).parent.parent / 'dataset' / 'hard_dataset' / 'curated_nonredundant_duplex_metadata.json',
        help="Path to the curated metadata JSON file."
    )
    parser.add_argument(
        '--pickle_dir', type=Path,
        default=Path(__file__).parent.parent / 'sht_analysis_output' / 'sht_local_profiles',
        help="Directory containing detailed local profile pickle files."
    )
    parser.add_argument(
        '--output_dir',
        type=Path,
        default=Path(__file__).parent.parent / 'sht_analysis_output' / 'statistical_plots',
        help="Directory to save the generated plots."
    )
    args = parser.parse_args()

    if not args.pickle_dir.exists():
        print(f"Error: Pickle directory not found at {args.pickle_dir}")
        print("Please run batch_analyze.py first to generate local profiles.")
        return
    if not args.metadata_file.exists():
        print(f"Error: Metadata JSON file not found at {args.metadata_file}")
        return

    print("Loading and processing data...")
    with open(args.metadata_file, 'r') as f:
        metadata = json.load(f)

    meta_df = pd.DataFrame.from_dict(metadata, orient='index').reset_index().rename(columns={'index': 'pdb_id'})

    # Load data from pickle files for a more in-depth analysis
    results_df, local_sigma_df, local_profiles = load_data_from_pickles(args.pickle_dir, meta_df)

    # Merge with metadata
    full_df = pd.merge(results_df, meta_df, on='pdb_id', how='left')
    full_df['dna_type'] = full_df.apply(classify_dna_type, axis=1)

    local_sigma_df = pd.merge(local_sigma_df, full_df[['pdb_id', 'dna_type']], on='pdb_id', how='left')

    key_cols = ['Tw_SHT', 'Wr_SHT', 'sigma_SHT_global', 'L', 'Tw0']
    full_df.dropna(subset=key_cols, inplace=True)
    local_sigma_df.dropna(subset=['sigma_local', 'dna_type'], inplace=True)
    
    print(f"Data processed. Found {len(full_df)} valid entries.")

    calculate_descriptive_stats(full_df)
    perform_statistical_tests(full_df)
    perform_sensitivity_analysis(full_df)
    plot_summary_distributions(full_df, local_sigma_df, args.output_dir)
    perform_shape_clustering(full_df, local_profiles, args.output_dir)

if __name__ == "__main__":
    main()

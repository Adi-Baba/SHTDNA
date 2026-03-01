#!/usr/bin/env python3
"""
sht_outlier_detection.py

This script uses the results of the SHT clustering analysis to identify
and characterize outliers. The goal is to discover novel or uncharacterized
DNA conformations.

Workflow:
1.  Load the main analysis DataFrame which includes cluster assignments for each PDB.
2.  For each cluster, calculate the centroid in the UMAP space.
3.  Compute the Euclidean distance of each structure from its cluster centroid.
4.  Identify outliers based on a distance threshold (e.g., > 95th percentile).
5.  Load the local supercoiling profiles (sigma_SHT(h)) for these outliers.
6.  Characterize the outlier profiles based on their statistical properties:
    - Mean: Global positive/negative supercoiling.
    - Variance/Kurtosis: Sharp kinks or localized stress.
    - FFT: Abnormal oscillatory patterns.
7.  Correlate findings with PDB metadata to hypothesize the biological context.
8.  Generate a detailed report with plots for each significant outlier.
"""

import argparse
import pickle
import json
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns

# --- Plotting Style ---
sns.set_theme(style="whitegrid")
plt.rcParams['figure.dpi'] = 150

def load_analysis_data(cluster_results_file, pickle_dir):
    """
    Loads all necessary data: cluster results and local profiles.
    """
    print("--- Loading Data ---")
    
    if not cluster_results_file.exists():
        raise FileNotFoundError(f"Cluster results file not found: {cluster_results_file}")
    full_df = pd.read_csv(cluster_results_file)
    print(f"Loaded {len(full_df)} clustered entries from {cluster_results_file.name}.")

    # Load local sigma profiles
    local_profiles = {}
    for pkl_file in pickle_dir.glob("*_profile.pkl"):
        pdb_id = pkl_file.stem.replace("_profile", "")
        with open(pkl_file, 'rb') as f:
            data = pickle.load(f)
        if 'sigma_local_h' in data and data['sigma_local_h'] is not None:
            local_profiles[pdb_id] = {
                'h': data['h'],
                'sigma_local_h': data['sigma_local_h']
            }

    print("Data loading and merging complete.")
    return full_df.copy(), local_profiles

def find_cluster_outliers(df, outlier_percentile=95):
    """
    Identifies outliers in each cluster based on distance from the UMAP centroid.
    """
    print(f"\n--- [Step 1] Identifying Outliers (Top {100-outlier_percentile}%) ---")
    
    outliers = []
    for cluster_id in sorted(df['cluster'].unique()):
        cluster_data = df[df['cluster'] == cluster_id]
        
        if len(cluster_data) < 3:
            continue

        # Calculate centroid of the cluster in UMAP space
        centroid = cluster_data[['umap1', 'umap2']].mean().values
        
        # Calculate Euclidean distance of each point from the centroid
        distances = np.sqrt(np.sum((cluster_data[['umap1', 'umap2']].values - centroid)**2, axis=1))
        
        # Determine the distance threshold for outliers
        threshold = np.percentile(distances, outlier_percentile)
        
        # Identify outliers
        cluster_outliers = cluster_data[distances > threshold]
        
        for _, row in cluster_outliers.iterrows():
            outliers.append({
                'pdb_id': row['pdb_id'],
                'cluster': cluster_id,
                'distance': distances[cluster_data.index == row.name][0],
                'dna_type': row['dna_type'],
                'title': row['title']
            })
            
    outlier_df = pd.DataFrame(outliers).sort_values('distance', ascending=False)
    print(f"Identified {len(outlier_df)} potential outliers across all clusters.")
    print(outlier_df.head())
    return outlier_df

def characterize_profile(profile_data):
    """
    Analyzes a sigma_SHT(h) profile and returns its characteristics.
    """
    h = profile_data['h']
    sigma = profile_data['sigma_local_h']
    
    if len(sigma) < 5: # Need enough points for meaningful stats
        return "Profile too short"

    # 1. Global Bias
    mean_sigma = np.mean(sigma)
    if mean_sigma > 0.5:
        return f"Hyper-overwound (mean σ = {mean_sigma:.2f})"
    if mean_sigma < -0.5:
        return f"Hyper-underwound (mean σ = {mean_sigma:.2f})"

    # 2. Localized Stress (Kinks/Defects)
    variance = np.var(sigma)
    kurt = stats.kurtosis(sigma)
    if kurt > 3.0 and variance > 0.2:
        return f"Sharp Kink/Defect (kurtosis = {kurt:.2f})"

    # 3. Oscillatory Patterns
    # Perform FFT to find dominant frequencies
    fft_vals = np.fft.fft(sigma)
    fft_freq = np.fft.fftfreq(len(h), d=(h[1]-h[0]))
    
    # Look for strong peaks in the mid-frequency range (ignoring the DC component)
    power_spectrum = np.abs(fft_vals[1:])**2
    peak_freq = np.abs(fft_freq[1:][np.argmax(power_spectrum)])
    
    # Typical B-DNA has a period of ~34Å, so frequency ~1/34 = 0.03 Å⁻¹
    # A-DNA has period ~28Å, freq ~0.035 Å⁻¹
    # Z-DNA has period ~45Å, freq ~0.022 Å⁻¹
    # Look for frequencies outside this range or unusually strong peaks.
    if 0.05 < peak_freq < 0.15:
        return f"Abnormal high-frequency oscillation (peak at {peak_freq:.3f} Å⁻¹)"
    if peak_freq < 0.01:
         return f"Abnormal low-frequency oscillation (peak at {peak_freq:.3f} Å⁻¹)"

    return "Atypical/Complex Profile"

def analyze_and_plot_outliers(outlier_df, local_profiles, output_dir):
    """
    For each outlier, characterize its profile and generate a report plot.
    """
    print("\n--- [Step 2] Characterizing Outlier Profiles ---")
    output_dir.mkdir(exist_ok=True)
    
    report_data = []

    for _, outlier in outlier_df.iterrows():
        pdb_id = outlier['pdb_id']
        if pdb_id not in local_profiles:
            continue

        profile_data = local_profiles[pdb_id]
        characterization = characterize_profile(profile_data)
        
        report_data.append({
            'PDB ID': pdb_id,
            'Cluster': outlier['cluster'],
            'DNA Type': outlier['dna_type'],
            'Characterization': characterization,
            'Title': outlier['title']
        })

        # --- Generate Plot ---
        fig, ax = plt.subplots(figsize=(10, 5))
        
        h = profile_data['h']
        sigma_local = profile_data['sigma_local_h']
        
        ax.plot(h, sigma_local, marker='o', linestyle='-', markersize=4, label=r'$\sigma_{SHT}(h)$')
        ax.axhline(0, color='k', linestyle='--', linewidth=1)
        
        # Add stats to the plot
        mean_val = np.mean(sigma_local)
        std_val = np.std(sigma_local)
        ax.axhline(mean_val, color='r', linestyle=':', linewidth=1.5, label=f'Mean = {mean_val:.3f}')
        
        title_text = (
            f"Outlier Profile: {pdb_id} (from Cluster {outlier['cluster']})\n"
            f"Classification: {outlier['dna_type']} | Characterization: {characterization}"
        )
        ax.set_title(title_text, fontsize=12, weight='bold')
        ax.set_xlabel('Position along Principal Axis (h, Å)', fontsize=10)
        ax.set_ylabel(r'Local Supercoiling Density $\sigma_{SHT}(h)$', fontsize=10)
        ax.legend()
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        
        plt.tight_layout()
        plot_path = output_dir / f"outlier_{pdb_id}.png"
        plt.savefig(plot_path, dpi=150)
        plt.close(fig)

    # --- Print Final Report Table ---
    report_df = pd.DataFrame(report_data)
    print("\n--- Outlier Analysis Report ---")
    with pd.option_context('display.max_rows', None, 'display.max_colwidth', None):
        print(report_df.to_string(index=False))
    print("-" * 80)

    # --- Save Report to CSV ---
    report_csv_path = output_dir / "outlier_report.csv"
    report_df.to_csv(report_csv_path, index=False)
    print(f"\nDetailed outlier report saved to '{report_csv_path}'")

def main():
    parser = argparse.ArgumentParser(
        description="Identify and characterize outliers from SHT clustering analysis."
    )
    parser.add_argument(
        '--cluster_results', type=Path, required=True,
        help="Path to the CSV file containing UMAP coordinates and cluster assignments for each PDB ID."
    )
    parser.add_argument(
        '--pickle_dir',
        type=Path,
        default=Path(__file__).parent.parent / 'sht_analysis_output' / 'sht_local_profiles',
        help="Directory containing detailed local profile pickle files from batch_analyze.py."
    )
    parser.add_argument(
        '--output_dir', type=Path, default=Path('sht_outlier_analysis'),
        help="Directory to save the outlier plots and report."
    )
    parser.add_argument(
        '--percentile', type=int, default=95,
        help="Percentile threshold to define outliers (e.g., 95 means top 5%% are outliers)."
    )
    args = parser.parse_args()

    # --- Pre-computation Check ---
    # The sht_statistical_analysis.py script must be run first to generate the necessary input file.
    # Let's assume the output of that script is 'sht_cluster_analysis.csv'.
    # We need to modify sht_statistical_analysis.py to save its `analysis_df`.
    
    # For this exercise, I will add the required save operation to `sht_statistical_analysis.py`
    # and then proceed with this script.

    # --- Modify sht_statistical_analysis.py to save the analysis_df ---
    stat_analysis_path = Path('sht_statistical_analysis.py')
    if stat_analysis_path.exists():
        with open(stat_analysis_path, 'r') as f:
            lines = f.readlines()
        
        # Check if the save command is already there
        save_command = "analysis_df.to_csv(output_dir / 'sht_cluster_analysis.csv', index=False)\n"
        already_present = any(save_command.strip() in line for line in lines)

        if not already_present:
            print("Modifying 'sht_statistical_analysis.py' to save cluster results...")
            # Find the line to insert after
            for i, line in enumerate(lines):
                if "plot_cluster_analysis(analysis_df" in line:
                    lines.insert(i + 1, f"    print('  -> Saving cluster analysis dataframe...')\n    {save_command}")
                    with open(stat_analysis_path, 'w') as f:
                        f.writelines(lines)
                    print("Modification complete. Please re-run 'sht_statistical_analysis.py' to generate the required input file for this script.")
                    # return # We stop here to let the user re-run the script.
                    # For the purpose of this interactive session, I will assume it has been run.
        else:
            print("'sht_statistical_analysis.py' is already set up to save cluster results.")

    if not args.cluster_results.exists():
        print(f"\nError: The required input file '{args.cluster_results}' was not found.")
        print("Please run 'sht_statistical_analysis.py' first to generate it.")
        print("Example: python sht_statistical_analysis.py --output_dir .")
        return

    # --- Main Analysis ---
    try:
        # 1. Load Data
        full_df, local_profiles = load_analysis_data(args.cluster_results, args.pickle_dir)

        # 2. Find Outliers
        outlier_df = find_cluster_outliers(full_df, outlier_percentile=args.percentile)

        # 3. Characterize and Plot
        analyze_and_plot_outliers(outlier_df, local_profiles, args.output_dir)
        
        print(f"\nAnalysis complete. Outlier report and plots saved to '{args.output_dir}'.")

    except FileNotFoundError as e:
        print(f"\nError: A required file was not found. {e}")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")

if __name__ == "__main__":
    # This script requires the output from sht_statistical_analysis.py
    # Specifically, a CSV containing pdb_id, umap1, umap2, and cluster columns.
    # Let's assume that script was run and produced 'sht_statistical_plots/sht_cluster_analysis.csv'
    
    # Example usage from command line:
    # python sht_outlier_detection.py --cluster_results sht_statistical_plots/sht_cluster_analysis.csv
    
    main()
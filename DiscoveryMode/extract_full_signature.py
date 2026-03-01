
import sys
import os
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import pearsonr, spearmanr

# Add legacy code path
sys.path.append(r"d:\OnlyHST\SHTDNA\PDB\code")
import sht_dna_analysis

def extract_slice_data(pdb_path):
    try:
        results = sht_dna_analysis.compute_sigma_sht_from_pdb(
            str(pdb_path), num_slices=100
        )
        
        # Get vectors
        kappa = results.get('kappa_h', []) # v0
        epsilon = results.get('epsilon_h', []) # eta_r
        omega = results.get('omega_h', []) # Omega
        
        # Ensure they are the same length and valid
        n = len(kappa)
        if len(epsilon) != n or len(omega) != n:
            return None
            
        # Stack into (N, 3) matrix
        stack = np.column_stack((kappa, epsilon, omega))
        
        # Filter NaNs
        valid_mask = ~np.isnan(stack).any(axis=1)
        valid_stack = stack[valid_mask]
        
        return valid_stack

    except Exception as e:
        # print(f"Error processing {pdb_path}: {e}")
        return None

def main():
    dataset_dir = Path(r"d:\OnlyHST\SHTDNA\PDB\dataset\hard_dataset\structures")
    output_dir = Path(r"d:\OnlyHST\SHTDNA\DiscoveryMode")
    output_csv = output_dir / "mechanical_state_cloud.csv"
    output_log = output_dir / "discovery_log.txt"
    
    all_pdbs = list(dataset_dir.glob("*.cif")) + list(dataset_dir.glob("*.pdb"))
    
    # The "Cloud": A massive list of slice points
    # Columns: [Kappa(v0), Eta_r, Omega]
    cloud_data = []
    
    count = 0
    limit = 200 
    
    print(f"Starting Discovery Extraction (Full Signature) on {limit} structures...")
    
    for i, pdb_file in enumerate(all_pdbs):
        if count >= limit: break
            
        slice_matrix = extract_slice_data(pdb_file)
        if slice_matrix is not None and slice_matrix.size > 0:
            cloud_data.append(slice_matrix)
            count += 1
            if count % 10 == 0: print(f".", end='', flush=True)
            
    print("\nExtraction Complete. Assembling the Cloud...")
    
    # Concatenate all slice matrices
    if not cloud_data:
        print("No data extracted.")
        return

    full_cloud = np.vstack(cloud_data)
    df = pd.DataFrame(full_cloud, columns=['v0_Curvature', 'Eta_r_Shambhian', 'Omega_Mamton'])
    
    # Save raw data
    df.to_csv(output_csv, index=False)
    print(f"Cloud Size: {len(df)} slices.")
    
    # --- perform Correlation Analysis ---
    print("Computing Correlations...")
    
    # Pearson (Linear)
    corr_pearson, p_pearson = pearsonr(df['v0_Curvature'], df['Eta_r_Shambhian'])
    
    # Spearman (Rank/Non-linear)
    corr_spearman, p_spearman = spearmanr(df['v0_Curvature'], df['Eta_r_Shambhian'])
    
    log_content = []
    log_content.append("# SHT Discovery Mode: Correlation Analysis Log\n")
    log_content.append(f"Total Slices Analysis: {len(df)}\n")
    log_content.append("-" * 30 + "\n")
    
    log_content.append("## 1. Global Correlations (v0 vs Eta_r)\n")
    log_content.append(f"Pearson R:  {corr_pearson:.5f} (P-value: {p_pearson:.5e})\n")
    log_content.append(f"Spearman R: {corr_spearman:.5f} (P-value: {p_spearman:.5e})\n")
    if p_spearman < 1e-9:
        log_content.append("-> Result is STATISTICALLY SIGNIFICANT (P ~ 0).\n")
    else:
        log_content.append("-> Result may be noise.\n")
    log_content.append("\n")

    # --- Binning Analysis (Concrete Evidence) ---
    log_content.append("## 2. Binning Test (Concrete Evidence)\n")
    
    # Split into Quartiles based on Flatness (Eta_r)
    q1 = df['Eta_r_Shambhian'].quantile(0.25)
    q4 = df['Eta_r_Shambhian'].quantile(0.75)
    
    group_round = df[df['Eta_r_Shambhian'] <= q1]['v0_Curvature']
    group_flat  = df[df['Eta_r_Shambhian'] >= q4]['v0_Curvature']
    
    mean_round = group_round.mean()
    mean_flat = group_flat.mean()
    
    from scipy.stats import ttest_ind
    t_stat, p_ttest = ttest_ind(group_round, group_flat, equal_var=False)
    
    log_content.append(f"**Group 1: Round DNA** (Eta_r <= {q1:.3f})\n")
    log_content.append(f"   Mean Stiffness (v0): {mean_round:.5f}\n")
    log_content.append(f"**Group 2: Flat DNA** (Eta_r >= {q4:.3f})\n")
    log_content.append(f"   Mean Stiffness (v0): {mean_flat:.5f}\n")
    log_content.append(f"\nDifference: {mean_flat - mean_round:.5f}\n")
    log_content.append(f"T-Test P-value: {p_ttest:.5e}\n")
    
    if p_ttest < 0.05:
         log_content.append("CONCLUSION: The Stiffness-Flatness Relationship is REAL (Distinct Populations).\n")
    
    # --- Visualization ---
    try:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10, 6))
        plt.scatter(df['Eta_r_Shambhian'], df['v0_Curvature'], alpha=0.3, s=10, c='blue', label='DNA Slices')
        
        # Add trendline
        z = np.polyfit(df['Eta_r_Shambhian'], df['v0_Curvature'], 1)
        p = np.poly1d(z)
        plt.plot(df['Eta_r_Shambhian'], p(df['Eta_r_Shambhian']), "r--", label=f"Trend (R={corr_spearman:.2f})")
        
        plt.xlabel('Shambhian Twist ($\eta_r$)', fontsize=12)
        plt.ylabel('Curvature / Bending ($v_0$)', fontsize=12)
        plt.title('The Stiffness-Flatness Scaling: High Flattening Induces Rigidity', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(output_dir / "stiffness_flatness_scaling.png", dpi=150)
        plt.close()
        log_content.append("\n[Chart Generated]: stiffness_flatness_scaling.png\n")
    except Exception as e:
        log_content.append(f"\n[Chart Error]: {e}\n")

    with open(output_log, "w") as f:
        f.writelines(log_content)
        
    print("\nDiscovery Analysis Complete. See 'discovery_log.txt'.")
    print(f"Log Path: {output_log}")

if __name__ == "__main__":
    main()

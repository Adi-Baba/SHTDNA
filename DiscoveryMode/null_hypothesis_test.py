
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr

def run_null_test():
    cloud_file = r"d:\OnlyHST\SHTDNA\DiscoveryMode\mechanical_state_cloud.csv"
    output_dir = Path(r"d:\OnlyHST\SHTDNA\DiscoveryMode")
    output_log = output_dir / "null_test_log.txt"
    
    df = pd.read_csv(cloud_file)
    features = ['v0_Curvature', 'Eta_r_Shambhian', 'Omega_Mamton']
    data = df[features].dropna().values
    
    n_samples = len(data)
    
    # 1. Real Analysis
    # Correlation
    real_corr, _ = spearmanr(data[:,0], data[:,1]) # v0 vs Eta
    
    # PCA Variance
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(data)
    pca = PCA(n_components=3)
    pca.fit(X_scaled)
    real_var_2d = np.sum(pca.explained_variance_ratio_[:2])
    
    # 2. Null Analysis (Bootstrap Shuffling)
    n_shuffles = 1000
    null_corrs = []
    null_vars = []
    
    print(f"Running {n_shuffles} Null Shuffles...")
    
    for i in range(n_shuffles):
        # Shuffle each column independently to destroy relationships
        null_data = data.copy()
        np.random.shuffle(null_data[:, 0]) # Shuffle v0
        np.random.shuffle(null_data[:, 1]) # Shuffle Eta
        np.random.shuffle(null_data[:, 2]) # Shuffle Omega
        
        # Null Correlation
        c, _ = spearmanr(null_data[:,0], null_data[:,1])
        null_corrs.append(c)
        
        # Null PCA
        X_null_scaled = scaler.fit_transform(null_data) # Re-scale each shuffle
        pca_null = PCA(n_components=3)
        pca_null.fit(X_null_scaled)
        v = np.sum(pca_null.explained_variance_ratio_[:2])
        null_vars.append(v)
        
    # 3. Statistics (Z-Score)
    null_corrs = np.array(null_corrs)
    null_vars = np.array(null_vars)
    
    mean_null_corr = np.mean(null_corrs)
    std_null_corr = np.std(null_corrs)
    z_score_corr = (real_corr - mean_null_corr) / std_null_corr
    
    mean_null_var = np.mean(null_vars)
    std_null_var = np.std(null_vars)
    z_score_var = (real_var_2d - mean_null_var) / std_null_var
    
    # 4. Report
    print("-" * 40)
    print("NULL HYPOTHESIS TEST RESULTS")
    print("-" * 40)
    print("TEST 1: Stiffness-Flatness Correlation")
    print(f"  Real Signal: {real_corr:.4f}")
    print(f"  Null Noise:  {mean_null_corr:.4f} +/- {std_null_corr:.4f}")
    print(f"  Z-Score:     {z_score_corr:.2f} (Sigma)")
    print("-" * 40)
    print("TEST 2: Dimensionality Reduction (PCA)")
    print(f"  Real Variance (2D): {real_var_2d:.4f}")
    print(f"  Null Variance (2D): {mean_null_var:.4f} +/- {std_null_var:.4f}")
    print(f"  Z-Score:            {z_score_var:.2f} (Sigma)")
    print("-" * 40)
    
    # 5. Visualization
    plt.figure(figsize=(10, 5))
    
    # Plot Correlation Dist
    plt.subplot(1, 2, 1)
    plt.hist(null_corrs, bins=30, color='gray', alpha=0.7, label='Null Noise')
    plt.axvline(real_corr, color='red', linewidth=3, label='REAL DATA')
    plt.title(f'Correlation Signal\nZ = {z_score_corr:.1f} $\sigma$')
    plt.xlabel('Spearman Correlation')
    plt.legend()
    
    # Plot PCA Dist
    plt.subplot(1, 2, 2)
    plt.hist(null_vars, bins=30, color='gray', alpha=0.7, label='Null Noise')
    plt.axvline(real_var_2d, color='blue', linewidth=3, label='REAL DATA')
    plt.title(f'2D Manifold Strength\nZ = {z_score_var:.1f} $\sigma$')
    plt.xlabel('Explained Variance (PC1+PC2)')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / "null_test_proof.png", dpi=150)
    
    with open(output_log, "w") as f:
        f.write(f"Z_SCORE_CORR: {z_score_corr:.2f}\n")
        f.write(f"Z_SCORE_VAR: {z_score_var:.2f}\n")
        if abs(z_score_corr) > 5 and abs(z_score_var) > 5:
            f.write("CONCLUSION: PASSED. Impossible to be noise.\n")
        else:
            f.write("CONCLUSION: FAILED. Indistinguishable from noise.\n")

if __name__ == "__main__":
    run_null_test()

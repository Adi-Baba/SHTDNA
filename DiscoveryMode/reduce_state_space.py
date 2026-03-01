
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def analyze_reduction(cloud_file, output_dir):
    df = pd.read_csv(cloud_file)
    output_dir = Path(output_dir)
    
    # 1. Define the State Vector
    # We want to see if [v0, Eta, Omega] collapses.
    # Note: Omega is large (thousands), Eta is 0-1. Scaling is crucial.
    
    features = ['v0_Curvature', 'Eta_r_Shambhian', 'Omega_Mamton']
    data = df[features].dropna()
    X = data.values
    
    # 2. Standardize Features (Mean=0, Std=1)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # 3. Principal Component Analysis (PCA)
    pca = PCA(n_components=3)
    X_pca = pca.fit_transform(X_scaled)
    
    # 4. Explained Variance
    var_ratio = pca.explained_variance_ratio_
    cum_var = np.cumsum(var_ratio)
    
    print("PCA Analysis Result:")
    print(f"PC1 Variance: {var_ratio[0]:.4f}")
    print(f"PC2 Variance: {var_ratio[1]:.4f}")
    print(f"PC3 Variance: {var_ratio[2]:.4f}")
    print(f"Cumulative (PC1+PC2): {cum_var[1]:.4f}")
    
    # 5. Interpretation
    # Loadings (How much each original variable contributes to PCs)
    loadings = pd.DataFrame(
        pca.components_.T, 
        columns=['PC1', 'PC2', 'PC3'], 
        index=features
    )
    print("\nLoadings (Eigenvectors):")
    print(loadings)
    
    # 6. Visualization 1: The Effective 2D Manifold (PC1 vs PC2)
    plt.figure(figsize=(10, 6))
    sc = plt.scatter(X_pca[:, 0], X_pca[:, 1], c=data['v0_Curvature'], cmap='viridis', s=5, alpha=0.5)
    plt.colorbar(sc, label='Curvature $v_0$ (Stiffness)')
    plt.xlabel(f'Principal Component 1 ({var_ratio[0]*100:.1f}%)')
    plt.ylabel(f'Principal Component 2 ({var_ratio[1]*100:.1f}%)')
    plt.title('Low-Dimensional Manifold of DNA Mechanics', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.savefig(output_dir / "pca_manifold.png", dpi=150)
    
    # 7. Visualization 2: The "Phase Diagram" (Eta vs Omega)
    # Showing that these two parameters sort the mechanics
    plt.figure(figsize=(10, 6))
    sc = plt.scatter(data['Eta_r_Shambhian'], data['Omega_Mamton'], c=data['v0_Curvature'], cmap='viridis', s=5, alpha=0.5)
    cbar = plt.colorbar(sc, label='Curvature $v_0$ (Stiffness)')
    
    # Annotate Regions
    plt.text(0.95, 50, "Stiff\nRibbon", color='white', ha='center', fontweight='bold', bbox=dict(facecolor='black', alpha=0.5))
    plt.text(0.75, 4000, "Flexible\nCoil", color='white', ha='center', fontweight='bold', bbox=dict(facecolor='black', alpha=0.5))
    
    plt.xlabel('Shambhian Twist $\eta_r$ (Flatness)', fontsize=12)
    plt.ylabel('Mamton $\Omega$ (Roughness)', fontsize=12)
    plt.title('Geometric State Map: Separation of Stiff/Flexible Modes', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.savefig(output_dir / "geometric_state_map.png", dpi=150)
    
    # Save Log
    with open(output_dir / "reduction_log.txt", "w") as f:
        f.write("# Dimensionality Reduction Analysis\n\n")
        f.write(f"Explained Variance (PC1+PC2): {cum_var[1]*100:.2f}%\n")
        f.write("\n## Loadings\n")
        f.write(loadings.to_string())
        f.write("\n\n## Conclusion\n")
        if cum_var[1] > 0.8:
            f.write("System is effectively 2D. (Success)\n")
        else:
            f.write("System requires 3 dimensions.\n")

if __name__ == "__main__":
    cloud_file = r"d:\OnlyHST\SHTDNA\DiscoveryMode\mechanical_state_cloud.csv"
    extract_dir = r"d:\OnlyHST\SHTDNA\DiscoveryMode"
    
    print("Running SHT Reduction Analysis...")
    analyze_reduction(cloud_file, extract_dir)
    print(f"Done. See {extract_dir}")

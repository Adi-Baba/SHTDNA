#!/usr/bin/env python3
"""
Coupling Analysis: Discovery Path 2
Analyze causal relationships and coupling between Vij, Shambhian, and Mamton units.
"""

import numpy as np
import pandas as pd
from pathlib import Path
import json
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score


class CouplingAnalyzer:
    def __init__(self, base_dir='.'):
        self.base_dir = Path(base_dir)
        self.results = {}
        
    def load_data(self):
        """Load and merge results from all three modes."""
        curvature_file = self.base_dir / 'CurvatureMode' / 'curvature_results.csv'
        eccentricity_file = self.base_dir / 'EccentricityMode' / 'eccentricity_results.csv'
        roughness_file = self.base_dir / 'RoughnessMode' / 'roughness_results.csv'
        
        # Load CSVs
        df_curv = pd.read_csv(curvature_file)
        df_ecc = pd.read_csv(eccentricity_file)
        df_rough = pd.read_csv(roughness_file)
        
        # Rename columns for clarity
        df_curv = df_curv.rename(columns={
            'mean_kappa': 'vij_mean',
            'std_kappa': 'vij_std',
            'max_kappa': 'vij_max'
        })
        df_ecc = df_ecc.rename(columns={
            'mean_epsilon': 'shambhian_mean',
            'std_epsilon': 'shambhian_std',
            'max_epsilon': 'shambhian_max'
        })
        df_rough = df_rough.rename(columns={
            'mean_omega': 'mamton_mean',
            'std_omega': 'mamton_std',
            'max_omega': 'mamton_max'
        })
        
        # Merge on pdb_id
        df_merged = df_curv.merge(df_ecc, on='pdb_id', suffixes=('_curv', '_ecc'))
        df_merged = df_merged.merge(df_rough, on='pdb_id', suffixes=('', '_rough'))
        
        # Keep only one len_h column
        if 'len_h_curv' in df_merged.columns:
            df_merged['len_h'] = df_merged['len_h_curv']
            df_merged = df_merged.drop(columns=[c for c in df_merged.columns if 'len_h' in c and c != 'len_h'])
        
        print(f"Merged data: {len(df_merged)} structures with all three measurements")
        self.df = df_merged
        return df_merged
    
    def compute_correlations(self):
        """Compute pairwise correlations between units."""
        # Define key variables
        vars_dict = {
            'Vij (Curvature)': 'vij_mean',
            'Shambhian (Flatness)': 'shambhian_mean',
            'Mamton (Roughness)': 'mamton_mean'
        }
        
        results = {}
        print("\n=== PAIRWISE CORRELATIONS ===\n")
        
        for (name1, var1), (name2, var2) in [
            (('Vij (Curvature)', 'vij_mean'), ('Shambhian (Flatness)', 'shambhian_mean')),
            (('Vij (Curvature)', 'vij_mean'), ('Mamton (Roughness)', 'mamton_mean')),
            (('Shambhian (Flatness)', 'shambhian_mean'), ('Mamton (Roughness)', 'mamton_mean'))
        ]:
            # Filter valid data
            valid = self.df[[var1, var2]].dropna()
            
            # Pearson
            r_pearson, p_pearson = pearsonr(valid[var1], valid[var2])
            
            # Spearman
            r_spearman, p_spearman = spearmanr(valid[var1], valid[var2])
            
            # Store results
            pair_key = f"{name1} vs {name2}"
            results[pair_key] = {
                'n': len(valid),
                'pearson_r': float(r_pearson),
                'pearson_p': float(p_pearson),
                'spearman_r': float(r_spearman),
                'spearman_p': float(p_spearman)
            }
            
            print(f"{pair_key}:")
            print(f"  N = {len(valid)}")
            print(f"  Pearson:  r = {r_pearson:+.4f}, p = {p_pearson:.2e}")
            print(f"  Spearman: ρ = {r_spearman:+.4f}, p = {p_spearman:.2e}")
            print()
        
        self.results['correlations'] = results
        return results
    
    def test_causal_arrows(self):
        """Test directional relationships using predictive power."""
        print("\n=== CAUSAL ARROW ANALYSIS ===\n")
        print("Testing predictive power: X → Y (can X predict Y?)\n")
        
        arrows = []
        
        # Test all directional pairs
        pairs = [
            ('vij_mean', 'shambhian_mean', 'Vij', 'Shambhian'),
            ('shambhian_mean', 'vij_mean', 'Shambhian', 'Vij'),
            ('vij_mean', 'mamton_mean', 'Vij', 'Mamton'),
            ('mamton_mean', 'vij_mean', 'Mamton', 'Vij'),
            ('shambhian_mean', 'mamton_mean', 'Shambhian', 'Mamton'),
            ('mamton_mean', 'shambhian_mean', 'Mamton', 'Shambhian'),
        ]
        
        for x_var, y_var, x_name, y_name in pairs:
            # Prepare data
            valid = self.df[[x_var, y_var]].dropna()
            X = valid[[x_var]].values
            y = valid[y_var].values
            
            # Random Forest predictor
            rf = RandomForestRegressor(n_estimators=100, max_depth=5, random_state=42)
            scores = cross_val_score(rf, X, y, cv=5, scoring='r2')
            r2_mean = scores.mean()
            r2_std = scores.std()
            
            # Fit to get feature importance
            rf.fit(X, y)
            
            arrow_result = {
                'from': x_name,
                'to': y_name,
                'r2_mean': float(r2_mean),
                'r2_std': float(r2_std),
                'predictive': bool(r2_mean > 0.1)
            }
            arrows.append(arrow_result)
            
            direction = "→" if r2_mean > 0.1 else "⇢"
            print(f"{x_name} {direction} {y_name}:  R² = {r2_mean:.3f} ± {r2_std:.3f}")
        
        self.results['causal_arrows'] = arrows
        return arrows
    
    def identify_outliers(self, threshold_sigma=3):
        """Identify structures with extreme geometric states."""
        print(f"\n=== OUTLIER DETECTION (>{threshold_sigma}σ) ===\n")
        
        outliers = []
        
        for var, name in [
            ('vij_mean', 'Vij'),
            ('shambhian_mean', 'Shambhian'),
            ('mamton_mean', 'Mamton')
        ]:
            mean = self.df[var].mean()
            std = self.df[var].std()
            
            # Find structures beyond threshold
            z_scores = (self.df[var] - mean) / std
            extreme = self.df[np.abs(z_scores) > threshold_sigma].copy()
            extreme['z_score'] = z_scores[np.abs(z_scores) > threshold_sigma]
            extreme['variable'] = name
            
            if len(extreme) > 0:
                print(f"{name} outliers ({len(extreme)} structures):")
                for _, row in extreme.nlargest(5, 'z_score').iterrows():
                    print(f"  {row['pdb_id']}: {name} = {row[var]:.4f} (z = {row['z_score']:+.2f}σ)")
                print()
                
                outliers.append(extreme[['pdb_id', 'variable', var, 'z_score']])
        
        if outliers:
            self.outliers_df = pd.concat(outliers, ignore_index=True)
            self.results['n_outliers'] = len(self.outliers_df)
        else:
            self.outliers_df = pd.DataFrame()
            self.results['n_outliers'] = 0
        
        return self.outliers_df
    
    def geometric_clustering(self):
        """Find natural clusters in 3D geometric space."""
        print("\n=== GEOMETRIC STATE CLUSTERING ===\n")
        
        # Normalize variables for clustering
        from sklearn.preprocessing import StandardScaler
        from sklearn.cluster import KMeans
        
        X = self.df[['vij_mean', 'shambhian_mean', 'mamton_mean']].dropna()
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        
        # Find optimal number of clusters using elbow method
        inertias = []
        K_range = range(2, 8)
        for k in K_range:
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
            kmeans.fit(X_scaled)
            inertias.append(kmeans.inertia_)
        
        # Use 3-4 clusters for interpretability
        n_clusters = 4
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        self.df['cluster'] = np.nan
        self.df.loc[X.index, 'cluster'] = kmeans.fit_predict(X_scaled)
        
        # Characterize clusters
        print(f"Found {n_clusters} geometric states:\n")
        for i in range(n_clusters):
            cluster_data = self.df[self.df['cluster'] == i]
            print(f"Cluster {i} (N={len(cluster_data)}):")
            print(f"  Mean Vij:       {cluster_data['vij_mean'].mean():.4f}")
            print(f"  Mean Shambhian: {cluster_data['shambhian_mean'].mean():.4f}")
            print(f"  Mean Mamton:    {cluster_data['mamton_mean'].mean():.1f}")
            print(f"  Example PDBs: {', '.join(cluster_data['pdb_id'].head(3).tolist())}")
            print()
        
        self.results['n_clusters'] = n_clusters
        return self.df
    
    def generate_plots(self, output_dir):
        """Generate visualization plots."""
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        
        # 1. Correlation matrix heatmap
        fig, ax = plt.subplots(figsize=(8, 6))
        corr_vars = ['vij_mean', 'shambhian_mean', 'mamton_mean']
        corr_matrix = self.df[corr_vars].corr(method='spearman')
        sns.heatmap(corr_matrix, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
                    xticklabels=['Vij', 'Shambhian', 'Mamton'],
                    yticklabels=['Vij', 'Shambhian', 'Mamton'],
                    vmin=-1, vmax=1, ax=ax)
        plt.title('Spearman Correlation Matrix\n(Geometric Units)', fontsize=14)
        plt.tight_layout()
        plt.savefig(output_dir / 'correlation_matrix.png', dpi=150)
        plt.close()
        
        # 2. Pairwise scatter plots
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        
        pairs = [
            ('vij_mean', 'shambhian_mean', 'Vij (Curvature)', 'Shambhian (Flatness)'),
            ('vij_mean', 'mamton_mean', 'Vij (Curvature)', 'Mamton (Roughness)'),
            ('shambhian_mean', 'mamton_mean', 'Shambhian (Flatness)', 'Mamton (Roughness)')
        ]
        
        for ax, (var1, var2, name1, name2) in zip(axes, pairs):
            valid = self.df[[var1, var2]].dropna()
            r, p = spearmanr(valid[var1], valid[var2])
            
            ax.scatter(valid[var1], valid[var2], alpha=0.5, s=20)
            ax.set_xlabel(name1)
            ax.set_ylabel(name2)
            ax.set_title(f'ρ = {r:+.3f}, p = {p:.2e}')
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'pairwise_scatter.png', dpi=150)
        plt.close()
        
        # 3. 3D geometric state space
        if 'cluster' in self.df.columns:
            from mpl_toolkits.mplot3d import Axes3D
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection='3d')
            
            for cluster_id in self.df['cluster'].dropna().unique():
                cluster_data = self.df[self.df['cluster'] == cluster_id]
                ax.scatter(cluster_data['vij_mean'],
                          cluster_data['shambhian_mean'],
                          np.log10(cluster_data['mamton_mean'] + 1),
                          label=f'State {int(cluster_id)}',
                          alpha=0.6, s=30)
            
            ax.set_xlabel('Vij (Curvature)')
            ax.set_ylabel('Shambhian (Flatness)')
            ax.set_zlabel('log₁₀(Mamton)')
            ax.set_title('3D Geometric State Space')
            ax.legend()
            plt.tight_layout()
            plt.savefig(output_dir / 'geometric_state_space.png', dpi=150)
            plt.close()
        
        print(f"\nPlots saved to {output_dir}")
    
    def save_results(self, output_dir):
        """Save analysis results."""
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        
        # Save JSON summary
        with open(output_dir / 'coupling_analysis_summary.json', 'w') as f:
            json.dump(self.results, f, indent=2)
        
        # Save merged data
        self.df.to_csv(output_dir / 'merged_geometric_data.csv', index=False)
        
        # Save outliers
        if not self.outliers_df.empty:
            self.outliers_df.to_csv(output_dir / 'geometric_outliers.csv', index=False)
        
        print(f"\nResults saved to {output_dir}")


def main():
    print("=" * 60)
    print("COUPLING ANALYSIS: Discovery Path 2")
    print("Analyzing relationships between Vij, Shambhian, and Mamton")
    print("=" * 60)
    
    # Initialize analyzer
    analyzer = CouplingAnalyzer()
    
    # Run analysis pipeline
    df = analyzer.load_data()
    
    correlations = analyzer.compute_correlations()
    
    arrows = analyzer.test_causal_arrows()
    
    outliers = analyzer.identify_outliers(threshold_sigma=2.5)
    
    clustered_df = analyzer.geometric_clustering()
    
    # Generate visualizations
    output_dir = Path('DiscoveryMode/coupling_results')
    analyzer.generate_plots(output_dir)
    
    # Save results
    analyzer.save_results(output_dir)
    
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"\nKey Findings:")
    print(f"  • {len(df)} structures analyzed")
    print(f"  • {len(outliers)} geometric outliers identified")
    if 'n_clusters' in analyzer.results:
        print(f"  • {analyzer.results['n_clusters']} distinct geometric states")
    print(f"\nResults: {output_dir}")


if __name__ == '__main__':
    main()

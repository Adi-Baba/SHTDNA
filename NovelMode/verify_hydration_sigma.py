
import pandas as pd
import numpy as np
from scipy.stats import norm

# Fisher Z-Transformation
# Z = 0.5 * ln((1+r)/(1-r))
# Sigma = Z * sqrt(N-3)

def compute_sigma():
    csv_path = r"d:\OnlyHST\SHTDNA\NovelMode\hydration_results.csv"
    
    # Load Data (No headers in preview, assuming columns)
    # Based on output: pdb, rough_mean, rough_corr, density_corr?, n_waters?, n_slices?
    # Wait, the preview had: `4IX7,48.69...,0.74...,-0.13...,20,26,0.76...`
    # Let's assume columns: PDB, MeanOmega, OmegaCorr?, WaterCorr(Target), NWaters, NSlices, Ratio
    # The report says 4RTK has R = -0.859.
    
    # Let's search for 4RTK in the file
    try:
        df = pd.read_csv(csv_path, header=None)
        # Rename for clarity
        df.columns = ['pdb', 'mean_omega', 'omega_autocorr', 'water_corr', 'n_dense_slices', 'n_total_slices', 'ratio']
        
        # Convert to numeric, errors='coerce' to handle bad strings
        for col in ['water_corr', 'n_total_slices']:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # 1. Analyze 4RTK
        row = df[df['pdb'] == '4RTK']
        if not row.empty:
            r = float(row['water_corr'].values[0])
            n = float(row['n_total_slices'].values[0])
            
            # Fisher Z
            # r is the correlation between Omega and Water Density
            # The report said R = -0.859
            
            if abs(r) > 0.999: # Safety
                fz = 3.0
            else:
                fz = 0.5 * np.log((1 + abs(r)) / (1 - abs(r)))
            
            sigma = fz * np.sqrt(n - 3)
            
            # P-value (2-tailed)
            p_val = 2 * (1 - norm.cdf(abs(sigma)))
            
            print(f"Structure: 4RTK")
            print(f"Correlation (R): {r:.3f}")
            print(f"Sample Size (N): {n} slices")
            print(f"Fisher Z-Score: {float(sigma):.2f} Sigma")
            print(f"P-Value: {p_val:.2e}")
            
            # 2. Analyze 'Top 8' Cluster
            # How unlikely is it to find K structures with R < -0.6 in a dataset of M structures by chance?
            # Assuming Null Distribution of R is centered at 0 with width 1/sqrt(N-3).
            # Let's count how many have Z > 3.0 sigma intrinsically.
            
            df_valid = df[df['water_corr'] != -1.0] # Valid ones
            
            # Calculate Z for all
            def get_z(r, n):
                if n <= 3: return 0
                if abs(r) >= 1: return 0
                return (0.5 * np.log((1+abs(r))/(1-abs(r)))) * np.sqrt(n-3)
                
            df_valid['z_score'] = df_valid.apply(lambda x: get_z(x['water_corr'], x['n_total_slices']), axis=1)
            
            top_z = df_valid[df_valid['water_corr'] < -0.5]
            print("\nTop Structures (R < -0.5):")
            print(top_z[['pdb', 'water_corr', 'n_total_slices', 'z_score']])
            
            # Meta-Statistic
            # Stouffer's Z-score method for the cluster?
            # Or just "Max Z".
            max_z = df_valid['z_score'].max()
            print(f"\nMax Observed Significance: {max_z:.2f} Sigma")
            
        else:
            print("4RTK not found in CSV.")
            
    except Exception as e:
        print(e)
        
if __name__ == "__main__":
    compute_sigma()

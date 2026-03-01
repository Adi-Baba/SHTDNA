
import pandas as pd
import numpy as np
import sys

# Load data
csv_path = r"d:\OnlyHST\SHTDNA\PDB\sht_analysis_output\sht_analysis_results.csv"
output_path = r"d:\OnlyHST\SHTDNA\PDB\dataset\analysis_summary.txt"

with open(output_path, 'w') as f:
    orig_stdout = sys.stdout
    sys.stdout = f
    
    try:
        df = pd.read_csv(csv_path)
        print(f"Loaded {len(df)} records.")

        # Clean data 
        def clean_res(x):
            try:
                return float(x.strip('[]'))
            except:
                return np.nan

        df['resolution_val'] = df['resolution'].apply(clean_res)

        # 1. verify Lk = Tw + Wr relation implicit in sigma
        df['Lk_SHT'] = df['Tw_SHT'] + df['Wr_SHT']
        df['sigma_recalc'] = (df['Lk_SHT'] - df['Tw0']) / df['Tw0']
        df['sigma_diff'] = df['sigma_SHT_global'] - df['sigma_recalc']

        print("\n--- Consistency Check ---")
        print(f"Mean difference between reported and recalced sigma: {df['sigma_diff'].mean():.6e}")
        print(f"Max difference: {df['sigma_diff'].abs().max():.6e}")

        # 2. Statistics of Sigma
        print("\n--- Sigma SHT Statistics ---")
        print(df['sigma_SHT_global'].describe())

        # 3. Method comparison
        print("\n--- By Method ---")
        print(df.groupby('method')['sigma_SHT_global'].agg(['count', 'mean', 'std']))

        # 4. Outliers
        threshold = 0.5
        outliers = df[df['sigma_SHT_global'].abs() > threshold]
        print(f"\n--- Outliers (|sigma| > {threshold}) ---")
        print(f"Count: {len(outliers)}")
        if len(outliers) > 0:
            print(outliers[['pdb_id', 'method', 'sigma_SHT_global', 'Tw_SHT', 'Wr_SHT']].head(10))

        # 5. Tw vs Wr
        print("\n--- Twist vs Writhe ---")
        print(f"Mean Tw: {df['Tw_SHT'].mean():.4f}, Mean Wr: {df['Wr_SHT'].mean():.4f}")
        print(f"Correlation: {df[['Tw_SHT', 'Wr_SHT']].corr().iloc[0,1]:.4f}")

        # 6. Resolution analysis
        print("\n--- Resolution Correlation ---")
        print(f"Correlation between Resolution and Sigma error (diff): {df[['resolution_val', 'sigma_diff']].corr().iloc[0,1]:.4f}")
        
    except Exception as e:
        print(f"Error: {e}")
    finally:
        sys.stdout = orig_stdout

print("Done. Output written to analysis_summary.txt")

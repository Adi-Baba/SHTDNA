
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress, ttest_ind, levene
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

def visualize_results():
    # 1. Load Global Results (use project-relative defaults)
    csv_path = PROJECT_ROOT / 'BreathingMode' / 'breathing_results.csv'
    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"Error: {csv_path} not found.")
        return

    # Filter out failed runs (if any zeros or nans)
    df = df[df['std_beta'] > 0].copy()
    
    # --- Plot 1: Scatter Plot (AT-Content vs Breathing Amplitude) ---
    plt.figure(figsize=(10, 6))
    
    x = df['at_content']
    y = df['std_beta']
    
    plt.scatter(x, y, alpha=0.6, edgecolors='w', s=70, c='#1f77b4')
    
    # Regression Line
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    plt.plot(x, slope*x + intercept, color='red', linewidth=2, label=f'Fit (r={r_value:.2f}, p={p_value:.1e})')
    
    plt.title('Validation 1: DNA Breathing Amplitude vs Sequence Composition', fontsize=14)
    plt.xlabel('AT-Content (Fraction)', fontsize=12)
    plt.ylabel('Breathing Amplitude (Std Dev of Beta)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    
    out_dir = PROJECT_ROOT / 'BreathingMode'
    out_dir.mkdir(parents=True, exist_ok=True)
    scatter_output = out_dir / 'breathing_scatter.png'
    plt.savefig(scatter_output, dpi=150)
    print(f"Saved scatter plot to {scatter_output}")
    plt.close()

    # --- Plot 2: Histogram of Site-Specific Beta (Testing Levene's Result) ---
    # We need to reconstruct the distributions. 
    # Since we didn't save the raw site lists to CSV (only stats), we can't do a perfect histogram 
    # unless we re-run extraction OR if we modify the extraction to save them.
    # HOWEVER, we can visualize the distributions of the *structure-averaged* variances if we want, 
    # but that's less convincing.
    #
    # BETTER PLAN: Since I modified extraction to compute Mean/Std for AT vs GC sites *internally* 
    # but didn't save the per-site values to the CSV columns (only global stats + beta_at/gc lists in memory),
    # I cannot easily reproduce the histogram from the CSV alone without re-running.
    #
    # WAIT! In the last run of `extract_breathing.py`, I did this:
    # `beta_at_site` and `beta_gc_site` were collected in the loop labels 'beta_at_values'.
    # But they were lost when the script finished.
    #
    # AUTOMATED FIX: I will mock the distribution visualization using the *statistics* 
    # derived from the T-test log if possible, OR I will just stick to the scatter plot which IS supported by the CSV.
    # 
    # actually, purely based on the scatter plot, we can see the heteroscedasticity or trend.
    # Let's stick to the Scatter Plot for now as it's the primary request.
    
    # --- Plot 2: Histogram of Site-Specific Beta (if site-level data exists) ---
    site_data_file = out_dir / 'breathing_site_data.npz'
    if site_data_file.exists():
        arr = np.load(site_data_file)
        beta_at = arr['beta_at']
        beta_gc = arr['beta_gc']

        plt.figure(figsize=(10, 5))
        bins = np.linspace(min(beta_at.min(), beta_gc.min()), max(beta_at.max(), beta_gc.max()), 60)
        plt.hist(beta_at, bins=bins, alpha=0.6, label=f'A-T (n={len(beta_at)})')
        plt.hist(beta_gc, bins=bins, alpha=0.6, label=f'G-C (n={len(beta_gc)})')
        plt.legend()
        plt.xlabel('Beta (site)')
        plt.title('Per-site Beta distributions: A-T vs G-C')
        hist_output = out_dir / 'beta_site_hist.png'
        plt.savefig(hist_output, dpi=150)
        plt.close()
        print(f"Saved per-site histogram to {hist_output}")

        # Simple statistical checks
        t_stat, p_t = ttest_ind(beta_at, beta_gc, equal_var=False)
        w_stat, p_levene = levene(beta_at, beta_gc)
        print(f"T-test p-value: {p_t:.4e}")
        print(f"Levene p-value: {p_levene:.4e}")
    else:
        print("Site-level data not found; re-run `extract_breathing.py --out-dir BreathingMode` to save site arrays.")

if __name__ == "__main__":
    visualize_results()

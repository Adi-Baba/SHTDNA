
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit

def rayleigh_pdf(x, sigma, scale):
    return scale * (x / sigma**2) * np.exp(-x**2 / (2 * sigma**2))

def compute_Lp_for_subset(v0_data):
    if len(v0_data) < 50: return None
    
    # Histogram
    bins = 30
    counts, bin_edges = np.histogram(v0_data, bins=bins, density=True)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    valid = counts > 0
    x = bin_centers[valid]
    P = counts[valid]
    
    try:
        popt, pcov = curve_fit(rayleigh_pdf, x, P, p0=[0.1, 1.0], maxfev=5000)
        sigma_fit = popt[0]
        # Lp = 1 / (sigma^2 * b)
        b_rise = 3.4
        Lp_nm = (1.0 / (sigma_fit**2 * b_rise)) / 10.0
        return Lp_nm
    except:
        return None

def main():
    cloud_file = r"d:\OnlyHST\SHTDNA\DiscoveryMode\mechanical_state_cloud.csv"
    output_dir = Path(r"d:\OnlyHST\SHTDNA\ThermodynamicsMode")
    output_img = output_dir / "stiffness_convergence.png"
    output_log = output_dir / "lost_stiffness_log.txt"
    
    df = pd.read_csv(cloud_file)
    
    # Thresholds to scan: Eta_r from 0.7 to 0.99
    thresholds = np.linspace(0.70, 0.98, 20)
    results_Lp = []
    results_Thresh = []
    sample_sizes = []
    
    print("Scanning Stiffness vs Flatness Threshold...")
    
    for th in thresholds:
        # Select slices STRONGER than this flatness
        subset = df[df['Eta_r_Shambhian'] > th]
        v0 = subset['v0_Curvature'].dropna().values
        
        Lp = compute_Lp_for_subset(v0)
        
        if Lp is not None:
            print(f"Eta_r > {th:.3f} : N={len(v0)} : Lp={Lp:.2f} nm")
            results_Lp.append(Lp)
            results_Thresh.append(th)
            sample_sizes.append(len(v0))
        else:
            print(f"Eta_r > {th:.3f} : Insufficient Data")
            
    # Plotting
    fig, ax1 = plt.subplots(figsize=(10, 6))

    color = 'tab:red'
    ax1.set_xlabel('Shambhian Twist Threshold ($\eta_r$)', fontsize=12)
    ax1.set_ylabel('Predicted Persistence Length ($L_p$) [nm]', color=color, fontsize=12)
    ax1.plot(results_Thresh, results_Lp, 'o-', color=color, linewidth=2, label='Stiffness ($L_p$)')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(True, alpha=0.3)
    
    # Reference Line 50nm
    ax1.axhline(y=50.0, color='green', linestyle='--', label='Theoretical Ideal (50nm)')
    
    ax2 = ax1.twinx()  
    color = 'tab:blue'
    ax2.set_ylabel('Sample Size (slices)', color=color, fontsize=12)  
    ax2.plot(results_Thresh, sample_sizes, 'x:', color=color, label='Sample Size')
    ax2.tick_params(axis='y', labelcolor=color)
    
    plt.title('Finding the "Lost" Stiffness: Convergence to 50nm', fontsize=14)
    fig.tight_layout()  
    plt.savefig(output_img, dpi=150)
    
    # Save Log
    with open(output_log, 'w') as f:
        f.write("# Stiffness Convergence Analysis\n")
        f.write("Does Lp approach 50nm as we select perfect crystals?\n\n")
        f.write("| Eta_r Threshold | Sample Size | Predicted Lp (nm) |\n")
        f.write("|---|---|---|\n")
        for t, n, l in zip(results_Thresh, sample_sizes, results_Lp):
            f.write(f"| > {t:.3f} | {n} | **{l:.2f}** |\n")
            
    print(f"Scan Complete. Max Lp found: {max(results_Lp):.2f} nm")
    print(f"See {output_img}")

if __name__ == "__main__":
    main()


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit

def compute_energy_landscape(data_file):
    # 1. Load the "Mechanical Cloud"
    df = pd.read_csv(data_file)
    # Filter for "Flat" DNA (Eta_r > 0.9) to remove intrinsic curvature noise
    df_flat = df[df['Eta_r_Shambhian'] > 0.9]
    v0_data = df_flat['v0_Curvature'].dropna().values
    
    # 2. Compute Probability Density P(v0)
    # Use bins to estimate the histogram
    bins = 50
    counts, bin_edges = np.histogram(v0_data, bins=bins, density=True)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    
    # Filter empty bins
    valid = counts > 0
    x = bin_centers[valid]
    P = counts[valid]
    
    # 4. Rayleigh Fit (Correct Physics for Magnitude of Curvature)
    # Theory: kappa vector k ~ N(0, sigma^2) in 2D.
    # Magnitude x = |k| follows Rayleigh: P(x) = (x / sigma^2) * exp(-x^2 / (2*sigma^2))
    # Relationship to Lp:
    # E = 1/2 * Lp * b * k^2 (in kT)
    # P(k) ~ exp(-1/2 * Lp * b * k^2)
    # Comparing to Gaussian form exp(-k^2 / 2sigma^2):
    # 1/2 * Lp * b = 1 / (2 * sigma^2)
    # Lp = 1 / (sigma^2 * b)
    
    b_rise = 3.4 # Angstroms

    def rayleigh_pdf(x, sigma, scale):
        return scale * (x / sigma**2) * np.exp(-x**2 / (2 * sigma**2))

    try:
        # Fit Probability Density P(x) directly, not V_eff
        # x are bin centers, P are density values
        popt, pcov = curve_fit(rayleigh_pdf, x, P, p0=[0.1, 1.0])
        sigma_fit = popt[0]
        
        # Calculate Persistence Length
        # Lp = 1 / (sigma^2 * b)
        Lp_angstrom = 1.0 / (sigma_fit**2 * b_rise)
        Lp_nm = Lp_angstrom / 10.0
        
        return {
            'x': x,
            'P': P, # Return density
            'Lp_nm': Lp_nm,
            'sigma': sigma_fit,
            'popt': popt,
            'n_samples': len(v0_data)
        }
        
    except Exception as e:
        print(f"Fitting failed: {e}")
        return None

def main():
    cloud_file = r"d:\OnlyHST\SHTDNA\DiscoveryMode\mechanical_state_cloud.csv"
    output_dir = Path(r"d:\OnlyHST\SHTDNA\ThermodynamicsMode")
    output_log = output_dir / "thermo_log.txt"
    output_img = output_dir / "energy_landscape.png"
    
    print("Deriving Energy Landscape for FLAT DNA (Eta_r > 0.9)...")
    res = compute_energy_landscape(cloud_file)
    
    if res is None:
        print("Analysis failed.")
        return

    # --- Results ---
    Lp_nm = res['Lp_nm']
    print(f"Predicted Persistence Length (Flat DNA): {Lp_nm:.4f} nm")
    
    # --- Visualization ---
    x = res['x']
    P = res['P']
    popt = res['popt']
    Lp_nm = res['Lp_nm']
    sigma = res['sigma']
    
    def rayleigh_pdf(x, sigma, scale):
        return scale * (x / sigma**2) * np.exp(-x**2 / (2 * sigma**2))
    
    plt.figure(figsize=(10, 6))
    plt.bar(x, P, width=(x[1]-x[0])*0.9, alpha=0.5, color='blue', label=f'Observed Density (N={res["n_samples"]})')
    
    # Plot Fit
    x_smooth = np.linspace(0, max(x), 100)
    y_smooth = rayleigh_pdf(x_smooth, *popt)
    plt.plot(x_smooth, y_smooth, 'r-', linewidth=2, label=f'Rayleigh Fit ($L_p$={Lp_nm:.1f} nm)')
    
    plt.xlabel('Curvature Magnitude $v_0$ ($Å^{-1}$)', fontsize=12)
    plt.ylabel('Probability Density', fontsize=12)
    plt.title('Persistence Length of Stiff/Flat DNA (Eta_r > 0.9)', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(output_img, dpi=150)
    
    # --- Log ---
    with open(output_log, "w") as f:
        f.write("# SHT-Thermodynamics Result (Subset: Flat DNA > 0.9)\n")
        f.write(f"Predicted Persistence Length (Lp): {Lp_nm:.4f} nm\n")
        f.write(f"Fitted Sigma: {sigma:.4f}\n")
        f.write("\n## Interpretation\n")
        
        if 40 <= Lp_nm <= 60:
            f.write("SUCCESS: Prediction matches experimental physics (~50nm).\n")
            f.write("This confirms 'Option 1': Static structure predicts Dynamics.\n")
        else:
            f.write(f"DEVIATION: Predicted {Lp_nm:.1f}nm vs Expected 50nm.\n")
            f.write("Possible reasons: Crystal packing forces, non-canonical ensemble.\n")

    print(f"Done. See {output_log}")

if __name__ == "__main__":
    main()


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit

# 1. Reuse Physics Model
def geometric_stiffness(eta, L_coil, gamma):
    # Avoid singularity
    return L_coil * (1.0 - np.minimum(eta, 0.999))**(-gamma)

def calculate_penalty(v0_data, eta_data):
    # a. bin data
    bins = 20
    # We need Lp vs Eta for this specific shuffled dataset
    # This is hard because "Lp" is derived from v0 distribution in specific bins
    # Simplified approach:
    # 1. Fit the Scaling Law directly to the (Eta, v0) cloud?
    # No, we need the "Variance of v0" in each Eta bin.
    # Because Lp ~ 1/Variance.
    
    # Divide into bins of Eta
    n_bins = 10
    eta_edges = np.linspace(0.7, 0.98, n_bins+1)
    
    eta_centers = []
    Lp_values = []
    
    for i in range(n_bins):
        mask = (eta_data >= eta_edges[i]) & (eta_data < eta_edges[i+1])
        if np.sum(mask) > 20:
            v0_subset = v0_data[mask]
            # Estimate Lp from variance of v0 assuming Rayleigh
            # Var(v0) = (4 - pi)/2 * sigma^2
            # Lp ~ 1/sigma^2
            # So Lp ~ 1/Var(v0)
            var_v0 = np.var(v0_subset)
            if var_v0 > 1e-9:
                # Absolute calibration doesn't matter for Z-score as much as relative
                # But let's try to be consistent. 
                # Lp = A / <v0^2> approx?
                # Rayleigh: <x^2> = 2 sigma^2.
                # Lp = 1 / (sigma^2 * b) = 1 / (<x^2>/2 * b)
                mean_sq = np.mean(v0_subset**2)
                b = 3.4
                Lp_est = 1.0 / ( (mean_sq/2.0) * b ) * 0.1 # nm
                
                eta_centers.append( (eta_edges[i] + eta_edges[i+1])/2 )
                Lp_values.append(Lp_est)
                
    if len(Lp_values) < 3: return 0.0 # Failed fit
    
    # b. Fit Scaling Law: Lp = A * (1-eta)^-g
    try:
        def model(x, A, g): return A * (1-x)**(-g)
        popt, _ = curve_fit(model, eta_centers, Lp_values, p0=[0.1, 1], maxfev=1000)
        A, g = popt
        
        # c. Predict for Generic (0.90) and Flat (0.96)
        Lp_gen = model(0.90, A, g)
        Lp_flat = model(0.96, A, g)
        
        # d. Calculate Energy
        # E = 0.5 * Lp * (L/R^2). Factor = 0.5 * 50 / 4.3^2 = 1.35
        factor = 1.35
        E_gen = Lp_gen * factor
        E_flat = Lp_flat * factor
        
        return E_flat - E_gen # The Penalty
        
    except:
        return 0.0

def run_sigma_test():
    cloud_file = r"d:\OnlyHST\SHTDNA\DiscoveryMode\mechanical_state_cloud.csv"
    output_dir = Path(r"d:\OnlyHST\SHTDNA\PredictionMode")
    output_img = output_dir / "prediction_significance.png"
    
    print("Loading Cloud...")
    df = pd.read_csv(cloud_file)
    v0_real = df['v0_Curvature'].values
    eta_real = df['Eta_r_Shambhian'].values
    valid = ~np.isnan(v0_real) & ~np.isnan(eta_real)
    v0_real = v0_real[valid]
    eta_real = eta_real[valid]
    
    # 1. Real Penalty
    real_penalty = calculate_penalty(v0_real, eta_real)
    print(f"Real Prediction Penalty: {real_penalty:.2f} kBT")
    
    # 2. Null Shuffles
    n_shuffles = 500
    null_penalties = []
    
    print(f"Running {n_shuffles} Shuffles...")
    for k in range(n_shuffles):
        # Shuffle v0 against Eta (break link)
        v0_shuffled = v0_real.copy()
        np.random.shuffle(v0_shuffled)
        
        p = calculate_penalty(v0_shuffled, eta_real)
        null_penalties.append(p)
        
    null_penalties = np.array(null_penalties)
    
    # 3. Statistics
    mean_null = np.mean(null_penalties)
    std_null = np.std(null_penalties)
    z_score = (real_penalty - mean_null) / (std_null + 1e-9)
    
    print("-" * 30)
    print(f"Real Penalty: {real_penalty:.2f} kBT")
    print(f"Null Baseline: {mean_null:.2f} +/- {std_null:.2f} kBT")
    print(f"Z-SCORE (SIGMA): {z_score:.2f}")
    print("-" * 30)
    
    # 4. Visualization
    plt.figure(figsize=(10, 6))
    plt.hist(null_penalties, bins=30, color='gray', alpha=0.7, label='Random Noise Model')
    plt.axvline(real_penalty, color='red', linewidth=3, label=f'REAL PREDICTION ({real_penalty:.1f} kBT)')
    plt.title(f'Statistical Significance of Nucleosome Exclusion Prediction\nZ = {z_score:.1f} $\sigma$', fontsize=14)
    plt.xlabel('Predicted Energy Penalty ($k_B T$)')
    plt.legend()
    plt.savefig(output_img)

if __name__ == "__main__":
    run_sigma_test()

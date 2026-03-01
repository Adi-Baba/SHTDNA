
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def predict_nucleosome_exclusion():
    output_dir = Path(r"d:\OnlyHST\SHTDNA\PredictionMode")
    output_log = output_dir / "nucleosome_prediction_log.txt"
    output_img = output_dir / "nucleosome_energy_landscape.png"
    
    # 1. The Geometric Constitutive Law (from Thermodynamics Mode)
    # Lp(eta) = L_coil * (1 - eta)^(-gamma)
    L_coil = 0.8946  # nm
    gamma = 1.5503   # dimensionless
    
    def geometric_stiffness(eta):
        # Avoid singularity at 1.0
        return L_coil * (1.0 - np.minimum(eta, 0.999))**(-gamma)
        
    # 2. Physics of Nucleosome Wrapping
    # Energy = 1/2 * Lp * kTs * Integral( curvature^2 ) ds
    # For a circle of Radius R and Length L:
    # Curvature k = 1/R. Integral = L * (1/R)^2 = L / R^2.
    # E_wrap = 1/2 * Lp * (L / R^2)  (in kB T units)
    
    R_histone = 4.3 # nm (Superhelix radius)
    L_wrap = 50.0   # nm (~147 bp)
    
    bending_factor = 0.5 * L_wrap / (R_histone**2) # ~ 1.35 nm^-1
    
    # 3. Compare Scenarios
    # Scenario A: Generic B-DNA
    # Mean Eta ~ 0.90 from our cloud
    eta_generic = 0.90
    
    # Scenario B: A-tract / Exclusion Sequence
    # Mean Eta ~ 0.96 from our cloud/trajectory
    eta_flat = 0.96
    
    Lp_gen = geometric_stiffness(eta_generic)
    E_gen = Lp_gen * bending_factor
    
    Lp_flat = geometric_stiffness(eta_flat)
    E_flat = Lp_flat * bending_factor
    
    dE = E_flat - E_gen
    occupancy_ratio = np.exp(-dE) # Boltzmann factor
    
    # 4. Generate "Occupancy Landscape" Plot
    eta_range = np.linspace(0.85, 0.98, 100)
    Lp_vals = geometric_stiffness(eta_range)
    E_vals = Lp_vals * bending_factor
    occ_vals = np.exp(-(E_vals - E_gen)) # Normalized to Generic = 1.0
    
    plt.figure(figsize=(10, 6))
    
    # Dual Axis Plot
    ax1 = plt.gca()
    ax1.plot(eta_range, E_vals, 'r-', linewidth=3, label='Wrapping Energy Cost')
    ax1.set_xlabel('Shambhian Twist $\eta_r$ (Flatness)', fontsize=12)
    ax1.set_ylabel('Energy Cost ($k_B T$)', color='red', fontsize=12)
    ax1.tick_params(axis='y', labelcolor='red')
    
    # Annotate Points
    ax1.scatter([eta_generic], [E_gen], color='black', s=100, zorder=5)
    ax1.text(eta_generic, E_gen + 5, "Generic B-DNA\n(Flexible)", ha='center')
    
    ax1.scatter([eta_flat], [E_flat], color='black', s=100, zorder=5)
    ax1.text(eta_flat, E_flat + 5, "Exclusion Seq\n(Rigid Ribbon)", ha='center')
    
    ax2 = ax1.twinx()
    ax2.plot(eta_range, occ_vals, 'b--', linewidth=2, label='Occupancy Probability')
    ax2.set_ylabel('Relative Occupancy ($P/P_{gen}$)', color='blue', fontsize=12)
    ax2.tick_params(axis='y', labelcolor='blue')
    ax2.set_yscale('log')
    
    plt.title('Prediction: Geometric Control of Nucleosome Occupancy', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.savefig(output_img, dpi=150)
    
    # 5. Log Results
    with open(output_log, "w") as f:
        f.write("# Prediction: Nucleosome Exclusion Mechanism\n\n")
        f.write(f"Generic DNA (Eta={eta_generic:.2f}):\n")
        f.write(f"  Stiffness Lp: {Lp_gen:.2f} nm\n")
        f.write(f"  Wrapping Energy: {E_gen:.2f} kB T\n\n")
        
        f.write(f"Flat/Exclusion DNA (Eta={eta_flat:.2f}):\n")
        f.write(f"  Stiffness Lp: {Lp_flat:.2f} nm\n")
        f.write(f"  Wrapping Energy: {E_flat:.2f} kB T\n\n")
        
        f.write(f"Energy Penalty (dE): {dE:.2f} kB T\n")
        f.write(f"Predicted Occupancy Ratio: {occupancy_ratio:.2e}\n")
        f.write("\nCONCLUSION:\n")
        if dE > 5.0:
            f.write("STRONG EXCLUSION PREDICTED.\n")
            f.write("Flat DNA is mechanically prohibited from forming nucleosomes.\n")
            f.write("This matches the biological function of Poly(dA:dT) tracts.\n")
        else:
            f.write("WEAK EFFECT. Geometry implies only minor bias.\n")
            
    print("Prediction Complete.")
    print(f"Generic Energy: {E_gen:.2f} kT")
    print(f"Flat Energy:    {E_flat:.2f} kT")
    print(f"Penalty:        {dE:.2f} kT")

if __name__ == "__main__":
    predict_nucleosome_exclusion()

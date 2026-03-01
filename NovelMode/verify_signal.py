
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from Bio.PDB import MMCIFParser, NeighborSearch

# Reuse existing logic
sys.path.append(r"d:\OnlyHST\SHTDNA\NovelMode")
import sht_dna_analysis

PDB_ID = "4RTK"
STRUCTURE_DIR = Path(r"d:\OnlyHST\SHTDNA\PDB\dataset\hard_dataset\structures")
CIF_PATH = STRUCTURE_DIR / f"{PDB_ID}.cif"

def get_water_coords(structure):
    waters = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip() in ['HOH', 'WAT', 'SOL']:
                    for atom in residue:
                        if atom.element == 'O':
                            waters.append(atom.get_coord())
    return np.array(waters)

def run_verification():
    print(f"Verifying signal for {PDB_ID}...")
    
    # 1. Compute Profiles
    res_sht = sht_dna_analysis.compute_sigma_sht_from_pdb(str(CIF_PATH), num_slices=100)
    h = res_sht['h']
    omega_h = res_sht['omega_h']
    ell = res_sht['ell']
    
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(PDB_ID, str(CIF_PATH))
    water_coords = get_water_coords(structure)
    
    water_h = water_coords @ ell
    water_counts = np.zeros_like(h)
    dz = (h[-1] - h[0]) / (len(h) - 1)
    
    for wh in water_h:
        idx = (np.abs(h - wh)).argmin()
        if np.abs(h[idx] - wh) < dz:
            water_counts[idx] += 1
            
    # Smooth water
    window = 3
    water_smooth = np.convolve(water_counts, np.ones(window)/window, mode='same')
    
    # Mask edges
    mask = (h > h[0] + 5) & (h < h[-1] - 5)
    
    x = omega_h[mask]
    y = water_smooth[mask]
    h_valid = h[mask]
    
    obs_corr = np.corrcoef(x, y)[0,1]
    print(f"Observed Correlation: {obs_corr:.5f}")
    
    # 2. Permutation Test
    n_perms = 100000
    better_count = 0
    
    print(f"Running {n_perms} permutations...")
    y_shuff = y.copy()
    for _ in range(n_perms):
        np.random.shuffle(y_shuff)
        r = np.corrcoef(x, y_shuff)[0,1]
        # Two-tailed test
        if abs(r) >= abs(obs_corr):
            better_count += 1
            
    p_value = better_count / n_perms
    print(f"P-value: {p_value:.5f}")
    
    if p_value < 0.001:
        print("RESULT: Statistically Significant (Not Noise)")
    else:
        print("RESULT: Could be noise")

    # 3. Visualization
    fig, ax1 = plt.subplots(figsize=(10, 5))
    
    color = 'tab:red'
    ax1.set_xlabel('Height h (A)')
    ax1.set_ylabel('Roughness (Mamton)', color=color)
    ax1.plot(h_valid, x, color=color, linewidth=2, label='Roughness')
    ax1.tick_params(axis='y', labelcolor=color)
    
    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel('Water Density', color=color)
    ax2.plot(h_valid, y, color=color, linewidth=2, linestyle='--', label='Hydration')
    ax2.tick_params(axis='y', labelcolor=color)
    
    plt.title(f"SHT Geometric Gating of Hydration ({PDB_ID})\nR = {obs_corr:.3f}, P = {p_value:.5f}")
    fig.tight_layout()
    plt.savefig("d:\\OnlyHST\\SHTDNA\\NovelMode\\verification_plot.png")
    print("Plot saved to verification_plot.png")

if __name__ == "__main__":
    run_verification()

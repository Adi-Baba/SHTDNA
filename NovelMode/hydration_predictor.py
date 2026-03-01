

import sys
import os
import numpy as np
import pandas as pd
from pathlib import Path
from Bio.PDB import MMCIFParser

import sht_dna_analysis

STRUCTURES_DIR = Path(r"d:\OnlyHST\SHTDNA\PDB\dataset\hard_dataset\structures")
PDB_LIST_FILE = Path(r"d:\OnlyHST\SHTDNA\NovelMode\high_res_pdbs.txt")
OUTPUT_CSV = Path(r"d:\OnlyHST\SHTDNA\NovelMode\local_hydration_results.csv")

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

def compute_local_correlation(pdb_id):
    cif_path = STRUCTURES_DIR / f"{pdb_id}.cif"
    if not cif_path.exists():
        return None

    try:
        # 1. SHT Analysis
        res_sht = sht_dna_analysis.compute_sigma_sht_from_pdb(str(cif_path), num_slices=100)
        
        h = res_sht['h'] # Slice centers
        omega_h = res_sht['omega_h'] # Roughness profile
        eta_h = res_sht['epsilon_h'] # Flatness profile
        ell = res_sht['ell'] # Principal axis
        
        if len(h) < 10: return None
        
        # 2. Water Analysis
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(pdb_id, str(cif_path))
        water_coords = get_water_coords(structure)
        
        if len(water_coords) == 0:
            return None
            
        # Project waters onto axis
        # We need "center" of DNA to shift waters? 
        # SHT heights are coords @ ell.
        # So water heights are water_coords @ ell.
        
        water_h = water_coords @ ell
        
        # Bin waters into slices
        # Slices are defined by bin edges approx midway between h
        # Let's just use simple nearest slice assignment for speed
        
        water_counts = np.zeros_like(h)
        dz = (h[-1] - h[0]) / (len(h) - 1)
        
        for wh in water_h:
            # Find nearest slice index
            # h is sorted? Usually yes.
            idx = (np.abs(h - wh)).argmin()
            # Check if within reasonable slab thickness (e.g., +/- dz)
            if np.abs(h[idx] - wh) < dz:
                water_counts[idx] += 1
                
        # Normalize water counts (density) -> not strictly necessary for correlation but good for stats
        # Smoothing water counts? Waters are discrete and sparse. 
        # Apply a rolling window to smoothing both signals to match scales?
        # SHT is already somewhat smooth (kernel_width=3.4A).
        # Let's smooth water counts with same window.
        
        # Simple moving average
        window = 3
        water_smooth = np.convolve(water_counts, np.ones(window)/window, mode='same')
        
        # 3. Compute Local Correlation for this molecule
        # Remove edges where density might be weird
        mask = (h > h[0] + 5) & (h < h[-1] - 5)
        if mask.sum() < 5: return None
        
        corr_omega = np.corrcoef(omega_h[mask], water_smooth[mask])[0, 1]
        corr_eta = np.corrcoef(eta_h[mask], water_smooth[mask])[0, 1]
        
        return {
            'pdb_id': pdb_id,
            'corr_roughness_water': corr_omega,
            'corr_flatness_water': corr_eta,
            'mean_roughness': np.mean(omega_h[mask]),
            'mean_water': np.mean(water_smooth[mask])
        }

    except Exception as e:
        # print(f"Err {pdb_id}: {e}")
        return None

def main():
    if not PDB_LIST_FILE.exists(): return
    with open(PDB_LIST_FILE, 'r') as f:
        pdb_ids = [line.strip() for line in f if line.strip()]

    print(f"Processing local profiles for {len(pdb_ids)} structures...")
    results = []
    
    for i, pid in enumerate(pdb_ids):
        res = compute_local_correlation(pid)
        if res:
            results.append(res)
        if i % 10 == 0: print(f"{i}/{len(pdb_ids)}", end='\r')
        
    df = pd.DataFrame(results)
    df.to_csv(OUTPUT_CSV, index=False)
    
    print("\n--- Local Profile Correlations (Meta-Analysis) ---")
    print(f"Mean Correlation (Roughness vs Water): {df['corr_roughness_water'].mean():.4f}")
    print(f"Mean Correlation (Flatness vs Water):  {df['corr_flatness_water'].mean():.4f}")
    
    # Save significant findings
    sig = df[df['corr_roughness_water'] < -0.5]
    print(f"Number of strong negative correlations (<-0.5): {len(sig)}")

if __name__ == "__main__":
    main()

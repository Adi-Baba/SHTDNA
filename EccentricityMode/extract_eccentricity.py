
import sys
import os
import numpy as np
import pandas as pd
from pathlib import Path

# Add legacy code path
sys.path.append(r"d:\OnlyHST\SHTDNA\PDB\code")
import sht_dna_analysis

def analyze_eccentricity(pdb_path):
    try:
        results = sht_dna_analysis.compute_sigma_sht_from_pdb(
            str(pdb_path), num_slices=100
        )
        
        epsilon_h = results['epsilon_h']
        h = results['h']
        
        # Load sequence (A-tract detection)
        from Bio.PDB import PDBParser, MMCIFParser
        from Bio.SeqUtils import seq1
        parser = MMCIFParser(QUIET=True) if Path(pdb_path).suffix == '.cif' else PDBParser(QUIET=True)
        structure = parser.get_structure('X', str(pdb_path))
        seq_record = ""
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] == " ":
                        seq_record += seq1(residue.get_resname())
            break 
            
        eps_atract = []
        eps_other = []
        
        n_bases = len(seq_record)
        if n_bases < 4: return None
        
        h_min, h_max = np.min(h), np.max(h)
        total_height = h_max - h_min
        
        for i in range(len(seq_record)):
            is_atract = False
            if seq_record[i:i+4] == "AAAA" or seq_record[i:i+4] == "TTTT": is_atract = True
            elif i>=1 and (seq_record[i-1:i+3] == "AAAA" or seq_record[i-1:i+3] == "TTTT"): is_atract = True
            elif i>=2 and (seq_record[i-2:i+2] == "AAAA" or seq_record[i-2:i+2] == "TTTT"): is_atract = True
            elif i>=3 and (seq_record[i-3:i+1] == "AAAA" or seq_record[i-3:i+1] == "TTTT"): is_atract = True
            
            relative_pos = i / (n_bases - 1)
            target_h = h_min + relative_pos * total_height
            idx = (np.abs(h - target_h)).argmin()
            val = epsilon_h[idx]
            
            if is_atract:
                eps_atract.append(val)
            else:
                eps_other.append(val)

        return {
            'global_stats': {
                'pdb_id': Path(pdb_path).stem,
                'mean_epsilon': np.mean(epsilon_h),
                'max_epsilon': np.max(epsilon_h),
                'std_epsilon': np.std(epsilon_h),
                'len_h': results['L']
            },
            'eps_atract': eps_atract,
            'eps_other': eps_other
        }

    except Exception as e:
        return None

def main():
    dataset_dir = Path(r"d:\OnlyHST\SHTDNA\PDB\dataset\hard_dataset\structures")
    output_csv = r"d:\OnlyHST\SHTDNA\EccentricityMode\eccentricity_results.csv"
    
    all_pdbs = list(dataset_dir.glob("*.cif")) + list(dataset_dir.glob("*.pdb"))
    
    data = []
    eps_atract_pool = []
    eps_other_pool = []
    
    count = 0
    limit = 200 
    
    print(f"Starting A-tract Validation for Eccentricity on {limit} structures...")
    
    for i, pdb_file in enumerate(all_pdbs):
        if count >= limit: break
            
        res = analyze_eccentricity(pdb_file)
        if res:
            data.append(res['global_stats'])
            eps_atract_pool.extend(res['eps_atract'])
            eps_other_pool.extend(res['eps_other'])
            count += 1
            if count % 10 == 0: print(f".", end='', flush=True)
            
    print("\nExtraction Complete.")
    df = pd.DataFrame(data)
    df.to_csv(output_csv, index=False)
    
    # --- Validation Stats ---
    print("\n--- Validation Evidence (Eccentricity) ---")
    print(f"Total A-tract sites: {len(eps_atract_pool)}")
    print(f"Total Non-A-tract sites: {len(eps_other_pool)}")
    
    if len(eps_atract_pool) > 10 and len(eps_other_pool) > 10:
        mean_atract = np.mean(eps_atract_pool)
        mean_other = np.mean(eps_other_pool)
        
        from scipy.stats import ttest_ind
        t_stat, p_val = ttest_ind(eps_atract_pool, eps_other_pool, equal_var=False)
        
        print(f"Mean Epsilon (A-tracts): {mean_atract:.4f}")
        print(f"Mean Epsilon (Others):   {mean_other:.4f}")
        print(f"Difference: {mean_atract - mean_other:.4f}")
        print(f"P-value: {p_val:.4e}")
        
        with open(r"d:\OnlyHST\SHTDNA\EccentricityMode\validation_epsilon.txt", "w") as f:
            f.write(f"Mean Epsilon (A-tracts): {mean_atract:.4f}\n")
            f.write(f"Mean Epsilon (Others):   {mean_other:.4f}\n")
            f.write(f"Difference: {mean_atract - mean_other:.4f}\n")
            f.write(f"P-value: {p_val:.4e}\n")
            
            if p_val < 0.05:
                print("EVIDENCE FOUND: A-tracts have statistically distinct shape.")
                f.write("Conclusion: EVIDENCE FOUND\n")
            else:
                print("NO EVIDENCE: Shape indistinguishable from noise.")
                f.write("Conclusion: NO EVIDENCE\n")

if __name__ == "__main__":
    main()

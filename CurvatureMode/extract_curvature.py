
import sys
import os
import numpy as np
import pandas as pd
from pathlib import Path

# Add legacy code path
sys.path.append(r"d:\OnlyHST\SHTDNA\PDB\code")
import sht_dna_analysis

def analyze_curvature(pdb_path):
    try:
        results = sht_dna_analysis.compute_sigma_sht_from_pdb(
            str(pdb_path), num_slices=100
        )
        
        kappa_h = results['kappa_h'] # Unit: Vij (nm^-1 or similar)
        h = results['h']
        sigma_global = results['sigma_SHT_global']
        
        # Load sequence (needed for A-tracts)
        # Using the same logic as breathing extraction (re-using parser would be better but this is quick)
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
            
        # Detect A-tracts (4+ A or T)
        # We need to map discrete sequence to continuous properties kappa(h) same as before.
        kappa_atract = []
        kappa_other = []
        
        n_bases = len(seq_record)
        if n_bases < 4: return None
        
        h_min, h_max = np.min(h), np.max(h)
        total_height = h_max - h_min
        
        for i in range(len(seq_record)):
            # Define if this base is part of an A-tract (simplified detection)
            # Check window of 4 around i
            is_atract = False
            # Check for AAAA or TTTT
            if seq_record[i:i+4] == "AAAA" or seq_record[i:i+4] == "TTTT":
                is_atract = True
            # Also check backward (if we are essentially inside one)
            elif i>=1 and (seq_record[i-1:i+3] == "AAAA" or seq_record[i-1:i+3] == "TTTT"): is_atract = True
            elif i>=2 and (seq_record[i-2:i+2] == "AAAA" or seq_record[i-2:i+2] == "TTTT"): is_atract = True
            elif i>=3 and (seq_record[i-3:i+1] == "AAAA" or seq_record[i-3:i+1] == "TTTT"): is_atract = True
            
            # Map to h
            relative_pos = i / (n_bases - 1)
            target_h = h_min + relative_pos * total_height
            idx = (np.abs(h - target_h)).argmin()
            val = kappa_h[idx]
            
            if is_atract:
                kappa_atract.append(val)
            else:
                kappa_other.append(val)

        return {
            'global_stats': {
                'pdb_id': Path(pdb_path).stem,
                'mean_vij': np.mean(kappa_h),
                'max_vij': np.max(kappa_h),
                'std_vij': np.std(kappa_h),
                'len_h': results['L']
            },
            'vij_atract': kappa_atract,
            'vij_other': kappa_other
        }

    except Exception as e:
        # print(f"Error processing {pdb_path}: {e}")
        return None

def main():
    dataset_dir = Path(r"d:\OnlyHST\SHTDNA\PDB\dataset\hard_dataset\structures")
    output_csv = r"d:\OnlyHST\SHTDNA\CurvatureMode\curvature_vij_results.csv"
    
    all_pdbs = list(dataset_dir.glob("*.cif")) + list(dataset_dir.glob("*.pdb"))
    
    data = []
    vij_atract_pool = []
    vij_other_pool = []
    
    count = 0
    limit = 200
    
    print(f"Starting A-tract Validation for Unit 'Vij' on {limit} structures...")
    
    for i, pdb_file in enumerate(all_pdbs):
        if count >= limit: break
            
        res = analyze_curvature(pdb_file)
        if res:
            data.append(res['global_stats'])
            vij_atract_pool.extend(res['vij_atract'])
            vij_other_pool.extend(res['vij_other'])
            count += 1
            if count % 10 == 0: print(f".", end='', flush=True)
            
    print("\nExtraction Complete.")
    df = pd.DataFrame(data)
    df.to_csv(output_csv, index=False)
    
    # --- Validation Stats ---
    print("\n--- Validation Evidence (Unit: Vij) ---")
    print(f"Total A-tract sites: {len(vij_atract_pool)}")
    print(f"Total Non-A-tract sites: {len(vij_other_pool)}")
    
    if len(vij_atract_pool) > 10 and len(vij_other_pool) > 10:
        mean_atract = np.mean(vij_atract_pool)
        mean_other = np.mean(vij_other_pool)
        
        from scipy.stats import ttest_ind
        t_stat, p_val = ttest_ind(vij_atract_pool, vij_other_pool, equal_var=False)
        
        print(f"Mean Curvature (A-tracts): {mean_atract:.4f} Vij")
        print(f"Mean Curvature (Others):   {mean_other:.4f} Vij")
        print(f"Difference: {mean_atract - mean_other:.4f} Vij")
        print(f"P-value: {p_val:.4e}")
        
        if p_val < 0.05:
            print("EVIDENCE FOUND: A-tracts have statistically distinct curvature.")
            conclusion = "EVIDENCE FOUND"
        else:
            print("NO EVIDENCE: Curvature indistinguishable from noise.")
            conclusion = "NO EVIDENCE"
            
        with open(r"d:\OnlyHST\SHTDNA\CurvatureMode\validation_vij.txt", "w") as f:
            f.write(f"Mean Curvature (A-tracts): {mean_atract:.4f} Vij\n")
            f.write(f"Mean Curvature (Others):   {mean_other:.4f} Vij\n")
            f.write(f"Difference: {mean_atract - mean_other:.4f} Vij\n")
            f.write(f"P-value: {p_val:.4e}\n")
            f.write(f"Conclusion: {conclusion}\n")

if __name__ == "__main__":
    main()

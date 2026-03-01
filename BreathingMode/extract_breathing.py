
import sys
import os
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import PDB
from Bio.SeqUtils import seq1

# Make import of local analysis code project-relative (cross-platform)
PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / 'PDB' / 'code'))
import sht_dna_analysis

def get_dna_sequence(structure):
    """
    Extract DNA sequence from a Biopython structure object.
    Returns a single string of the combined sequence (if duplex).
    Effectively just concatenates all nucleotide residues.
    """
    seq = ""
    # Standard DNA residue names
    dna_res = {'DA', 'DT', 'DC', 'DG', 'A', 'T', 'C', 'G'}
    
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                if resname in dna_res:
                    # Convert 3-letter (DA, DT) to 1-letter
                    if len(resname) == 2 and resname[0] == 'D':
                        seq += resname[1]
                    else:
                        seq += resname
        break # Only first model
    return seq

def analyze_breathing_mode(pdb_path):
    """
    Run SHT analysis and compute Breathing Mode statistics.
    Returns a dict of metrics.
    """
    try:
        # Run SHT
        results = sht_dna_analysis.compute_sigma_sht_from_pdb(
            str(pdb_path), num_slices=100
        )
        
        beta_h = results['beta_h']
        h_vals = results['h']
        
        if len(beta_h) == 0:
            return None

        # --- Sequence Extraction ---
        if str(pdb_path).lower().endswith('.cif'):
            parser = PDB.MMCIFParser(QUIET=True)
        else:
            parser = PDB.PDBParser(QUIET=True)
            
        structure = parser.get_structure('temp', str(pdb_path))
        sequence = get_dna_sequence(structure)
        
        if len(sequence) == 0:
            return None

        # --- Site-Specific Mapping ---
        # We need to map the continuous h-profile beta(h) to discrete base pairs.
        # h spans from ~min(h) to ~max(h). Sequence has N bases.
        # We assume linear mapping for B-DNA (3.4A rise).
        
        beta_at_site = []
        beta_gc_site = []
        
        n_bases = len(sequence)
        
        at_count = sequence.count('A') + sequence.count('T')
        total_len = len(sequence)
        at_content = at_count / total_len if total_len > 0 else 0
        
        # Approximate h-range covered by measuring atoms
        h_min, h_max = np.min(h_vals), np.max(h_vals)
        total_height = h_max - h_min
        
        # Check if length matches roughly (allow some tolerance for end effects)
        expected_len = (n_bases - 1) * 3.4
        
        # Iterate through sequence indices
        for i, char in enumerate(sequence):
            # Estimated height of this base pair
            # We assume the first base is at the bottom (or top? PCA direction is arbitrary).
            # But Beta corresponds to Width, so detailed alignment isn't critical for *distributions*,
            # provided we sample vaguely evenly.
            # Better approach: Just grid the beta_h profile into N bins.
            
            # Simple interpolation
            relative_pos = i / (n_bases - 1) if n_bases > 1 else 0.5
            target_h = h_min + relative_pos * total_height
            
            # Find closest h index
            idx = (np.abs(h_vals - target_h)).argmin()
            val = beta_h[idx]
            
            if char in ['A', 'T']:
                beta_at_site.append(val)
            elif char in ['G', 'C']:
                beta_gc_site.append(val)
        
        return {
            'global_stats': {
                'pdb_id': Path(pdb_path).stem,
                'mean_beta': np.mean(beta_h),
                'std_beta': np.std(beta_h),
                'max_beta': np.max(beta_h),
                'min_beta': np.min(beta_h),
                'at_content': at_content,
                'seq_len': total_len
            },
            'beta_at': beta_at_site,
            'beta_gc': beta_gc_site
        }

    except Exception as e:
        print(f"Error processing {pdb_path}: {e}", flush=True)
        import traceback
        traceback.print_exc()
        return None

def main():
    parser = argparse.ArgumentParser(description='Extract breathing-mode statistics from PDB/CIF set')
    parser.add_argument('--dataset-dir', type=str,
                        default=str(PROJECT_ROOT / 'PDB' / 'dataset' / 'hard_dataset' / 'structures'),
                        help='Directory containing PDB/CIF structures')
    parser.add_argument('--out-dir', type=str, default=str(PROJECT_ROOT / 'BreathingMode'),
                        help='Output directory for CSV and site data')
    parser.add_argument('--limit', type=int, default=200, help='Max number of structures to process')
    args = parser.parse_args()

    dataset_dir = Path(args.dataset_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    output_csv = out_dir / 'breathing_results.csv'
    site_data_file = out_dir / 'breathing_site_data.npz'
    significance_file = out_dir / 'significance_test.txt'

    # Process a subset or all
    all_pdbs = list(dataset_dir.glob("*.cif"))
    print(f"Found {len(all_pdbs)} PDBs/CIFs.")

    data = []

    # Run on up to `limit` structures
    count = 0
    limit = args.limit

    # Collect site-specific statistics
    beta_at_values = []
    beta_gc_values = []

    for i, pdb_file in enumerate(all_pdbs):
        if count >= limit:
            break

        # Verbose progress
        print(f"[{count+1}/{limit}] Processing {pdb_file.stem}...", end='', flush=True)

        res = analyze_breathing_mode(pdb_file)
        if res:
            data.append(res['global_stats'])
            beta_at_values.extend(res['beta_at'])
            beta_gc_values.extend(res['beta_gc'])
            # Save per-structure site arrays for downstream bootstrapping and mixed-models
            site_dir = out_dir / 'site_data'
            site_dir.mkdir(parents=True, exist_ok=True)
            pdb_key = res['global_stats']['pdb_id']
            np.savez_compressed(site_dir / f"{pdb_key}_sites.npz",
                                beta_at=np.array(res['beta_at']),
                                beta_gc=np.array(res['beta_gc']))
            count += 1
            print(" Done.", flush=True)
        else:
            print(" Skipped (No valid data/Beta).", flush=True)

    df = pd.DataFrame(data)
    df.to_csv(output_csv, index=False)
    print(f"Saved global results to {output_csv}")

    # --- Site-Specific T-Test ---
    print("\n--- Site-Specific Analysis ---")
    print(f"Collected {len(beta_at_values)} A-T sites and {len(beta_gc_values)} G-C sites.")
    from scipy.stats import ttest_ind

    if len(beta_at_values) > 0 and len(beta_gc_values) > 0:
        mean_at = np.mean(beta_at_values)
        mean_gc = np.mean(beta_gc_values)
        std_at = np.std(beta_at_values)
        std_gc = np.std(beta_gc_values)

        # Test for difference in MEANS
        t_stat, p_val_t = ttest_ind(beta_at_values, beta_gc_values, equal_var=False)

        # Test for difference in VARIANCES (Levene's Test)
        from scipy.stats import levene
        stat_levene, p_val_levene = levene(beta_at_values, beta_gc_values)

        # Save a small significance summary and also persist the raw site-level arrays
        np.savez_compressed(site_data_file, beta_at=np.array(beta_at_values), beta_gc=np.array(beta_gc_values))

        with open(significance_file, 'w') as f:
            f.write(f"Mean Beta (A-T): {mean_at:.4f} +/- {std_at:.4f}\n")
            f.write(f"Mean Beta (G-C): {mean_gc:.4f} +/- {std_gc:.4f}\n")
            f.write(f"Mean Diff P-value (T-test): {p_val_t:.4e}\n")
            f.write(f"Variance Diff P-value (Levene): {p_val_levene:.4e}\n")
            if p_val_levene < 0.05:
                f.write("RESULT: Statistically Significant Difference in VARIANCE found!\n")
            else:
                f.write("RESULT: No significant difference in variance.\n")

        print(f"Variance P-value (Levene): {p_val_levene:.4e}")
        print(f"Mean Beta (G-C): {mean_gc:.4f} +/- {std_gc:.4f}")
        print(f"Difference: {mean_at - mean_gc:.4f}")
        print(f"T-statistic: {t_stat:.4f}")
        print(f"P-value (Means): {p_val_t:.4e}")

        print(f"Saved site-level arrays to {site_data_file}")
        print(f"Saved significance summary to {significance_file}")

        if p_val_t < 0.05:
            print("RESULT: Statistically Significant MEAN Difference found!")
        else:
            print("RESULT: No significant MEAN difference.")

        if p_val_levene < 0.05:
            print("RESULT: Statistically Significant VARIANCE Difference found!")
        else:
            print("RESULT: No significant VARIANCE difference.")

    # Quick Correlation Check
    if not df.empty:
        corr = df['at_content'].corr(df['std_beta'])
        print(f"\nCorrelation (AT-Content vs Breathing Amplitude): {corr:.4f}")


if __name__ == "__main__":
    main()

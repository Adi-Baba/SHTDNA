#!/usr/bin/env python3
"""
batch_analyze.py

Run SHT analysis on a curated list of DNA duplex PDB files and
aggregate the results into a single data file.

This script forms the final step of the automated workflow:
1. Download (pdb_downloader.py)
2. Filter (filter_and_prepare.py)
3. Analyze (this script)
"""

import os
import json
import pickle
import argparse
import pandas as pd
from pathlib import Path

from sht_dna_analysis import compute_sigma_sht_from_pdb

def main():
    # --- Define robust paths based on script location ---
    SCRIPT_DIR = Path(__file__).parent.resolve()
    PROJECT_ROOT = SCRIPT_DIR.parent
    DATA_DIR = PROJECT_ROOT / "dataset" / "hard_dataset"
    # Define a default output directory at the project root for analysis results
    DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "sht_analysis_output"

    parser = argparse.ArgumentParser(
        description="Batch run SHT analysis on a curated set of PDB files."
    )
    parser.add_argument(
        '--metadata_file',
        type=Path,
        default=DATA_DIR / 'curated_nonredundant_duplex_metadata.json',
        help="Path to the curated metadata JSON file from filter_and_prepare.py"
    )
    parser.add_argument(
        '--pdb_dir',
        type=Path,
        default=DATA_DIR / 'structures',
        help="Directory where the mmCIF files are stored."
    )
    parser.add_argument(
        '--output_csv',
        type=Path,
        default=DEFAULT_OUTPUT_DIR / 'sht_analysis_results.csv',
        help="Path to save the aggregated results CSV file."
    )
    parser.add_argument(
        '--output_pickle_dir',
        type=Path,
        default=DEFAULT_OUTPUT_DIR / 'sht_local_profiles',
        help="Directory to save detailed local profile pickle files."
    )
    args = parser.parse_args()

    if not args.metadata_file.exists():
        print(f"Error: Curated metadata file not found at {args.metadata_file}")
        return

    # Ensure output directories exist before we start writing files
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_pickle_dir.mkdir(parents=True, exist_ok=True)

    with open(args.metadata_file, 'r') as f:
        curated_metadata = json.load(f)

    all_results = []
    processed_pdbs = set()

    # --- Resume from Checkpoint ---
    if args.output_csv.exists():
        print(f"Found existing results file at {args.output_csv}. Resuming...")
        try:
            existing_df = pd.read_csv(args.output_csv)
            if not existing_df.empty and 'pdb_id' in existing_df.columns:
                all_results = existing_df.to_dict('records')
                processed_pdbs = set(existing_df['pdb_id'])
                print(f"Loaded {len(processed_pdbs)} previously analyzed structures.")
        except (pd.errors.EmptyDataError, KeyError):
            print(f"  [WARN] Existing results file at {args.output_csv} is empty or malformed. Starting from scratch.")

    print(f"Starting SHT analysis for {len(curated_metadata)} structures...")

    # --- Configuration for Checkpointing ---
    CHECKPOINT_INTERVAL = 250
    new_results_count = 0

    for i, (pdb_id, meta) in enumerate(curated_metadata.items(), 1):
        if pdb_id in processed_pdbs:
            continue

        print(f"[{i}/{len(curated_metadata)}] Analyzing {pdb_id}...")
        pdb_path = args.pdb_dir / f"{pdb_id}.cif"
        chain_ids = meta.get('chain_ids')

        if not pdb_path.exists() or not chain_ids:
            print(f"  [WARN] Skipping {pdb_id}: File or chain IDs missing.")
            continue

        try:
            sht_results = compute_sigma_sht_from_pdb(
                pdb_path=str(pdb_path),
                chain_ids=chain_ids
            )

            # Aggregate the important scalar results
            record = {
                'pdb_id': pdb_id,
                'resolution': meta.get('resolution'),
                'method': meta.get('method'),
                **{k: sht_results[k] for k in ['sigma_SHT_global', 'Tw_SHT', 'Wr_SHT', 'Tw0', 'L']}
            }
            all_results.append(record)
            new_results_count += 1

            # Save the full, detailed results dictionary (including local profiles)
            # to a pickle file. This is better than CSV for storing numpy arrays.
            pickle_path = args.output_pickle_dir / f"{pdb_id}_profile.pkl"
            with open(pickle_path, 'wb') as pf:
                pickle.dump(sht_results, pf)
            print(f"  [INFO] Saved local profile to {pickle_path}")

        except Exception as e:
            print(f"  [ERROR] Failed to analyze {pdb_id}: {e}")

        # --- Checkpoint Saving Logic ---
        # Periodically save the results to disk to prevent data loss on long runs.
        if new_results_count > 0 and new_results_count % CHECKPOINT_INTERVAL == 0:
            print(f"\n--- Checkpoint: Saving {len(all_results)} total results to {args.output_csv} ---\n")
            results_df = pd.DataFrame(all_results)
            results_df.to_csv(args.output_csv, index=False)

    # Save aggregated results to a CSV file
    if new_results_count > 0:
        results_df = pd.DataFrame(all_results)
        results_df.to_csv(args.output_csv, index=False)
        print(f"\nAnalysis complete. Aggregated results saved to {args.output_csv}")
    else:
        print("\nAnalysis complete. No new structures were analyzed.")

if __name__ == "__main__":
    main()
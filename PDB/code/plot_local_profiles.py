#!/usr/bin/env python3
"""
plot_local_profiles.py

Load detailed SHT analysis results from pickle files and plot the local
supercoiling profiles for specified PDB IDs. This allows for a detailed
inspection of torsional hotspots and other local geometric features.
"""

import argparse
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# --- Plotting Style ---
sns.set_theme(style="whitegrid")

def plot_local_profile(pdb_id, results_data, output_dir):
    """
    Generates and saves a plot of the local supercoiling profile for a single PDB.
    """
    h = results_data.get("h")
    sigma_local = results_data.get("sigma_local_h")
    rho_twist = results_data.get("rho_twist_h")
    rho_writhe = results_data.get("rho_sh_h")
    rho_0 = results_data.get("rho_0")

    if h is None or sigma_local is None or rho_twist is None or rho_writhe is None:
        print(f"  [WARN] Skipping {pdb_id}: Missing local profile data.")
        return

    # Check for sufficient data points to plot
    if len(h) < 2:
        print(f"  [WARN] Skipping {pdb_id}: Not enough data points to create a meaningful plot.")
        return

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    fig.suptitle(f'Local Supercoiling Analysis for {pdb_id}', fontsize=18, weight='bold')

    # --- Plot 1: Local Supercoiling Density (sigma_local) ---
    ax1.plot(h, sigma_local, label=r'$\sigma_{SHT}(h)$', color='r', linewidth=2)
    ax1.axhline(0, color='grey', linestyle='--', linewidth=1)
    ax1.set_ylabel(r'Local Supercoiling Density ($\sigma_{SHT}(h)$)', fontsize=12)
    ax1.set_title('Local Supercoiling Profile', fontsize=14)
    ax1.legend()
    ax1.grid(True, which='both', linestyle='--', linewidth=0.5)

    # --- Plot 2: Twist and Writhe Density Components ---
    ax2.plot(h, rho_twist, label=r'Twist Density ($\rho_{twist}$)', color='b', linestyle='-')
    ax2.plot(h, rho_writhe, label=r'Writhe Density ($\rho_{sh}$)', color='g', linestyle='-')
    if rho_0 is not None and not np.isnan(rho_0):
        ax2.axhline(rho_0, color='k', linestyle=':', linewidth=1.5, label=fr'Relaxed Twist Density ($\rho_0 \approx {rho_0:.3f}$)')
    ax2.set_xlabel('Position along Principal Axis (h, Å)', fontsize=12)
    ax2.set_ylabel('Density (turns / Å)', fontsize=12)
    ax2.set_title('Constituent Density Profiles', fontsize=14)
    ax2.legend()
    ax2.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust for suptitle
    output_path = output_dir / f"{pdb_id}_local_profile.png"
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"  -> Saved local profile plot to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Plot local SHT profiles from saved pickle files."
    )
    parser.add_argument(
        'pdb_ids',
        nargs='+',
        help="One or more PDB IDs to plot (e.g., 1BNA 1A2E 1DCG)."
    )
    parser.add_argument(
        '--pickle_dir',
        type=Path,
        # Default path points to the output location of batch_analyze.py
        default=Path(__file__).parent.parent / 'sht_analysis_output' / 'sht_local_profiles',
        help="Directory where the local profile pickle files are stored."
    )
    parser.add_argument(
        '--output_dir',
        type=Path,
        # Save plots to a dedicated folder within the main analysis output directory
        default=Path(__file__).parent.parent / 'sht_analysis_output' / 'local_profile_plots',
        help="Directory to save the generated plots."
    )
    args = parser.parse_args()

    # --- Validation and Setup ---
    if not args.pickle_dir.exists():
        print(f"Error: Pickle directory not found at {args.pickle_dir}")
        return
    args.output_dir.mkdir(exist_ok=True)

    # --- Plotting Loop ---
    print(f"Generating local profile plots for {len(args.pdb_ids)} PDB(s)...")
    for pdb_id in args.pdb_ids:
        pdb_id_upper = pdb_id.upper()
        pickle_path = args.pickle_dir / f"{pdb_id_upper}_profile.pkl"

        if not pickle_path.exists():
            print(f"  [ERROR] No pickle file found for {pdb_id_upper} at {pickle_path}")
            continue

        print(f"Processing {pdb_id_upper}...")
        try:
            with open(pickle_path, 'rb') as f:
                results_data = pickle.load(f)
            plot_local_profile(pdb_id_upper, results_data, args.output_dir)
        except Exception as e:
            print(f"  [ERROR] Failed to plot {pdb_id_upper}: {e}")

    print("\nPlotting complete.")


if __name__ == "__main__":
    main()
"""
analyze_single.py (formerly test_local_sigma.py)

A utility to run SHT analysis on specific PDB IDs on-demand and
save the results in a format compatible with the main analysis pipeline.
This is useful for analyzing canonical structures or specific outliers that
may have been excluded from the main non-redundant dataset.

Run with:
    python test_local_sigma.py 1BNA 1ANA 1AOI
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pickle
from pathlib import Path

from sht_dna_analysis import compute_sigma_sht_from_pdb


# ---------- CONFIG ----------

# Define paths relative to this script's location for robustness
SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = SCRIPT_DIR.parent
PDB_STRUCTURE_DIR = PROJECT_ROOT / "dataset" / "hard_dataset" / "structures"
OUTPUT_PICKLE_DIR = PROJECT_ROOT / "sht_analysis_output" / "sht_local_profiles"
OUTPUT_PLOT_DIR = PROJECT_ROOT / "sht_analysis_output" / "single_analysis_plots"
OUTPUT_PLOT_DIR.mkdir(exist_ok=True)
OUTPUT_PICKLE_DIR.mkdir(exist_ok=True)

# For convenience, pre-define chain IDs for common canonical structures.
# The main pipeline discovers these automatically, but this is for on-demand analysis.
CANONICAL_CHAINS = {
    "1BNA": ["A", "B"],
    "1ANA": ["A", "B"],
    "1AOI": ["I", "J"],
    "2DCG": ["A", "B"], # Z-DNA example
}


def plot_local_profiles(label, results, outdir):
    """
    Make a 4-panel plot:
      1) theta(h),
      2) phi(h),
      3) local densities rho_twist(h), rho_sh(h),
      4) local supercoiling sigma_SHT(h).
    """
    h = results["h"]
    theta = results["theta"]
    phi = results["phi"]
    rho_twist_h = results["rho_twist_h"]
    rho_sh_h = results["rho_sh_h"]
    sigma_local_h = results["sigma_local_h"]
    rho_0 = results["rho_0"]

    fig, axes = plt.subplots(4, 1, figsize=(7, 10), sharex=True)
    fig.suptitle(f"SHT Local Supercoiling Analysis - {label}", fontsize=14)

    # Panel 1: theta(h)
    ax = axes[0]
    ax.plot(h, theta, "o-", ms=3, lw=1)
    ax.set_ylabel("θ(h) [rad]")
    ax.set_title("Helical Phase vs. Height")
    ax.grid(True)

    # Panel 2: phi(h)
    ax = axes[1]
    ax.plot(h, phi, "o-", ms=3, lw=1, color="green")
    ax.set_ylabel("φ(h) [rad]")
    ax.set_title("Axis Winding vs. Height")
    ax.grid(True)

    # Panel 3: local densities
    ax = axes[2]
    ax.plot(h, rho_twist_h, "o-", ms=3, lw=1, label="Twist density ρ_twist(h)")
    ax.plot(h, rho_sh_h, "o-", ms=3, lw=1, label="Writhe-like density ρ_sh(h)", color="red")
    ax.axhline(rho_0, color="k", ls="--", label="Relaxed twist density ρ₀")
    ax.set_ylabel("Density [turns/Å]")
    ax.set_title("Local Twist and Writhe Densities")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True)

    # Panel 4: local supercoiling profile
    ax = axes[3]
    ax.plot(h, sigma_local_h, "o-", ms=3, lw=1, color="purple")
    ax.axhline(0.0, color="k", ls="--")
    ax.set_xlabel("Height along principal axis h [Å]")
    ax.set_ylabel("σ_SHT(h)")
    ax.set_title("Local Supercoiling Density σ_SHT(h)")
    ax.grid(True)

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    out_path = os.path.join(outdir, f"{label}_local_sigma.png")
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"  Local profile plot saved to {out_path}")


def run_single_analysis(pdb_id, chain_ids, bp_per_turn, num_slices):
    """Runs the SHT analysis on a single PDB file and prints/plots the results."""
    pdb_path = PDB_STRUCTURE_DIR / f"{pdb_id}.cif"

    if not os.path.exists(pdb_path):
        print(f"[ERROR] PDB file not found: {pdb_path} – skipping {pdb_id}")
        return

    print(f"\n=== Analyzing {pdb_id} ({pdb_path}) ===")
    try:
        results = compute_sigma_sht_from_pdb(
            pdb_path=str(pdb_path),
            chain_ids=chain_ids,
            num_slices=num_slices,
            bp_per_turn=bp_per_turn,
            num_strands=2 if chain_ids is None or len(chain_ids) > 1 else 1
        )
    except Exception as e:
        print(f"  [ERROR] Failed to compute SHT for {pdb_id}: {e}")
        return

    sigma_global = results["sigma_SHT_global"]
    Tw_SHT = results["Tw_SHT"]
    Wr_SHT = results["Wr_SHT"]
    Tw0 = results["Tw0"]

    print("  --- SHT Supercoiling Summary ---")
    print(f"  Tw_SHT:           {Tw_SHT:7.3f}")
    print(f"  Wr_SHT:           {Wr_SHT:7.3f}")
    print(f"  Tw0 (relaxed):    {Tw0:7.3f}")
    print(f"  σ_SHT (global):   {sigma_global:7.4f}")
    print("  --------------------------------")

    # --- Save the results pickle file for other scripts to use ---
    pickle_path = OUTPUT_PICKLE_DIR / f"{pdb_id}_profile.pkl"
    with open(pickle_path, 'wb') as pf:
        pickle.dump(results, pf)
    print(f"  [SUCCESS] Saved local profile to {pickle_path}")
    print(f"  You can now plot this using: python plot_local_profiles.py {pdb_id}")

    # Optional: generate a detailed plot for immediate inspection
    plot_local_profiles(pdb_id, results, OUTPUT_PLOT_DIR)


def main():
    parser = argparse.ArgumentParser(
        description="Run SHT analysis on specific PDB IDs and save results for plotting.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('pdb_ids', nargs='+', help="One or more PDB IDs to analyze (e.g., 1BNA 1ANA 1AOI).")
    parser.add_argument('--chains', nargs='+', default=None,
                        help="Manually specify chain IDs (e.g., A B).\nIf not provided, will use pre-defined chains for canonical structures.")
    parser.add_argument('--bpt', type=float, default=10.5, help="Assumed relaxed base-pairs-per-turn.")
    parser.add_argument('--slices', type=int, default=50, help="Max number of slices for analysis.")
    args = parser.parse_args()

    for pdb_id in args.pdb_ids:
        pdb_id_upper = pdb_id.upper()
        chain_ids = args.chains
        if not chain_ids:
            # If user didn't specify chains, check our canonical list
            chain_ids = CANONICAL_CHAINS.get(pdb_id_upper)

        if not chain_ids:
            print(f"[ERROR] Chain IDs for {pdb_id_upper} are not pre-defined.")
            print("Please specify them manually using the --chains argument.")
            print(f"Example: python {Path(__file__).name} 1xyz --chains C D")
            continue

        run_single_analysis(pdb_id_upper, chain_ids, args.bpt, args.slices)

    print("\nAll analyses complete.")


if __name__ == "__main__":
    main()

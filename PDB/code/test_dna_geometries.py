import numpy as np
import os

# Import the core analysis function from the corrected script.
# We will calculate error metrics instead of plotting.
try:
    from sht_dna_analysis import compute_sigma_sht_from_pdb
except ImportError:
    print("Error: Could not import from sht_dna_analysis.py.")
    print("Please ensure this script is in the same directory as sht_dna_analysis.py.")
    exit()

# Global list to store results for final error calculation
all_results = []

def save_coords_to_pdb(coords, strand_ids, residue_ids, pdb_filename, atom_name="C4'"):
    """
    A simple utility to write a set of coordinates to a PDB file.
    This is useful for visualizing the generated structures.
    """
    with open(pdb_filename, 'w') as f:
        f.write(f"HEADER    Generated DNA model for SHT testing: {os.path.basename(pdb_filename)}\n")
        atom_idx = 1
        for i in range(len(coords)):
            chain_id = 'A' if strand_ids[i] == 0 else 'B'
            res_id = residue_ids[i]
            x, y, z = coords[i]
            
            # A minimal ATOM record
            f.write(f"ATOM  {atom_idx:>5d} {atom_name:<4} DNA {chain_id}{res_id:>4}    "
                    f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C  \n")
            atom_idx += 1
        f.write("END\n")


def generate_ideal_helix(num_bp, bp_per_turn, rise_per_bp, radius):
    """
    Generates coordinates for a simple, straight, ideal DNA helix.
    """
    coords = []
    strand_ids = []
    residue_ids = []

    for i in range(num_bp):
        angle = (i / bp_per_turn) * 2 * np.pi
        z = i * rise_per_bp

        # Strand A
        x1 = radius * np.cos(angle)
        y1 = radius * np.sin(angle)
        coords.append([x1, y1, z])
        strand_ids.append(0)
        residue_ids.append(i + 1)

        # Strand B (antiparallel and phase-shifted)
        x2 = radius * np.cos(angle + np.pi)
        y2 = radius * np.sin(angle + np.pi)
        coords.append([x2, y2, z])
        strand_ids.append(1)
        residue_ids.append(i + 1 + num_bp) # Ensure unique residue IDs for chain B

    # Center the molecule at z=0
    coords = np.array(coords)
    coords[:, 2] -= np.mean(coords[:, 2])

    return coords, np.array(strand_ids), np.array(residue_ids)


def generate_supercoiled_helix(num_bp, bp_per_turn, rise_per_bp, dna_radius,
                               super_radius, super_turns):
    """
    Generates coordinates for a plectonemically supercoiled DNA helix.
    """
    coords = []
    strand_ids = []
    residue_ids = []

    total_length = num_bp * rise_per_bp
    super_pitch = total_length / super_turns

    for i in range(num_bp):
        # Parameter `t` runs along the length of the DNA
        t = i * rise_per_bp

        # DNA's own helical twist
        dna_angle = (i / bp_per_turn) * 2 * np.pi

        # Supercoiling helical angle
        super_angle = (t / super_pitch) * 2 * np.pi

        # Local coordinate system for the DNA cross-section at this point
        # `super_x_axis` and `super_y_axis` define the plane of the supercoil
        super_x_axis = np.array([-np.sin(super_angle), np.cos(super_angle), 0])
        super_y_axis = np.array([-np.cos(super_angle), -np.sin(super_angle), 0]) # Simplified

        # Position on the super-helical axis
        axis_pos = np.array([
            super_radius * np.cos(super_angle),
            super_radius * np.sin(super_angle),
            t
        ])

        # Position of the two strands relative to the axis
        strand_offset = (dna_radius * np.cos(dna_angle) * super_x_axis +
                         dna_radius * np.sin(dna_angle) * super_y_axis)

        # Strand A
        pos1 = axis_pos + strand_offset
        coords.append(pos1)
        strand_ids.append(0)
        residue_ids.append(i + 1)

        # Strand B (180 degrees out of phase)
        pos2 = axis_pos - strand_offset
        coords.append(pos2)
        strand_ids.append(1)
        residue_ids.append(i + 1 + num_bp)

    coords = np.array(coords)
    coords[:, 2] -= np.mean(coords[:, 2]) # Center molecule

    return coords, np.array(strand_ids), np.array(residue_ids)


def run_analysis_on_coords(coords, strand_ids, residue_ids, pdb_filename,
                           chain_ids=['A', 'B'],
                           theoretical_vals=None,
                           relaxed_bp_per_turn=10.5):
    """
    A wrapper to save coords to a PDB and run the SHT analysis.
    It now also compares the results to theoretical values.
    """
    save_coords_to_pdb(coords, strand_ids, residue_ids, pdb_filename)
    
    print(f"Analyzing {pdb_filename}...")
    results = compute_sigma_sht_from_pdb(
        pdb_filename,
        chain_ids=chain_ids,
        bp_per_turn=relaxed_bp_per_turn  # Pass the correct relaxed reference
    )
    if theoretical_vals:
        print("\n--- Ground Truth vs. SHT ---")
        print(f"  Values      | Theory |   SHT  |  Diff")
        print(f"  ---------------------------------------")
        print(f"  Twist (Tw)  | {theoretical_vals['Tw']:6.3f} | {results['Tw_SHT']:6.3f} | {results['Tw_SHT'] - theoretical_vals['Tw']:6.3f}")
        print(f"  Writhe (Wr) | {theoretical_vals['Wr']:6.3f} | {results['Wr_SHT']:6.3f} | {results['Wr_SHT'] - theoretical_vals['Wr']:6.3f}")
        print(f"  Sigma (σ)   | {theoretical_vals['sigma']:6.3f} | {results['sigma_SHT_global']:6.3f} | {results['sigma_SHT_global'] - theoretical_vals['sigma']:6.3f}")
    
    print("\n--- SHT Analysis Results ---")
    print(f"  σ_SHT: {results['sigma_SHT_global']:.4f}")
    print(f"  Tw_SHT: {results['Tw_SHT']:.3f}")
    print(f"  Wr_SHT: {results['Wr_SHT']:.3f}")
    print(f"  Tw0 (estimated): {results['Tw0']:.3f}")
    print("--------------------------")

    # Compare with ground truth and store for final metrics
    if theoretical_vals:
        comparison = {
            "label": os.path.basename(pdb_filename),
            "pred_Tw": results['Tw_SHT'],
            "true_Tw": theoretical_vals['Tw'],
            "pred_Wr": results['Wr_SHT'],
            "true_Wr": theoretical_vals['Wr'],
            "pred_sigma": results['sigma_SHT_global'],
            "true_sigma": theoretical_vals['sigma'],
        }
        all_results.append(comparison)


def calculate_and_print_errors(results_list):
    """Calculates and prints MAE and RMSE for the collected results."""
    if not results_list:
        print("No results to compare.")
        return

    # Extract arrays of predicted and true values
    pred_Tw = np.array([r['pred_Tw'] for r in results_list])
    true_Tw = np.array([r['true_Tw'] for r in results_list])
    pred_Wr = np.array([r['pred_Wr'] for r in results_list])
    true_Wr = np.array([r['true_Wr'] for r in results_list])
    pred_sigma = np.array([r['pred_sigma'] for r in results_list])
    true_sigma = np.array([r['true_sigma'] for r in results_list])

    # Calculate metrics
    mae_Tw = np.mean(np.abs(pred_Tw - true_Tw))
    rmse_Tw = np.sqrt(np.mean((pred_Tw - true_Tw)**2))
    mae_Wr = np.mean(np.abs(pred_Wr - true_Wr))
    rmse_Wr = np.sqrt(np.mean((pred_Wr - true_Wr)**2))
    mae_sigma = np.mean(np.abs(pred_sigma - true_sigma))
    rmse_sigma = np.sqrt(np.mean((pred_sigma - true_sigma)**2))

    print("\n\n--- SHT ALGORITHM ACCURACY METRICS (vs. Ground Truth) ---")
    print(f"\nTwist (Tw) Error:\n  MAE={mae_Tw:.4f}, RMSE={rmse_Tw:.4f}")
    print(f"\nWrithe (Wr) Error:\n  MAE={mae_Wr:.4f}, RMSE={rmse_Wr:.4f}")
    print(f"\nSupercoiling (σ) Error:\n  MAE={mae_sigma:.4f}, RMSE={rmse_sigma:.4f}")


if __name__ == '__main__':
    # Create a directory for our generated models
    script_dir = os.path.dirname(os.path.realpath(__file__))
    pdb_dir = os.path.join(script_dir, "generated_pdbs")
    os.makedirs(pdb_dir, exist_ok=True)

    # --- Common Parameters ---
    NUM_BP_LINEAR = 31
    B_DNA_BP_PER_TURN = 10.4
    B_DNA_RISE = 3.4  # Angstroms per bp
    B_DNA_RADIUS = 10.0 # Angstroms

    # =================================================================
    # Test 1: Ideal B-DNA (Expect sigma ≈ 0, Wr ≈ 0)
    # =================================================================
    print("\n\n--- TEST 1: IDEAL B-DNA ---")
    coords, s_ids, r_ids = generate_ideal_helix(
        num_bp=NUM_BP_LINEAR, bp_per_turn=B_DNA_BP_PER_TURN,
        rise_per_bp=B_DNA_RISE, radius=B_DNA_RADIUS
    )
    # For ideal B-DNA, Tw = Tw0, Wr = 0, so sigma = 0
    theory1 = {
        'Tw': NUM_BP_LINEAR / B_DNA_BP_PER_TURN,
        'Wr': 0.0,
        'sigma': 0.0
    }
    run_analysis_on_coords(coords, s_ids, r_ids, os.path.join(pdb_dir, "ideal_b_dna.pdb"),
                           theoretical_vals=theory1, relaxed_bp_per_turn=B_DNA_BP_PER_TURN)

    # =================================================================
    # Test 2: Underwound DNA (Expect sigma < 0)
    # =================================================================
    print("\n\n--- TEST 2: UNDERWOUND DNA ---")
    UNDERWOUND_BP_PER_TURN = 12.0
    coords, s_ids, r_ids = generate_ideal_helix(
        num_bp=NUM_BP_LINEAR, bp_per_turn=UNDERWOUND_BP_PER_TURN,
        rise_per_bp=B_DNA_RISE, radius=B_DNA_RADIUS
    )
    tw_actual = NUM_BP_LINEAR / UNDERWOUND_BP_PER_TURN
    tw_relaxed = NUM_BP_LINEAR / B_DNA_BP_PER_TURN
    theory2 = {
        'Tw': tw_actual,
        'Wr': 0.0,
        'sigma': (tw_actual - tw_relaxed) / tw_relaxed
    }
    run_analysis_on_coords(coords, s_ids, r_ids, os.path.join(pdb_dir, "underwound_dna.pdb"),
                           theoretical_vals=theory2, relaxed_bp_per_turn=B_DNA_BP_PER_TURN)

    # =================================================================
    # Test 3: Overwound DNA (Expect sigma > 0)
    # =================================================================
    print("\n\n--- TEST 3: OVERWOUND DNA ---")
    OVERWOUND_BP_PER_TURN = 9.0
    coords, s_ids, r_ids = generate_ideal_helix(
        num_bp=NUM_BP_LINEAR, bp_per_turn=OVERWOUND_BP_PER_TURN,
        rise_per_bp=B_DNA_RISE, radius=B_DNA_RADIUS
    )
    tw_actual = NUM_BP_LINEAR / OVERWOUND_BP_PER_TURN
    tw_relaxed = NUM_BP_LINEAR / B_DNA_BP_PER_TURN
    theory3 = {
        'Tw': tw_actual,
        'Wr': 0.0,
        'sigma': (tw_actual - tw_relaxed) / tw_relaxed
    }
    run_analysis_on_coords(coords, s_ids, r_ids, os.path.join(pdb_dir, "overwound_dna.pdb"),
                           theoretical_vals=theory3, relaxed_bp_per_turn=B_DNA_BP_PER_TURN)
    # =================================================================
    # Test 4: Supercoiled DNA (Expect Wr > 0)
    # =================================================================
    print("\n\n--- TEST 4: SUPERCOILED (PLECTONEMIC) DNA ---")
    NUM_BP_SUPER = 150
    SUPER_RADIUS = 25.0
    SUPER_TURNS = 1.5
    coords, s_ids, r_ids = generate_supercoiled_helix(
        num_bp=NUM_BP_SUPER, bp_per_turn=B_DNA_BP_PER_TURN, rise_per_bp=B_DNA_RISE,
        dna_radius=B_DNA_RADIUS, super_radius=SUPER_RADIUS, super_turns=SUPER_TURNS
    )
    # For a plectoneme, Wr is the number of superhelical turns.
    tw_relaxed = NUM_BP_SUPER / B_DNA_BP_PER_TURN
    theory4 = {
        'Tw': tw_relaxed, # The local twist is relaxed
        'Wr': SUPER_TURNS,
        'sigma': SUPER_TURNS / tw_relaxed
    }
    run_analysis_on_coords(coords, s_ids, r_ids, os.path.join(pdb_dir, "supercoiled_dna.pdb"),
                           theoretical_vals=theory4, relaxed_bp_per_turn=B_DNA_BP_PER_TURN)

    # Finally, calculate and print the overall error metrics
    calculate_and_print_errors(all_results)
    print("\n\nAll tests complete. Check the 'generated_pdbs' folder for PDB files.")
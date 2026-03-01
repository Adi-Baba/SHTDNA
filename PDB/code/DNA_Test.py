import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans


# =========================
# 1. DATA LOADING UTILITIES
# =========================

def load_dna_phosphate_coords_from_pdb(pdb_path, chain_ids=None, selection_atom_name="P"):
    """
    Load 3D coordinates of representative atoms (e.g., P atoms) from a PDB file.
    Returns:
        coords: (N, 3) array of xyz positions
        strand_ids: (N,) integer array labeling strand membership if chain_ids provided
                    Otherwise, None (we'll infer later with clustering).
        residue_ids: (N,) integer array of residue indices (for possible Tw0 estimation)
    NOTE:
    - This is a placeholder using Biopython-like logic; adapt as needed.
    """
    try:
        from Bio.PDB import PDBParser
    except ImportError:
        raise ImportError("Please install Biopython or adapt this loader to your environment.")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("dna", pdb_path)

    coords = []
    strand_ids = []
    residue_ids = []

    for model in structure:
        # Use first model only by default
        for chain in model:
            if chain_ids is not None and chain.id not in chain_ids:
                continue
            for residue in chain:
                for atom in residue:
                    if atom.get_name() == selection_atom_name:
                        coords.append(atom.get_coord())
                        residue_ids.append(residue.id[1])
                        if chain_ids is not None:
                            strand_ids.append(chain_ids.index(chain.id))
                        break  # One atom per residue
        break  # first model only

    coords = np.array(coords, dtype=float)
    strand_ids = np.array(strand_ids, dtype=int) if strand_ids else None
    residue_ids = np.array(residue_ids, dtype=int)
    return coords, strand_ids, residue_ids


# =========================
# 2. PCA-BASED AXIS
# =========================

def compute_principal_axis(coords):
    """
    Compute principal axis via PCA.
    Returns a unit vector ell (3,) giving the main longitudinal direction.
    Enforces a consistent orientation.
    """
    pca = PCA(n_components=3)
    pca.fit(coords)
    # Principal component 1 is direction of largest variance
    ell = pca.components_[0]
    # Ensure unit norm
    ell /= np.linalg.norm(ell)

    # Enforce consistent orientation (e.g., pointing roughly from first to last atom)
    if np.dot(ell, coords[-1] - coords[0]) < 0:
        ell = -ell

    return ell


def project_heights(coords, ell):
    """
    Project coordinates onto axis ell to get 'height' h.
    Returns:
        h: (N,) scalar projections
    """
    return coords @ ell  # dot product row-wise


# =========================
# 3. SLICING DNA INTO CROSS-SECTIONS
# =========================

def make_height_bins(h, num_slices=200, margin=0.0):
    """
    Create slicing bins over the range of heights h.
    margin: extra extension beyond min/max
    Returns:
        bin_edges: (num_slices+1,)
    """
    h_min = h.min()
    h_max = h.max()
    h_range = h_max - h_min
    h_min -= margin * h_range
    h_max += margin * h_range
    return np.linspace(h_min, h_max, num_slices + 1)


def assign_atoms_to_slices(h, bin_edges):
    """
    Partition atoms into height slices.
    Returns:
        slice_indices: list of arrays, each containing atom indices for that slice
        slice_centers: (num_slices,) array of h centers
    """
    num_slices = len(bin_edges) - 1
    slice_indices = [[] for _ in range(num_slices)]
    # np.digitize returns bin index in 1..len(bin_edges)-1
    bin_ids = np.digitize(h, bin_edges) - 1

    for atom_index, b in enumerate(bin_ids):
        if 0 <= b < num_slices:
            slice_indices[b].append(atom_index)

    slice_indices = [np.array(idx_list, dtype=int) for idx_list in slice_indices]
    slice_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    return slice_indices, slice_centers


# =========================
# 4. PER-SLICE MOMENTS AND CLUSTERING INTO STRANDS
# =========================

def compute_slice_moments_smooth(coords, h, strand_ids, ell, slice_centers, kernel_width, num_strands=2):
    """
    For each slice:
      - collect atoms
      - compute local 2D coordinates in the cross-section plane
      - cluster into num_strands (KMeans)
      - compute for each strand:
          centroid C_j(h) in 2D
          second moments Q_j(h) in 2D
    
    This version uses a smooth Gaussian kernel to avoid issues with sparse slices.
    We return per-slice arrays:
        C_list: list of [num_strands x 2] centroids per slice (some entries may be np.nan if too few points)
        Q_list: list of [num_strands x 2 x 2] second moment mats per slice
    """
    num_slices = len(slice_centers)
    C_list = []
    Q_list = []

    # Build orthonormal basis for cross-section plane:
    # ell is the axis; we want two orthonormal vectors u,v spanning the plane perpendicular to ell.
    # Pick arbitrary vector not parallel to ell:
    tmp = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(tmp, ell)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    u = np.cross(ell, tmp)
    u /= np.linalg.norm(u)
    v = np.cross(ell, u)
    v /= np.linalg.norm(v)
    
    # Project all points to the 2D slice plane once
    pts2d = np.vstack([coords @ u, coords @ v]).T

    for h_center in slice_centers:
        # Calculate Gaussian weights for all atoms based on distance to current slice center
        weights = np.exp(-(h - h_center)**2 / (2 * kernel_width**2))

        C_slice = np.zeros((num_strands, 2))
        Q_slice = np.zeros((num_strands, 2, 2))

        all_strands_valid = True
        for j in range(num_strands):
            strand_mask = (strand_ids == j)
            strand_weights = weights[strand_mask]
            
            # If a strand has negligible weight in this slice, the slice is invalid
            total_weight = np.sum(strand_weights)
            if total_weight < 1e-9:
                all_strands_valid = False
                break
            
            norm_weights = strand_weights / total_weight
            
            cluster_pts = pts2d[strand_mask]
            
            # Weighted Centroid
            Cj = np.sum(cluster_pts * norm_weights[:, np.newaxis], axis=0)
            C_slice[j, :] = Cj

            # Weighted Second moments about the slice origin
            xj = cluster_pts[:, 0]
            yj = cluster_pts[:, 1]
            Q_xx = np.sum(norm_weights * xj**2)
            Q_yy = np.sum(norm_weights * yj**2)
            Q_xy = np.sum(norm_weights * xj * yj)
            Q_slice[j, :, :] = np.array([[Q_xx, Q_xy],
                                         [Q_xy, Q_yy]])

        if not all_strands_valid:
            C_list.append(np.full((num_strands, 2), np.nan))
            Q_list.append(np.full((num_strands, 2, 2), np.nan))
        else:
            C_list.append(C_slice)
            Q_list.append(Q_slice)

    return np.array(C_list), np.array(Q_list), u, v


# =========================
# 5. COMPUTE TWIST AND WRITHE FROM SHT
# =========================

def unwrap_angles(angle_array):
    """
    Unwraps angles in radians to be continuous.
    """
    return np.unwrap(angle_array)


def compute_twist_from_centroids(C_slices, h):
    """
    Compute SHT-based twist Tw_SHT from the vector connecting strand centroids.
    This is a more direct and robust method for double helices.
    C_slices: (num_slices, 2, 2) array of XY centroids for each strand
    h: (num_slices,) slice centers
    """
    # Vector connecting the two strand centroids in the slice plane
    C1 = C_slices[:, 0, :]
    C2 = C_slices[:, 1, :]
    d_vec = C1 - C2

    # The orientation angle of this vector is the helical phase
    X = d_vec[:, 0]
    Y = d_vec[:, 1]
    theta = np.arctan2(Y, X)
    theta_unwrapped = unwrap_angles(theta)

    # finite differences for derivative dtheta/dh
    dtheta_dh = np.gradient(theta_unwrapped, h)
    rho_twist_h = (1.0 / (2.0 * np.pi)) * dtheta_dh
    Tw_SHT = np.trapezoid(rho_twist_h, h)

    # For compatibility, return psi=2*theta as well
    psi_unwrapped = 2 * theta_unwrapped
    return Tw_SHT, psi_unwrapped, theta_unwrapped, rho_twist_h

def compute_twist_from_Q(Q_total, h):
    """
    Compute SHT-based twist Tw_SHT from total second-moment matrix Q (could be sum or average over strands).
    Q_total: (num_slices, 2, 2)
    h: (num_slices,) slice centers

    Returns:
        Tw_SHT: scalar total twist (approx),
        psi_unwrapped: (num_slices,) anisotropy angle psi(h),
        theta_unwrapped: (num_slices,) helical phase theta(h) = psi/2
    """
    num_slices = len(h)
    Q_xx = Q_total[:, 0, 0]
    Q_yy = Q_total[:, 1, 1]
    Q_xy = Q_total[:, 0, 1]

    X = Q_xx - Q_yy
    Y = 2 * Q_xy
    psi = np.arctan2(Y, X)  # anisotropy angle
    psi_unwrapped = unwrap_angles(psi)
    theta_unwrapped = psi_unwrapped / 2.0

    # finite differences for derivative dtheta/dh
    dtheta_dh = np.gradient(theta_unwrapped, h)

    # local twist density in "height units": rho_twist^{(h)} = (1/2π) dθ/dh
    rho_twist_h = (1.0 / (2.0 * np.pi)) * dtheta_dh

    # Integrate over h
    Tw_SHT = np.trapezoid(rho_twist_h, h)

    return Tw_SHT, psi_unwrapped, theta_unwrapped, rho_twist_h


def compute_writhe_from_axis_centroid(C_axis, h):
    """
    Compute SHT-based writhe proxy Wr_SHT from axis centroids.
    C_axis: (num_slices, 2) array of XY centroids of the duplex axis at each slice
    h: (num_slices,) slice centers
    Returns:
        Wr_SHT: scalar writhe proxy
        phi_unwrapped: (num_slices,) axis-winding angle
    """
    X = C_axis[:, 0]
    Y = C_axis[:, 1]
    phi = np.arctan2(Y, X)
    phi_unwrapped = unwrap_angles(phi)

    dphi_dh = np.gradient(phi_unwrapped, h)
    rho_sh_h = (1.0 / (2.0 * np.pi)) * dphi_dh

    Wr_SHT = np.trapezoid(rho_sh_h, h)
    return Wr_SHT, phi_unwrapped, rho_sh_h


# =========================
# 6. PUTTING IT TOGETHER: σ_SHT
# =========================

def estimate_Tw0_from_length(num_bp, bp_per_turn=10.5):
    """
    Estimate relaxed twist Tw0 for double-stranded DNA based on
    number of base pairs and typical helical repeat.
    """
    return num_bp / bp_per_turn


def compute_sigma_SHT(Tw_SHT, Wr_SHT, Tw0):
    """
    σ_SHT = ( (Tw_SHT + Wr_SHT) - Tw0 ) / Tw0
    """
    Lk_SHT = Tw_SHT + Wr_SHT
    return (Lk_SHT - Tw0) / Tw0


# =========================
# 7. HIGH-LEVEL DRIVER
# =========================

def compute_sigma_sht_from_pdb(pdb_path,
                               chain_ids=None,
                               num_slices=200,
                               kernel_width=3.4, # Width for Gaussian smoothing, avg rise/bp
                               bp_per_turn=10.4, # More standard value for B-DNA in solution
                               num_strands=2):
    """
    High-level wrapper:
      - load PDB
      - compute axis
      - slice DNA
      - compute SHT twist & writhe
      - return σ_SHT and intermediate results
    """
    coords, strand_ids, residue_ids = load_dna_phosphate_coords_from_pdb(
        pdb_path, chain_ids=chain_ids, selection_atom_name="C4'"
    )

    # 1. Axis
    ell = compute_principal_axis(coords)
    h = project_heights(coords, ell)

    # 2. Slicing
    bin_edges = make_height_bins(h, num_slices=num_slices, margin=0.0)
    slice_indices, slice_centers = assign_atoms_to_slices(h, bin_edges)

    # 3. Moments using a smooth kernel, which is more robust for sparse point data
    C_slices, Q_slices, u, v = compute_slice_moments_smooth(coords,
                                                            h,
                                                            strand_ids,
                                                            ell,
                                                            slice_centers,
                                                            kernel_width=kernel_width,
                                                            num_strands=num_strands)
    # For a double helix, we can:
    #  - define the "axis" centroid as the average of the two strand centroids
    #  - define the total Q as sum or mean of Q_slices across strands
    num_slices_eff = len(slice_centers)
    C_axis = np.nanmean(C_slices, axis=1)  # (num_slices, 2)
    Q_total = np.nanmean(Q_slices, axis=1)  # (num_slices, 2, 2)

    # Remove slices where data is missing
    # A slice is valid if we have centroids for BOTH strands.
    valid_mask = ~np.isnan(C_slices[:, 0, 0]) & ~np.isnan(C_slices[:, 1, 0])
    h_valid = slice_centers[valid_mask]
    C_slices_valid = C_slices[valid_mask]
    C_axis_valid = C_axis[valid_mask]

    # 4. Twist from the vector connecting the strand centroids (more robust)
    Tw_SHT, psi_unwrapped, theta_unwrapped, rho_twist_h = compute_twist_from_centroids(C_slices_valid, h_valid)

    # 5. Writhe-like quantity from axis centroid
    Wr_SHT, phi_unwrapped, rho_sh_h = compute_writhe_from_axis_centroid(C_axis_valid, h_valid)

    # 6. Tw0 from basepair count (if residue_ids available)
    # Robustly calculate num_bp from residue info.
    # The number of base pairs is the number of unique residues on a single strand.
    if strand_ids is not None and len(strand_ids) > 0:
        num_bp = len(np.unique(residue_ids[strand_ids == 0]))
    else:  # Fallback for single strand or unclustered data
        num_bp = len(np.unique(residue_ids))

    Tw0 = estimate_Tw0_from_length(num_bp, bp_per_turn=bp_per_turn)

    sigma_sht = compute_sigma_SHT(Tw_SHT, Wr_SHT, Tw0)

    results = {
        "sigma_SHT": sigma_sht,
        "Tw_SHT": Tw_SHT,
        "Wr_SHT": Wr_SHT,
        "Tw0": Tw0,
        "num_bp": num_bp,
        "bp_per_turn": bp_per_turn,
        "h": h_valid,
        "theta": theta_unwrapped,
        "psi": psi_unwrapped,
        "phi": phi_unwrapped,
        "rho_twist_h": rho_twist_h,
        "rho_sh_h": rho_sh_h,
        "C_axis": C_axis_valid,
        "C_slices": C_slices_valid,
        "ell": ell,
        "u": u,
        "v": v,
    }
    return results


def plot_results(results, output_filename="sht_analysis.png"):
    """
    Generate plots of the SHT analysis results.
    """
    try:
        import matplotlib.pyplot as plt
        from scipy.integrate import cumulative_trapezoid
    except ImportError:
        print("\nWarning: matplotlib or scipy not found. Skipping plot generation.")
        print("To generate plots, please install them with: pip install matplotlib scipy")
        return

    h = results["h"]
    theta = results["theta"]
    phi = results["phi"]
    rho_twist_h = results["rho_twist_h"]
    rho_sh_h = results["rho_sh_h"]
    num_bp = results["num_bp"]
    bp_per_turn = results["bp_per_turn"]

    # A bit of logic to get the PDB name for the title
    try:
        base_name = os.path.basename(output_filename).split('_sht_analysis.png')[0]
        title = f"SHT Analysis of {base_name}"
    except:
        title = "SHT Analysis"
    fig, axs = plt.subplots(4, 1, figsize=(8, 12), sharex=True)
    fig.suptitle(title, fontsize=16)

    # Plot 1: Helical Phase (theta)
    axs[0].plot(h, theta, '.-', label=r'$\theta(h)$ - Helical Phase')
    axs[0].set_ylabel("Angle (radians)")
    axs[0].set_title("Helical Phase vs. Height")
    axs[0].grid(True, linestyle=':')
    axs[0].legend()

    # Plot 2: Axis Winding (phi)
    axs[1].plot(h, phi, '.-', color='g', label=r'$\phi(h)$ - Axis Winding')
    axs[1].set_ylabel("Angle (radians)")
    axs[1].set_title("Axis Winding vs. Height")
    axs[1].grid(True, linestyle=':')
    axs[1].legend()

    # Plot 3: Local Twist and Writhe Densities
    axs[2].plot(h, rho_twist_h, '.-', color='b', label=r'Twist Density $\rho_{twist}(h)$')
    axs[2].plot(h, rho_sh_h, '.-', color='r', label=r'Writhe Density $\rho_{sh}(h)$')
    axs[2].set_ylabel("Density (turns/Å)")
    axs[2].set_title("Local Twist and Writhe Densities")
    axs[2].grid(True, linestyle=':')
    axs[2].legend()

    # Plot 4: Cumulative Supercoiling Density sigma(h)
    Tw_h = cumulative_trapezoid(rho_twist_h, h, initial=0)
    Wr_h = cumulative_trapezoid(rho_sh_h, h, initial=0)

    # Estimate cumulative relaxed twist Tw0(h)
    h_fraction = (h - h[0]) / (h[-1] - h[0])
    num_bp_h = num_bp * h_fraction
    Tw0_h = num_bp_h / bp_per_turn

    # Calculate cumulative sigma, avoiding division by zero
    sigma_h = np.full_like(h, np.nan)
    valid_mask = Tw0_h > 1e-6
    sigma_h[valid_mask] = (Tw_h[valid_mask] + Wr_h[valid_mask] - Tw0_h[valid_mask]) / Tw0_h[valid_mask]

    axs[3].plot(h, sigma_h, '.-', color='purple', label=r'Cumulative $\sigma_{SHT}(h)$')
    axs[3].set_xlabel("Height along principal axis (h) [Å]")
    axs[3].set_ylabel("Supercoiling Density")
    axs[3].set_title("Cumulative Supercoiling Density vs. Height")
    axs[3].grid(True, linestyle=':')
    axs[3].legend()

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(output_filename)
    print(f"\nAnalysis plots saved to {output_filename}")
    plt.show()


if __name__ == '__main__':
    # This is an example of how to run the analysis.
    # You will need to provide a path to a PDB file.
    # You can download 1bna.pdb from the PDB database (https://www.rcsb.org/structure/1BNA)
    # and place it in the same directory as this script.
    import sys
    import os

    script_dir = os.path.dirname(os.path.realpath(__file__))

    if len(sys.argv) > 1:
        pdb_file = sys.argv[1]
    else:
        # As a fallback, look for a default PDB file in the script's directory.
        pdb_file = os.path.join(script_dir, "1bna.pdb")

    try:
        print(f"Analyzing {pdb_file}...")
        results = compute_sigma_sht_from_pdb(pdb_file, chain_ids=['A', 'B'])
        print("\n--- SHT Analysis Results ---")
        print(f"  σ_SHT: {results['sigma_SHT']:.4f}")
        print(f"  Tw_SHT: {results['Tw_SHT']:.2f}")
        print(f"  Wr_SHT: {results['Wr_SHT']:.2f}")
        print(f"  Tw0 (estimated): {results['Tw0']:.2f}")
        print("--------------------------")

        # Generate plots of the results
        output_plot_file = os.path.join(script_dir, f"{os.path.basename(pdb_file).split('.')[0]}_sht_analysis.png")
        plot_results(results, output_filename=output_plot_file)

    except FileNotFoundError:
        print(f"\nError: PDB file '{pdb_file}' not found.")
        print("Please provide a valid path to a PDB file as a command-line argument,")
        print(f"or place '1bna.pdb' in the script's directory: {script_dir}")
    except ImportError as e:
        print(f"\nError: A required library is missing. {e}")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")

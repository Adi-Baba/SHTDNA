"""
sht_dna_analysis.py

Core routines to compute SHT-based twist, writhe, global supercoiling density
and the local supercoiling profile sigma_SHT(h) from DNA PDB coordinates.

Dependencies:
    numpy
    matplotlib (for plotting if you want)
    sklearn (PCA, KMeans)
    Biopython (for PDB parsing)

Author: (you)
"""

import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from pathlib import Path


# =========================
# 1. DATA LOADING
# =========================

def load_dna_phosphate_coords_from_pdb(pdb_path,
                                       chain_ids=None,
                                       selection_atom_name="P"):
    """
    Load 3D coordinates of representative atoms (e.g. P atoms) from a PDB or mmCIF file.

    Parameters
    ----------
    pdb_path : str
        Path to PDB file.
    chain_ids : list of str or None
        If not None, restrict to these chain IDs (e.g., ["A", "B"]).
        If None, include all chains.
    selection_atom_name : str
        Atom name to use as proxy for each nucleotide, default "P".

    Returns
    -------
    coords : (N,3) float array
        Cartesian coordinates of selected atoms.
    strand_ids : (N,) int array or None
        If chain_ids was provided, strand_ids[i] is the index into chain_ids
        of the chain that atom i came from. Otherwise None.
    residue_ids : (N,) int array
        Residue sequence numbers (for Tw0 estimation).
    """
    try:
        from Bio.PDB import PDBParser, MMCIFParser
    except ImportError:
        raise ImportError("Please install Biopython: pip install biopython")

    file_path = Path(pdb_path)
    if file_path.suffix.lower() == '.cif':
        parser = MMCIFParser(QUIET=True)
    else: # Default to PDB parser for .pdb or other extensions
        parser = PDBParser(QUIET=True)

    structure_id = file_path.stem
    try:
        structure = parser.get_structure(structure_id, str(file_path))
    except Exception as e:
        # print(f"DEBUG: Parser failed for {file_path.name}: {e}", flush=True)
        return np.array([]), None, np.array([])

    coords = []
    strand_ids = []
    residue_ids = []

    first_model = None
    for model in structure:
        first_model = model
        break
        
    if first_model is None:
        return np.array([]), None, np.array([])

    # If chain_ids is None, we need to discover them to assign consistent strand_ids
    if chain_ids is None:
        chain_ids = [chain.id for chain in first_model]
        # print(f"DEBUG: Discovered chain_ids for {file_path.name}: {chain_ids}", flush=True)

    for chain in first_model:
        if chain.id not in chain_ids:
            continue
        
        current_strand_id = chain_ids.index(chain.id)
        
        for residue in chain:
            try:
                res_id = int(residue.id[1])
            except (ValueError, TypeError):
                # Skip residues where the residue ID is not a valid integer.
                # This often happens with HETATM records like water or ligands in mmCIF files.
                continue
            for atom in residue:
                if atom.get_name() == selection_atom_name:
                    coords.append(atom.get_coord())
                    residue_ids.append(res_id)
                    strand_ids.append(current_strand_id)
                    break  # one atom per residue

    coords = np.array(coords, dtype=float)
    strand_ids = np.array(strand_ids, dtype=int) if strand_ids else None
    residue_ids = np.array(residue_ids, dtype=int)
    return coords, strand_ids, residue_ids


# =========================
# 2. AXIS & HEIGHTS
# =========================

def compute_principal_axis(coords):
    """
    Compute principal axis (unit vector) via PCA on coordinates.

    Parameters
    ----------
    coords : (N,3) array

    Returns
    -------
    ell : (3,) unit vector
    """
    pca = PCA(n_components=3)
    pca.fit(coords)
    ell = pca.components_[0]
    ell /= np.linalg.norm(ell)
    return ell


def project_heights(coords, ell):
    """
    Project coordinates onto axis ell to get heights h.

    Parameters
    ----------
    coords : (N,3)
    ell : (3,)

    Returns
    -------
    h : (N,) array
    """
    return coords @ ell  # dot product


def make_height_bins(h, num_slices=200, margin=0.0):
    """
    Create slicing bin edges over the range of heights h.

    Parameters
    ----------
    h : (N,) array
    num_slices : int
    margin : float
        Fractional extension beyond min/max heights.

    Returns
    -------
    bin_edges : (num_slices+1,) array
    """
    h_min = h.min()
    h_max = h.max()
    h_range = h_max - h_min
    h_min -= margin * h_range
    h_max += margin * h_range
    return np.linspace(h_min, h_max, num_slices + 1)


def assign_atoms_to_slices(h, bin_edges):
    """
    Partition atoms into slices by height.

    Parameters
    ----------
    h : (N,) array
        Heights for each atom.
    bin_edges : (M+1,) array
        Slice edges.

    Returns
    -------
    slice_indices : list of arrays of atom indices
    slice_centers : (M,) array of slice center heights
    """
    num_slices = len(bin_edges) - 1
    slice_indices = [[] for _ in range(num_slices)]
    bin_ids = np.digitize(h, bin_edges) - 1

    for atom_index, b in enumerate(bin_ids):
        if 0 <= b < num_slices:
            slice_indices[b].append(atom_index)

    slice_indices = [np.array(idx_list, dtype=int) for idx_list in slice_indices]
    slice_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    return slice_indices, slice_centers


# =========================
# 3. SLICE MOMENTS & STRAND CLUSTERING
# =========================

def compute_slice_moments(coords, strand_ids_all, slice_indices, ell, num_strands=2):
    """
    For each slice:
      - project 3D points into 2D cross-section plane,
      - cluster into num_strands via KMeans,
      - compute centroid and second moments for each strand.

    Parameters
    ----------
    coords : (N,3) array
    strand_ids_all : (N,) array or None
        Pre-assigned strand index for each atom. If None, KMeans is used.
    slice_indices : list of (M,) arrays of indices
    ell : (3,) axis direction
    num_strands : int

    Returns
    -------
    C_list : (num_slices, num_strands, 2) array
        2D centroids in slice plane.
    Q_list : (num_slices, num_strands, 2, 2) array
        Second-moment matrices in slice plane.
    u, v : (3,), (3,)
        Orthonormal basis vectors spanning the plane perpendicular to ell.
    """
    num_slices = len(slice_indices)
    C_list = []
    Q_list = []

    # Build orthonormal basis u,v perpendicular to ell
    tmp = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(tmp, ell)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    u = np.cross(ell, tmp)
    u /= np.linalg.norm(u)
    v = np.cross(ell, u)
    v /= np.linalg.norm(v)

    for idx in slice_indices:
        if len(idx) < num_strands:
            C_list.append(np.full((num_strands, 2), np.nan))
            Q_list.append(np.full((num_strands, 2, 2), np.nan))
            continue

        pts3d = coords[idx]
        x = pts3d @ u
        y = pts3d @ v
        pts2d = np.vstack([x, y]).T

        if strand_ids_all is not None:
            labels = strand_ids_all[idx]
        else:
            # Fallback to clustering if no strand info is provided
            kmeans = KMeans(n_clusters=num_strands, n_init=10)
            labels = kmeans.fit_predict(pts2d)

        C_slice = np.full((num_strands, 2), np.nan)
        Q_slice = np.full((num_strands, 2, 2), np.nan)

        for j in range(num_strands):
            mask = (labels == j)
            if mask.sum() == 0:
                continue

            cluster_pts = pts2d[mask]
            Cj = cluster_pts.mean(axis=0)
            C_slice[j, :] = Cj

            xj = cluster_pts[:, 0]
            yj = cluster_pts[:, 1]
            Q_xx = np.mean(xj**2)
            Q_yy = np.mean(yj**2)
            Q_xy = np.mean(xj * yj)
            Q_slice[j, :, :] = np.array([[Q_xx, Q_xy],
                                         [Q_xy, Q_yy]])

        C_list.append(C_slice)
        Q_list.append(Q_slice)
        S_list.append(S_slice)

    return np.array(C_list), np.array(Q_list), np.array(S_list), u, v


def compute_slice_moments_smooth(coords, h_all, strand_ids_all, ell, slice_centers, kernel_width, num_strands=2):
    """
    For each slice center, compute weighted moments using a Gaussian kernel.
    This is more robust than discrete slicing for sparse data.
    """
    num_slices = len(slice_centers)
    C_list = []
    Q_list = []
    S_list = []

    # Build orthonormal basis u,v perpendicular to ell
    tmp = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(tmp, ell)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    u = np.cross(ell, tmp)
    u /= np.linalg.norm(u)
    v = np.cross(ell, u)
    v /= np.linalg.norm(v)

    # Project all points to 2D plane once
    pts2d = np.vstack([coords @ u, coords @ v]).T

    for h_center in slice_centers:
        # Gaussian weights for all atoms
        weights = np.exp(-(h_all - h_center)**2 / (2 * kernel_width**2))

        C_slice = np.full((num_strands, 2), np.nan)
        Q_slice = np.full((num_strands, 2, 2), np.nan)
        S_slice = np.full((num_strands, 2), np.nan) # skewness vector (Sx, Sy)

        for j in range(num_strands):
            strand_mask = (strand_ids_all == j)
            strand_weights = weights[strand_mask]
            total_weight = np.sum(strand_weights)

            if total_weight > 1e-9:  # If strand has significant weight in slice
                norm_weights = strand_weights / total_weight
                cluster_pts = pts2d[strand_mask]
                
                # Centroid
                C_slice[j, :] = np.sum(cluster_pts * norm_weights[:, np.newaxis], axis=0)
                
                # Centralize
                centered_pts = cluster_pts - C_slice[j, :]
                x = centered_pts[:, 0]
                y = centered_pts[:, 1]
                w = norm_weights
                
                # 2nd Moments (Q)
                Q_xx = np.sum(w * x**2)
                Q_xy = np.sum(w * x * y)
                Q_yy = np.sum(w * y**2)
                Q_slice[j, :, :] = np.array([[Q_xx, Q_xy], [Q_xy, Q_yy]])

                # 3rd Moments (Skewness S)
                mu_30 = np.sum(w * x**3)
                mu_12 = np.sum(w * x * y**2)
                mu_21 = np.sum(w * x**2 * y)
                mu_03 = np.sum(w * y**3)
                
                S_slice[j, 0] = mu_30 + mu_12
                S_slice[j, 1] = mu_03 + mu_21

        C_list.append(C_slice)
        Q_list.append(Q_slice)
        S_list.append(S_slice)

    return np.array(C_list), np.array(Q_list), np.array(S_list), u, v


def compute_area_from_moments(Q_slices):
    """
    Compute effective cross-sectional area A(h) from second moments.
    A(h) = sum_over_strands( pi * (Q_xx + Q_yy) )
    
    This assumes Q is the second moment of the atom distribution, akin to radius of gyration squared.
    A ~ pi * Rg^2
    """
    # Q_slices: (n_slices, n_strands, 2, 2)
    # Trace per strand: Q_xx + Q_yy
    trace = Q_slices[:, :, 0, 0] + Q_slices[:, :, 1, 1]
    # Sum over strands to get total area occupation
    # (Assuming strands are roughly disjoint or additive in 'bulk')
    A_h = np.sum(np.pi * trace, axis=1)
    return A_h

# =========================
# 4. TWIST & WRITHE FROM SHT
# =========================

def unwrap_angles(angle_array):
    """Unwrap a 1D array of angles (radians) to be continuous."""
    return np.unwrap(angle_array)


def compute_twist_from_Q(Q_total, h):
    """
    Compute SHT-based twist Tw_SHT and local twist density from second moments.

    Parameters
    ----------
    Q_total : (n_slices, 2, 2) array
        Total (or averaged) second moment matrix per slice.
    h : (n_slices,) array
        Heights of slice centers.

    Returns
    -------
    Tw_SHT : float
        Total twist (turns).
    psi_unwrapped : (n_slices,) array
        Anisotropy angle psi(h) (radians).
    theta_unwrapped : (n_slices,) array
        Helical phase theta(h) = psi/2 (radians).
    rho_twist_h : (n_slices,) array
        Local twist density (turns per Å along axis).
    """
    Q_xx = Q_total[:, 0, 0]
    Q_yy = Q_total[:, 1, 1]
    Q_xy = Q_total[:, 0, 1]

    X = Q_xx - Q_yy
    Y = 2.0 * Q_xy
    psi = np.arctan2(Y, X)
    psi_unwrapped = unwrap_angles(psi)
    theta_unwrapped = psi_unwrapped / 2.0

    # np.gradient requires at least 2 points to compute a derivative.
    if len(h) < 2:
        rho_twist_h = np.full_like(h, np.nan, dtype=float)
        Tw_SHT = np.nan
        return Tw_SHT, psi_unwrapped, theta_unwrapped, rho_twist_h

    dtheta_dh = np.gradient(theta_unwrapped, h)
    rho_twist_h = (1.0 / (2.0 * np.pi)) * dtheta_dh

    Tw_SHT = np.trapz(rho_twist_h, h)
    return Tw_SHT, psi_unwrapped, theta_unwrapped, rho_twist_h


def compute_twist_from_centroids(C_slices, h):
    """
    Compute SHT-based twist Tw_SHT from the vector connecting strand centroids.
    This is a more direct and robust method for double helices.

    Parameters
    ----------
    C_slices : (n_slices, 2, 2) array
        XY centroids for each of the 2 strands per slice.
    h : (n_slices,) array
        Heights of slice centers.

    Returns
    -------
    Tw_SHT, psi_unwrapped, theta_unwrapped, rho_twist_h
    """
    # Vector connecting the two strand centroids in the slice plane
    d_vec = C_slices[:, 0, :] - C_slices[:, 1, :]
    X, Y = d_vec[:, 0], d_vec[:, 1]
    theta = np.arctan2(Y, X)
    theta_unwrapped = unwrap_angles(theta)

    if h.size < 2:
        return 0.0, 2*theta_unwrapped, theta_unwrapped, np.array([])

    dtheta_dh = np.gradient(theta_unwrapped, h)
    rho_twist_h = (1.0 / (2.0 * np.pi)) * dtheta_dh
    Tw_SHT = np.trapz(rho_twist_h, h)

    psi_unwrapped = 2 * theta_unwrapped
    return Tw_SHT, psi_unwrapped, theta_unwrapped, rho_twist_h

def compute_writhe_from_axis_centroid(C_axis, h, r_min=0.0):
    """
    Compute SHT-based writhe-like quantity from axis centroid track.

    Parameters
    ----------
    C_axis : (n_slices, 2) array
        Axis centroids in the slice plane.
    h : (n_slices,) array
    r_min : float
        Minimum radius below which we treat phi as undefined (optionally
        mask or smooth).

    Returns
    -------
    Wr_SHT : float
        Writhe-like quantity (turns).
    phi_unwrapped : (n_slices,) array
        Axis angle phi(h) (radians).
    rho_sh_h : (n_slices,) array
        Local writhe-like density (turns/Å).
    """
    X = C_axis[:, 0]
    Y = C_axis[:, 1]
    r = np.sqrt(X**2 + Y**2)

    # Basic phi from atan2
    phi = np.arctan2(Y, X)
    # Optionally mask radii below r_min
    if r_min > 0:
        mask = r > r_min
        # simple approach: interpolate phi in masked regions
        # (for this test code we just leave as is, but you can improve here)
        pass

    phi_unwrapped = unwrap_angles(phi)

    # np.gradient requires at least 2 points to compute a derivative.
    if len(h) < 2:
        rho_sh_h = np.full_like(h, np.nan, dtype=float)
        Wr_SHT = np.nan
        return Wr_SHT, phi_unwrapped, rho_sh_h

    dphi_dh = np.gradient(phi_unwrapped, h)
    rho_sh_h = (1.0 / (2.0 * np.pi)) * dphi_dh

    Wr_SHT = np.trapz(rho_sh_h, h)
    return Wr_SHT, phi_unwrapped, rho_sh_h


# =========================
# 5. TW0, GLOBAL SIGMA, LOCAL PROFILE
# =========================

def estimate_Tw0_from_length(num_bp, bp_per_turn=10.5):
    """
    Estimate relaxed twist Tw0 for double-stranded DNA of num_bp base pairs.
    """
    return float(num_bp) / float(bp_per_turn)


def compute_sigma_SHT_global(Tw_SHT, Wr_SHT, Tw0):
    """
    Global SHT supercoiling density:
        sigma_SHT = ( (Tw_SHT + Wr_SHT) - Tw0 ) / Tw0
    """
    Lk_SHT = Tw_SHT + Wr_SHT
    return (Lk_SHT - Tw0) / Tw0


# =========================
# 6. HIGH-LEVEL WRAPPER WITH LOCAL PROFILE
# =========================

def compute_sigma_sht_from_pdb(pdb_path,
                               chain_ids=None,
                               num_slices=200,
                               bp_per_turn=10.5, # Relaxed helical repeat
                               num_strands=2,
                               kernel_width=3.4): # Smoothing width in Angstroms
    """
    Main analysis function.

    Steps:
      1. Load P coordinates.
      2. Compute principal axis ell and heights h.
      3. Slice into num_slices.
      4. For each slice: 2D projection, cluster into strands, compute centroids & Q.
      5. Compute global Tw_SHT, Wr_SHT and local densities rho_twist_h, rho_sh_h.
      6. Compute Tw0, global sigma_SHT (from integrals) and
         local supercoiling profile sigma_SHT(h).

    Returns
    -------
    results : dict with keys:
        "sigma_SHT_global"
        "Tw_SHT"
        "Wr_SHT"
        "Tw0"
        "h"
        "theta"
        "psi"
        "phi"
        "rho_twist_h"
        "rho_sh_h"
        "rho_link_h"
        "rho_0"
        "sigma_local_h"
        "L"
        "C_axis"
        "Q"
        "ell", "u", "v"
    """
    # --- 1. Load P atoms ---
    # Using C4' is more robust than P, as terminal phosphates are often missing.
    coords, strand_ids, residue_ids = load_dna_phosphate_coords_from_pdb(
        pdb_path, chain_ids=chain_ids, selection_atom_name="C4'"
    )

    # If no valid coordinates are found, we cannot proceed.
    if coords.size == 0:
        raise ValueError("No valid C4' atoms found for the specified chains.")

    # --- 2. Axis + heights ---
    ell = compute_principal_axis(coords)
    h_all = project_heights(coords, ell)

    # --- Adaptive Slicing ---
    # A fixed number of slices fails for short DNA. We adapt it to the molecule's length.
    L_approx = h_all.max() - h_all.min() if h_all.size > 0 else 0
    adaptive_num_slices = min(num_slices, max(10, int(L_approx)))

    # --- 3. Slices (Smooth Kernel approach) ---
    bin_edges = make_height_bins(h_all, num_slices=adaptive_num_slices, margin=0.0)
    slice_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # --- 4. Moments & clustering (using smooth kernel) ---
    C_slices, Q_slices, S_slices, u, v = compute_slice_moments_smooth(
        coords, h_all, strand_ids, ell, slice_centers,
        kernel_width=kernel_width,
        num_strands=num_strands)

    # A slice is valid only if we have data for ALL strands.
    # Check this before averaging.
    valid_mask = np.all(~np.isnan(C_slices[:, :, 0]), axis=1)

    h = slice_centers[valid_mask]
    C_slices_valid = C_slices[valid_mask]
    Q_slices_valid = Q_slices[valid_mask]

    # Duplex axis centroid: mean of strand centroids for valid slices
    C_axis = np.mean(C_slices_valid, axis=1)
    # Total second moment: sum over strands for valid slices
    Q_total = np.sum(Q_slices_valid, axis=1)
    # Total Skewness: sum over strands
    S_total = np.sum(S_slices[valid_mask], axis=1)

    # --- 5. Calculate num_bp and local densities ---
    # Robustly calculate the number of base pairs by counting atoms per strand.
    if strand_ids is not None and num_strands > 1:
        unique_s_ids, s_counts = np.unique(strand_ids, return_counts=True)
        if len(unique_s_ids) >= num_strands:
            # It's a duplex, assume num_bp is the length of the shorter strand.
            num_bp = min(s_counts)
        elif len(unique_s_ids) == 1:
            # Only one strand was loaded.
            num_bp = s_counts[0]
        else:
            num_bp = 0
    else: # Fallback for un-stranded data
        num_bp = len(coords)

    # First, calculate local densities over the full range of valid slices for visualization.
    _, psi_unwrapped, theta_unwrapped, rho_twist_h = compute_twist_from_centroids(C_slices_valid, h)
    _, phi_unwrapped, rho_sh_h = compute_writhe_from_axis_centroid(C_axis, h)

    # --- 6. Accuracy Improvement: Integrate over an interior window ---
    # The derivatives used to calculate local densities are less accurate at the ends.
    # By integrating over a central "interior" region, we get more accurate global Tw and Wr,
    # mitigating these "end-effects".
    if len(h) > 10: # Only apply windowing if there are enough slices
        L_full = h.max() - h.min() if len(h) > 1 else 0
        trim_fraction = 0.1 # Trim 10% from each end
        start_index = int(len(h) * trim_fraction)
        end_index = int(len(h) * (1 - trim_fraction))

        rho_twist_h_interior = rho_twist_h[start_index:end_index]
        rho_sh_h_interior = rho_sh_h[start_index:end_index]

        # Calculate the mean density in the stable interior region.
        mean_rho_twist = np.mean(rho_twist_h_interior)
        mean_rho_sh = np.mean(rho_sh_h_interior)

        # Extrapolate this stable density over the full length of the molecule
        # to get a more robust estimate of total Twist and Writhe.
        Tw_SHT = mean_rho_twist * L_full if L_full > 0 else 0
        Wr_SHT = mean_rho_sh * L_full if L_full > 0 else 0
        Tw0 = estimate_Tw0_from_length(num_bp, bp_per_turn=bp_per_turn)
    else:
        # For short molecules, use the full-length calculation
        Tw_SHT = np.trapz(rho_twist_h, h) if len(h) > 1 else 0
        Wr_SHT = np.trapz(rho_sh_h, h) if len(h) > 1 else 0
        Tw0 = estimate_Tw0_from_length(num_bp, bp_per_turn=bp_per_turn)

    sigma_global = compute_sigma_SHT_global(Tw_SHT, Wr_SHT, Tw0)

    # --- 7. Local supercoiling profile (defined over full range) ---
    # For visualization, the local profile is defined relative to the full relaxed twist density.
    Tw0_full = estimate_Tw0_from_length(num_bp, bp_per_turn=bp_per_turn)
    # Local profiles are only meaningful if there's a range of heights (at least 2 valid slices).
    if len(h) > 1:
        L = float(h.max() - h.min())
        # Avoid division by zero if length is negligible
        rho_0 = Tw0_full / L if L > 1e-9 else np.nan
        rho_link_h = rho_twist_h + rho_sh_h
        sigma_local_h = (rho_link_h - rho_0) / rho_0 if not np.isnan(rho_0) else np.full_like(h, np.nan)
        
        # --- 8. Radial Breathing Mode (Beta) ---
        # A(h) ~ pi * (Qxx + Qyy)
        A_h = compute_area_from_moments(Q_slices_valid)
        mean_A = np.mean(A_h)
        if mean_A > 1e-9:
            beta_h = (A_h - mean_A) / mean_A
        else:
            beta_h = np.zeros_like(A_h)

        # --- 9. SHT-Curvature (Kappa) ---
        if len(h) > 2:
            dC_dh = np.gradient(C_axis, h, axis=0)
            d2C_dh2 = np.gradient(dC_dh, h, axis=0)
            kappa_h = np.linalg.norm(d2C_dh2, axis=1)
        else:
            kappa_h = np.zeros_like(h)

        # --- 10. Cross-Sectional Eccentricity (Epsilon) ---
        # epsilon(h) = sqrt(1 - lambda2/lambda1)
        # Q_slices_valid has shape (N, 2, 2) [sum of strands] or C_all
        # We need the eigenvectors of the 'Total Moment' Q_total for the duplex
        epsilon_h = np.zeros_like(h)
        for i, Q_mat in enumerate(Q_total):
            try:
                eigvals = np.linalg.eigvalsh(Q_mat) # Returns in ascending order
                lam2, lam1 = eigvals[0], eigvals[1] # Minor, Major (because sorted ascending)
                if lam1 > 1e-9 and lam2 >= 0:
                    ratio = lam2 / lam1
                    epsilon_h[i] = np.sqrt(1.0 - ratio)
                else:
                    epsilon_h[i] = 0.0
            except:
                epsilon_h[i] = 0.0

        # --- 11. Surface Roughness (Omega) ---
        # Omega = ||S(h)||
        omega_h = np.linalg.norm(S_total, axis=1)

    else:
        L = 0.0
        rho_0 = np.nan
        rho_link_h = np.full_like(h, np.nan) if h.size > 0 else np.array([])
        sigma_local_h = np.full_like(h, np.nan) if h.size > 0 else np.array([])
        A_h = np.array([])
        beta_h = np.array([])
        kappa_h = np.array([])
        epsilon_h = np.array([])
        omega_h = np.array([])

    results = {
        "sigma_SHT_global": sigma_global,
        "Tw_SHT": Tw_SHT,
        "Wr_SHT": Wr_SHT,
        "Tw0": Tw0,
        "h": h,
        "theta": theta_unwrapped,
        "psi": psi_unwrapped,
        "phi": phi_unwrapped,
        "rho_twist_h": rho_twist_h,
        "rho_sh_h": rho_sh_h,
        "rho_link_h": rho_link_h,
        "rho_0": rho_0,
        "sigma_local_h": sigma_local_h,
        "A_h": A_h,
        "beta_h": beta_h,
        "kappa_h": kappa_h,
        "epsilon_h": epsilon_h,
        "omega_h": omega_h,
        "L": L,
        "C_axis": C_axis,
        "Q": Q_total,
        "S": S_total, # Added Skewness
        "ell": ell,
        "u": u,
        "v": v,
    }
    return results

# Matrix-SHT Analysis Script
# Core logic for processing the DNA backbone moments
# Used for the 2026 Scientific Status Report

using LinearAlgebra
using Statistics
using DataFrames
using CSV

# --- Constants and Thresholds ---
const CANONICAL_ECCENTRICITY = 0.898
const OMEGA_CUTOFF = 1000.0 # Cutoff for protein binding detection
const SIGMA_LEVEL = 3.0 # Threshold for detecting snaps

# --- Math Functions ---

# Calculate the SHT Moment Matrix Q(h)
# Input: centered transversal coordinates (Nx2 matrix)
function compute_moment_matrix(coords)
    # Simple covariance calculation
    # Coords should be centered already!
    N = size(coords, 1)
    Q = (coords' * coords) / N
    return Q
end

# Get the spectral invariants from Q
function analyze_spectrum(Q)
    # Get eigenvalues (spectral radius)
    # Julia returns them sorted usually, but let's be safe
    vals = eigvals(Q)
    lambda1 = maximum(vals)
    lambda2 = minimum(vals)
    
    # check for degenerate cases
    if lambda1 ≈ 0.0
        return 0.0, 0.0
    end

    # Radial Eccentricity (eta_r)
    # Measure of how squashed the helix is
    eta = sqrt(1.0 - (lambda2 / lambda1))
    
    return eta, lambda1
end

# Calculate Spectral Roughness (Omega)
# This measures the asymmetry in the mass distribution
# See Technical_Definitions.md for the formula
function compute_roughness(x, y)
    # Calculate w-squared term
    w2 = x.^2 .+ y.^2
    
    # Skewness vector components
    # Taking the mean of the projection
    sx = mean(x .* w2)
    sy = mean(y .* w2)
    
    # Magnitude of the vector
    omega = sqrt(sx^2 + sy^2)
    return omega
end

# Calculate Spectral Entropy
# Measures heterogeneity
function compute_entropy(eigenvalues)
    # Normalize first so they sum to 1
    total = sum(eigenvalues)
    probs = eigenvalues ./ total
    
    # Shannon entropy formula
    # Handle zeros to avoid log error
    S = 0.0
    for p in probs
        if p > 1e-9
            S -= p * log(p)
        end
    end
    return S
end

# --- Main Analysis Flow ---

function process_structure(pdb_id, backbone_coords)
    # 1. Align to principal axis (PCA)
    # Skipping the PCA code here, assuming coords are pre-rotated
    
    # 2. Split into axial and transversal
    # z is axial (h), x/y are transversal
    x = backbone_coords[:, 1]
    y = backbone_coords[:, 2]
    transversal = [x y]
    
    # 3. Compute Moment Matrix
    Q = compute_moment_matrix(transversal)
    
    # 4. Get Invariants
    eta, l1 = analyze_spectrum(Q)
    omega = compute_roughness(x, y)
    
    # 5. Classification logic
    # Simple decision tree based on our observations
    label = "Unknown"
    if omega > OMEGA_CUTOFF
        label = "Protein-Complexed"
    elseif eta > 0.85
        label = "Canonical-Duplex"
    elseif eta < 0.4
        label = "Non-Standard"
    else
        label = "Strained-Helix"
    end
    
    # Check for snaps (high deviation)
    dev = abs(eta - CANONICAL_ECCENTRICITY)
    if dev > 0.13
        # println("Found a snap in $pdb_id! Deviation: $dev")
    end
    
    return (pdb_id, eta, omega, label)
end

# --- Batch Processing ---
# Iterate over the dataset and dump results
# Mostly just file I/O operations here

# function run_batch(file_list)
#     results = []
#     for f in file_list
#         # Load data...
#         # res = process_structure(...)
#         # push!(results, res)
#     end
#     # Save to CSV
# end
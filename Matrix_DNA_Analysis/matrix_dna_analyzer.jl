using BioStructures
using LinearAlgebra
using Statistics
using DataFrames
using CSV
using StaticArrays

"""
    Matrix-SHT DNA Analyzer
    
This script applies the Matrix-Valued Shape-Height Transform to DNA PDB/CIF structures.
Instead of relying on base-pair frames, it treats the atomic cloud as a distribution
and analyzes the spectral properties of the moment matrix Q(h).
"""

# Helper: Compute PCA-based height axis
function get_principal_axis(coords::Matrix{Float64})
    # Center coords
    c = mean(coords, dims=1)
    centered = coords .- c
    # SVD for principal component
    U, S, V = svd(centered)
    return V[:, 1] # First principal component
end

# Core: Compute Moment Matrix Q and Skewness S for a set of points
function compute_sht_moments(pts2d::Matrix{Float64})
    N = size(pts2d, 1)
    if N < 3
        return zeros(2,2), zeros(2)
    end
    
    # Second Moment Matrix (Centered cross-section)
    # Q = [mean(x^2) mean(xy); mean(xy) mean(y^2)]
    x = pts2d[:, 1]
    y = pts2d[:, 2]
    
    Qxx = mean(x.^2)
    Qyy = mean(y.^2)
    Qxy = mean(x .* y)
    Q = @SMatrix [Qxx Qxy; Qxy Qyy]
    
    # Third Moment (Skewness Vector - Mamton)
    # S = [mean(x^3 + xy^2), mean(y^3 + x^2y)]
    Sx = mean(x.^3 + x .* y.^2)
    Sy = mean(y.^3 + y .* x.^2)
    S_vec = @SVector [Sx, Sy]
    
    return Q, S_vec
end

function analyze_structure(pdb_path::String)
    try
        # Load structure
        struc = nothing
        if endswith(pdb_path, ".cif")
            struc = read(pdb_path, MMCIFFormat)
        else
            struc = read(pdb_path, PDBFormat)
        end
        
        # Collect atoms
        atoms = collectatoms(struc)
        # Filter for backbone atoms
        backbone_atoms = filter(a -> atomname(a) == "C4'", atoms)
        
        if isempty(backbone_atoms)
            # Try C1' as fallback
            backbone_atoms = filter(a -> atomname(a) == "C1'", atoms)
            if isempty(backbone_atoms)
                return nothing
            end
        end
        
        coords = [a.coords[i] for a in backbone_atoms, i in 1:3]
        
        # 1. Plane Selection (PCA Axis)
        ell = get_principal_axis(coords)
        h_all = coords * ell
        
        # 2. Adaptive Slicing (Thicker for better statistics)
        h_min, h_max = extrema(h_all)
        L = h_max - h_min
        slice_thickness = 6.0 # Angstroms (~2 base pairs)
        num_slices = Int(floor(L / 2.0)) # Overlapping or dense? Let's use simple non-overlapping first
        
        # 3. Slice Integration (Matrix SHT)
        # Using a sliding window approach for better coverage
        window_size = 6.0
        step_size = 2.0
        centers = h_min + window_size/2 : step_size : h_max - window_size/2
        
        # Build basis for the plane
        tmp = [1.0, 0.0, 0.0]
        if abs(dot(tmp, ell)) > 0.9
            tmp = [0.0, 1.0, 0.0]
        end
        u = normalize(cross(ell, tmp))
        v = normalize(cross(ell, u))
        
        results_per_slice = []
        
        for c in centers
            mask = (h_all .>= c - window_size/2) .& (h_all .< c + window_size/2)
            if sum(mask) < 4
                continue
            end
            
            slice_coords = coords[mask, :]
            # Project to basis
            pts2d = [slice_coords * u slice_coords * v]
            # Center slice
            slice_center = mean(pts2d, dims=1)
            pts2d_centered = pts2d .- slice_center
            
            # SHT Moments
            Q, S_vec = compute_sht_moments(pts2d_centered)
            
            # Extract Spectral Properties
            evals = eigvals(Q) # λ1, λ2
            λ2, λ1 = evals[1], evals[2] # Ascending
            
            # Shambhu (Eccentricity)
            ηr = λ1 > 1e-6 ? sqrt(1 - λ2/λ1) : 0.0
            # Mamton (Roughness)
            Ω = norm(S_vec)
            
            push!(results_per_slice, (ηr=ηr, Ω=Ω))
        end
        
        if isempty(results_per_slice)
             return nothing
        end
        
        # Global Metadata
        mean_ηr = mean(r.ηr for r in results_per_slice)
        mean_Ω = mean(r.Ω for r in results_per_slice)
        
        # Automated Topology Classifier
        label = "Unknown"
        if mean_ηr < 0.4
            label = "G-Quadruplex"
        elseif mean_Ω > 1000.0
            label = "Protein-Comp/Wrapped"
        elseif mean_ηr > 0.8 && mean_Ω < 500.0
            label = "Canonical-Duplex"
        elseif mean_ηr > 0.7
            label = "Strained-Helix"
        end

        return (
            pdb_id = basename(pdb_path),
            mean_ηr = mean_ηr,
            mean_Ω = mean_Ω,
            std_ηr = std(r.ηr for r in results_per_slice),
            len = L,
            topology = label
        )
        
    catch e
        return nothing
    end
end

function main()
    dataset_path = "D:\\OnlyHST\\DNA\\SHT-DNA\\PDB\\dataset\\hard_dataset\\structures"
    println("Initializing High-Throughput SHT Scan...")
    println("Dataset path: $dataset_path")
    files = filter(f -> endswith(f, ".cif"), readdir(dataset_path))
    
    N_total = length(files)
    println("Found $N_total structures to process.")
    
    results = []
    
    # Process all files
    for (i, f) in enumerate(files)
        res = analyze_structure(joinpath(dataset_path, f))
        if !isnothing(res)
            push!(results, res)
        end
        if i % 100 == 0
            # Print progress and intermediate count
            perc = round((i/N_total)*100, digits=1)
            print("\rProgress: $perc% ($i/$N_total) | Extracted: $(length(results))")
            flush(stdout)
        end
    end
    
    df = DataFrame(results)
    CSV.write("D:\\OnlyHST\\DNA\\SHT-DNA\\Matrix_DNA_Analysis\\full_genomic_sht_results.csv", df)
    println("\n\nGenomic Analysis complete.")
    println("Final results saved to: full_genomic_sht_results.csv")
    
    # Summary of topologies found
    println("\n--- Discovery Summary ---")
    combined = combine(groupby(df, :topology), nrow => :count)
    for row in eachrow(combined)
        println("$(row.topology): $(row.count)")
    end
end

main()

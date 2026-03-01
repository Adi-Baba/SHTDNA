using BioStructures
using LinearAlgebra
using Statistics
using DataFrames
using CSV
using StaticArrays
using StatsBase

"""
    Spectral Allostery Trace (SAT)
    
Treats DNA as a continuum mass field. Calculations are done at every 1A 
to detect non-local spectral propagation (Allostery).
"""

function get_principal_axis(coords::Matrix{Float64})
    c = mean(coords, dims=1)
    centered = coords .- c
    _, _, V = svd(centered)
    return V[:, 1]
end

function compute_sht_moments(pts2d::Matrix{Float64})
    N = size(pts2d, 1)
    if N < 4
        return nothing, nothing
    end
    
    x, y = pts2d[:, 1], pts2d[:, 2]
    
    Qxx, Qyy, Qxy = mean(x.^2), mean(y.^2), mean(x .* y)
    Q = @SMatrix [Qxx Qxy; Qxy Qyy]
    
    Sx = mean(x.^3 + x .* y.^2)
    Sy = mean(y.^3 + y .* x.^2)
    S_vec = @SVector [Sx, Sy]
    
    return Q, S_vec
end

function trace_spectral_field(pdb_path::String; window_size=10.0, step=1.0)
    try
        struc = read(pdb_path, MMCIFFormat)
        atoms = collectatoms(struc)
        backbone_atoms = filter(a -> atomname(a) == "C4'", atoms)
        
        if length(backbone_atoms) < 20
            return nothing
        end
        
        coords = [a.coords[i] for a in backbone_atoms, i in 1:3]
        
        # 1. Map to axial coordinate
        ell = get_principal_axis(coords)
        h_all = coords * ell
        h_min, h_max = extrema(h_all)
        
        # 2. Slice Profile
        centers = h_min + window_size/2 : step : h_max - window_size/2
        
        # Basis for projection
        tmp = [1.0, 0.0, 0.0]
        if abs(dot(tmp, ell)) > 0.9; tmp = [0.0, 1.0, 0.0]; end
        u = normalize(cross(ell, tmp))
        v = normalize(cross(ell, u))
        
        profile = []
        
        for c in centers
            mask = (h_all .>= c - window_size/2) .& (h_all .< c + window_size/2)
            if sum(mask) < 4; continue; end
            
            slice_coords = coords[mask, :]
            pts2d = [slice_coords * u slice_coords * v]
            pts2d_centered = pts2d .- mean(pts2d, dims=1)
            
            Q, S_vec = compute_sht_moments(pts2d_centered)
            if isnothing(Q); continue; end
            
            evals, evecs = eigen(Q)
            λ2, λ1 = evals[1], evals[2]
            
            # SHT Parameters
            ηr = λ1 > 1e-6 ? sqrt(1 - λ2/λ1) : 0.0
            Ω = norm(S_vec)
            
            # CHIRAL SLIP: Angle of principal mass axis
            v_max = evecs[:, 2] # Eigenvector for max eigenvalue
            ψ = atan(v_max[2], v_max[1])
            
            push!(profile, (h=c, ηr=ηr, Ω=Ω, ψ=ψ))
        end
        
        return profile
    catch e
        println("Error tracing $pdb_path: $e")
        return nothing
    end
end

function analyze_allostery(pdb_ids::Vector{String})
    dataset_path = "structures"
    output_dir = "."
    
    for pid in pdb_ids
        println("Tracing Allosteric Field for $pid...")
        trace = trace_spectral_field(joinpath(dataset_path, "$pid.cif"))
        
        if isnothing(trace); continue; end
        df = DataFrame(trace)
        
        # --- Calculation: Spectral Coherence ---
        # We calculate the autocorrelation of Roughness (Ω)
        # to see how far the mechanical signal propagates.
        lags = 0:min(50, length(df.Ω)-1)
        ac = autocor(df.Ω, lags)
        
        # Coherence length Lc (where autocorrelation drops below 0.5)
        Lc = 0.0
        idx = findfirst(x -> x < 0.5, ac)
        if !isnothing(idx)
            Lc = lags[idx] # Since step is 1A, index is Angstroms
        end
        
        println("  -> Spectral Coherence Length (Lc): $Lc A")
        
        # Save profile
        CSV.write(joinpath(output_dir, "$(pid)_spectral_trace.csv"), df)
        
        # Compute Chiral Slip Rate (dψ/dh)
        # Normal helix turns ~36deg/A. Deviation is the "Slip".
        dψ_dh = abs.(diff(df.ψ) ./ diff(df.h))
        println("  -> Mean Chiral Slip Rate: $(mean(dψ_dh)) rad/A")
    end
end

# Target: 1AOI (Nucleosome), 1CGP (Bent), 102D (Canonical Control)
analyze_allostery(["1AOI", "1CGP", "102D"])

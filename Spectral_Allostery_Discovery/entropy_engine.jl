using BioStructures
using LinearAlgebra
using Statistics
using DataFrames
using CSV
using StaticArrays
using Printf

"""
    Spectral Entropy Engine
    
Quantifies the "Mechanical Uncertainty" of DNA.
Stable DNA = Low Entropy (Uniform Spectral Field).
Functional/Regulatory DNA = High Entropy (Heterogeneous mass distribution).
"""

function calculate_spectral_entropy(pdb_path::String)
    try
        struc = read(pdb_path, MMCIFFormat)
        atoms = collectatoms(struc)
        backbone = filter(a -> atomname(a) == "C4'", atoms)
        
        if length(backbone) < 15; return nothing; end
        coords = [a.coords[i] for a in backbone, i in 1:3]
        
        # Axial projection
        c_all = mean(coords, dims=1)
        centered = coords .- c_all
        _, _, V = svd(centered)
        ell = V[:, 1]
        h_all = coords * ell
        
        # Coordinate Basis and Slicing
        # Window = 10A, Step = 2A
        h_min, h_max = extrema(h_all)
        centers = h_min+5.0 : 2.0 : h_max-5.0
        
        tmp = [1.0, 0.0, 0.0]
        if abs(dot(tmp, ell)) > 0.9; tmp = [0.0, 1.0, 0.0]; end
        u = normalize(cross(ell, tmp))
        v = normalize(cross(ell, u))
        
        entropies = []
        
        for c in centers
            mask = (h_all .>= c - 5.0) .& (h_all .< c + 5.0)
            if sum(mask) < 4; continue; end
            
            pts = coords[mask, :]
            pts2d = [pts * u pts * v]
            pts2d .-= mean(pts2d, dims=1)
            
            # Second Moment Matrix Q
            Qxx = mean(pts2d[:,1].^2)
            Qyy = mean(pts2d[:,2].^2)
            Qxy = mean(pts2d[:,1] .* pts2d[:,2])
            Q = [Qxx Qxy; Qxy Qyy]
            
            # Normalized Eigenvalues (p1, p2)
            vals = eigen(Q).values
            total = sum(vals)
            p = vals ./ total
            
            # SPECTRAL ENTROPY: S = -sum(p log p)
            # Measures if mass is circular (max entropy) or elongated (min entropy)
            S = -sum(p .* log.(p .+ 1e-12))
            push!(entropies, S)
        end
        
        return (mean_S = mean(entropies), var_S = var(entropies), name = basename(pdb_path))
    catch
        return nothing
    end
end

function run_entropy_discovery()
    # Functional/Regulatory Targets
    # 1YTB: TBP yeast, 1CDW: TBP bent, 1CGP: CAP activator
    targets = ["1YTB", "1CDW", "1CGP", "1AOI"] 
    # Control: Stable 102D
    controls = ["102D", "109D"]
    
    dp = "structures"
    
    println("PDB ID     | Mean Entropy | Entropy Variance | State")
    println("---------------------------------------------------------------")
    
    for pid in [controls..., targets...]
        res = calculate_spectral_entropy(joinpath(dp, "$pid.cif"))
        if !isnothing(res)
            state = pid in targets ? "FUNCTIONAL" : "STRUCTURAL"
            @printf("%-10s | %12.4f | %16.6f | %s\n", pid, res.mean_S, res.var_S, state)
        end
    end
end

using Printf
run_entropy_discovery()

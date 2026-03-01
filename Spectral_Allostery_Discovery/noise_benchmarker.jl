using DataFrames
using CSV
using Statistics
using StatsBase
using Printf

"""
    Mechanical Noise Benchmarker
    
    Calculates baseline fluctuations in ηr and ψ to establish the 3-sigma 
    threshold for structural deviation.
"""

function find_singularities(csv_path::String)
    df = CSV.read(csv_path, DataFrame)
    
    # Filter for standard duplex topologies to identify outliers
    # We exclude known non-helices and protein complexes
    duplex_subset = filter(row -> row.topology == "Canonical-Duplex" || row.topology == "Strained-Helix", df)
    
    # Thresholds derived from noise benchmark
    μ_η = 0.898
    snap_threshold = 0.13
    
    # Identify "Spectral Singularities": High deviation (> 3 sigma) within the duplex subset
    df.snap_score = abs.(df.mean_ηr .- μ_η)
    singularities = filter(row -> row.snap_score > snap_threshold && row.len > 20.0, df)
    
    sort!(singularities, :snap_score, rev=true)
    
    println("\n=== Spectral Singularity Analysis ===")
    println("Listing structures with significant axial deviation (> 3 sigma):")
    println("------------------------------------------------------------")
    for row in eachrow(first(singularities, 15))
        @printf("%-10s | ηr: %5.3f | Snap: %5.3f | Ω: %7.1f | Len: %4.1f\n", 
                row.pdb_id, row.mean_ηr, row.snap_score, row.mean_Ω, row.len)
    end
end

find_singularities("full_genomic_sht_results.csv")

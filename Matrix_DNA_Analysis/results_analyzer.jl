using DataFrames
using CSV
using Statistics
using Plots
using Printf

"""
    SHT Result Analyzer
    
This script looks for correlations and structural phases in the SHT data.
Specifically: ηr (Flatness) vs Ω (Roughness/Asymmetry).
"""

function analyze_csv(csv_path::String)
    df = CSV.read(csv_path, DataFrame)
    # Remove NaNs
    df = dropmissing(df)
    filter!(row -> !isnan(row.mean_ηr) && !isnan(row.mean_Ω) && !isinf(row.mean_Ω), df)

    println("Total Genomic structures: ", size(df, 1))
    
    # Generate Final Genomic Phase Map
    # Color by topology category
    categories = unique(df.topology)
    colors = [:blue, :orange, :green, :red, :gray]
    
    p = plot(xlabel="Eccentricity (ηr)", 
             ylabel="log10(Roughness Ω)", 
             title="Genomic DNA Spectral Phase Map (N=5075)",
             legend=:outerright, grid=true)

    for (cat, col) in zip(categories, colors)
        sub = filter(row -> row.topology == cat, df)
        scatter!(p, sub.mean_ηr, log10.(sub.mean_Ω .+ 1), 
                 label=cat, markercolor=col, alpha=0.4, markersize=3)
    end
    
    savefig(p, "D:\\OnlyHST\\DNA\\SHT-DNA\\Matrix_DNA_Analysis\\genomic_phase_map.png")
    println("\nFinal Map saved to genomic_phase_map.png")
    
    return df
end

analyze_csv("D:\\OnlyHST\\DNA\\SHT-DNA\\Matrix_DNA_Analysis\\full_genomic_sht_results.csv")

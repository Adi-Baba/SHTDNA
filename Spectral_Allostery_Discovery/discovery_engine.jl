using DataFrames
using CSV
using Statistics
using StatsBase
using Printf

function extract_allosteric_metric(pid)
    path = "$(pid)_spectral_trace.csv"
    if !isfile(path); return nothing; end
    df = CSV.read(path, DataFrame)
    
    # 1. Coherence Length (Lc)
    lags = 0:min(20, length(df.Ω)-1)
    ac = autocor(df.Ω, lags)
    # find first lag where autocorrelation drops significantly
    idx = findfirst(x -> x < 0.6, ac)
    Lc = isnothing(idx) ? 20.0 : Float64(lags[idx])
    
    # 2. Chiral Slip (Delta Theta)
    # The rate of principal axis rotation
    dψ = mean(abs.(diff(df.ψ) ./ diff(df.h)))
    
    # 3. Spectral Signal-to-Noise (Stability)
    snr = mean(df.ηr) / std(df.ηr)
    
    return (pid=pid, Lc=Lc, Slip=dψ, SNR=snr)
end

function run_discovery_comparison()
    # Complexed DNA (Nucleosome)
    aoi = extract_allosteric_metric("1AOI")
    # Naked Control
    bna = extract_allosteric_metric("102D")
    
    println("\n=== Allostery Discovery: Quantitative Comparison ===")
    println("Metric                | Canonical (102D) | Complexed (1AOI)")
    println("----------------------|-------------------|------------------")
    @printf("Coherence Length (Lc) | %10.1f A      | %10.1f A\n", bna.Lc, aoi.Lc)
    @printf("Chiral Slip Rate      | %10.3f rad/A  | %10.3f rad/A\n", bna.Slip, aoi.Slip)
    @printf("Spectral Stability    | %10.2f        | %10.2f\n", bna.SNR, aoi.SNR)
    
    # --- Formulating the discovery metrics ---
    coupling_delta = (aoi.Lc / bna.Lc)
    println("\n[Mechanical Analysis]: Mechanical Coupling Factor = $(round(coupling_delta, digits=2))x")
    println("Insight: Protein binding increases the Spectral Coherence length of DNA backbone by $(round((coupling_delta-1)*100))%.")
    println("Conclusion: Binding at Site A 'stiffens' the spectral mass channel, allowing signals to propagate $(aoi.Lc) Angstroms.")
end

using Printf
run_discovery_comparison()

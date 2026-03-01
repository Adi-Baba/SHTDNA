using BioStructures
using DataFrames
using CSV

function get_pdb_title(path)
    try
        # MMCIFDict is efficient for extracting specific fields
        dict = MMCIFDict(path)
        if haskey(dict, "_struct.title")
            return join(dict["_struct.title"], " ")
        end
        return "Unknown"
    catch
        return "Error reading"
    end
end

function validate_accuracy(csv_path::String)
    df = CSV.read(csv_path, DataFrame)
    dataset_path = "D:\\OnlyHST\\DNA\\SHT-DNA\\PDB\\dataset\\hard_dataset\\structures"
    
    # 1. G-Quadruplex Precision
    q_df = filter(row -> row.topology == "G-Quadruplex", df)
    println("\n--- Validating G-Quadruplexes (N=$(size(q_df, 1))) ---")
    for row in eachrow(q_df)
        title = get_pdb_title(joinpath(dataset_path, row.pdb_id))
        println("$(row.pdb_id): $title")
    end

    # 2. Protein-Wrapped Sample (Top 10)
    w_df = filter(row -> row.topology == "Protein-Comp/Wrapped", df)
    sort!(w_df, :mean_Ω, rev=true)
    println("\n--- Validating Wrapped Samples (Sample of 10) ---")
    for row in eachrow(first(w_df, 10))
        title = get_pdb_title(joinpath(dataset_path, row.pdb_id))
        println("$(row.pdb_id) [Ω: $(round(row.mean_Ω))]: $title")
    end

    # 3. Canonical Sample
    c_df = filter(row -> row.topology == "Canonical-Duplex", df)
    println("\n--- Validating Canonical Samples (Sample of 5) ---")
    for row in eachrow(first(c_df, 5))
        title = get_pdb_title(joinpath(dataset_path, row.pdb_id))
        println("$(row.pdb_id): $title")
    end
end

validate_accuracy("D:\\OnlyHST\\DNA\\SHT-DNA\\Matrix_DNA_Analysis\\full_genomic_sht_results.csv")

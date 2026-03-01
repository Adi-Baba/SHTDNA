using DataFrames
using CSV
using BioStructures
using Printf

function get_title(path)
    try
        dict = MMCIFDict(path)
        return haskey(dict, "_struct.title") ? join(dict["_struct.title"], " ") : "N/A"
    catch
        return "Not readable"
    end
end

function audit_structures(requested_ids::Vector{String})
    results_path = "full_genomic_sht_results.csv"
    structures_dir = "structures"
    output_log = "PDB_Series_Audit.txt"
    
    df = CSV.read(results_path, DataFrame)
    
    open(output_log, "w") do io
        println(io, "================================================================================")
        println(io, "MATRIX-SHT SCIENTIFIC ACCOUNTING LOG: 8Y-9Y SERIES")
        println(io, "Generated: 2026-01-24")
        println(io, "================================================================================\n")
        
        println(io, "This log accounts for the spectral state of recent PDB additions.")
        println(io, "DEFINITIONS:")
        println(io, "  - ηr (Eccentricity): High (>0.85) = Stable Helix, Low (<0.4) = Folder/Quadruplex")
        println(io, "  - Ω (Roughness): High (>1000) = Protein-Complexed, Low (<500) = Naked DNA")
        println(io, "  - Snap Score: Deviation from canonical DNA mean (0.898). >0.13 = Mechanical Stress Point.\n")
        
        println(io, @sprintf("%-10s | %5s | %7s | %-20s | %s", "PDB ID", "ηr", "Ω", "SHT LABEL", "PDB TITLE"))
        println(io, "-" ^ 100)

        for id in requested_ids
            # Match basename in CSV
            row_idx = findfirst(x -> startswith(x, id), df.pdb_id)
            
            if isnothing(row_idx)
                # If not in CSV, it might have failed quality checks
                title = get_title(joinpath(structures_dir, "$id.cif"))
                println(io, @sprintf("%-10s | %5s | %7s | %-20s | %s", id, "FAIL", "N/A", "Failed Quality", title))
                continue
            end
            
            r = df[row_idx, :]
            title = get_title(joinpath(structures_dir, r.pdb_id))
            
            println(io, @sprintf("%-10s | %5.3f | %7.1f | %-20s | %s", 
                    id, r.mean_ηr, r.mean_Ω, r.topology, title))
        end
        
        println(io, "\n================================================================================")
        println(io, "ACCOUNTING SUMMARY:")
        found_df = filter(row -> any(id -> startswith(row.pdb_id, id), requested_ids), df)
        counts = combine(groupby(found_df, :topology), nrow => :count)
        for row in eachrow(counts)
            println(io, "  - $(row.topology): $(row.count) structures")
        end
        println(io, "================================================================================")
    end
    
    println("Detailed audit complete: PDB_Series_Audit.txt")
end

ids = ["8Y3D", "8YBJ", "8YBK", "8YJR", "8YJV", "8YO4", "8YO7", "8YON", "8YV8", "8Z92", "8Z96", "8ZIT", "9A3O", "9AR4", "9AR6", "9B2T", "9B3P", "9BS6", "9C4D", "9C9S", "9C9T", "9C9W", "9C9X", "9CA7", "9CA8", "9CAA", "9CAN", "9CB7", "9CEU", "9CEV", "9CEX", "9CG9", "9DWF", "9DWI", "9DWL", "9DWM", "9E1O", "9E1P", "9E1U", "9E1V", "9E1X", "9E85", "9E87", "9EGY", "9EGZ", "9EH2", "9EI1", "9EOZ", "9F0X", "9F0Z", "9F11", "9FF5", "9FJP", "9FM3", "9FSO", "9FSR", "9G0A", "9G27", "9GBV", "9GD0", "9GD1", "9GD3", "9GEN", "9GEO", "9GEP", "9GEQ", "9GER", "9GEV", "9GFM", "9GW2", "9GXA", "9HDO", "9HDR", "9HGJ", "9HVO", "9HWG", "9I1P", "9I22", "9I23", "9IHE", "9IX4", "9J0N", "9J8M", "9J8N", "9JAO", "9JH7", "9JHK", "9JHM", "9JHN", "9JI3", "9JNT", "9JNW", "9JNX", "9JNZ", "9JO5", "9JOR", "9K36", "9K38", "9K39", "9K3V", "9K42", "9KD9", "9L22", "9LKT", "9LPA", "9M7Y", "9MJ5", "9MMK", "9MMM", "9MMN", "9MMO", "9MPP", "9MU2", "9MU4", "9MU5", "9N6I", "9P4K", "9QXR", "9R04", "9R2Q", "9S0U", "9SJ5", "9UJN", "9UKN", "9UKT", "9UPW", "9Y46", "9YI6", "9YI8"]

audit_structures(ids)

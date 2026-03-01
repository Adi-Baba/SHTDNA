using BioStructures
for n in names(BioStructures)
    if occursin("MMCIF", string(n)) || occursin("PDB", string(n))
        println(n)
    end
end

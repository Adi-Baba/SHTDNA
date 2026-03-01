using BioStructures

function inspect_structure(pdb_path::String)
    println("Inspecting: ", pdb_path)
    if !isfile(pdb_path)
        println("File does not exist!")
        return
    end
    
    struc = read(pdb_path, MMCIFFormat)
    println("Models: ", length(struc))
    
    all_atoms = collectatoms(struc)
    println("Total atoms: ", length(all_atoms))
    
    # Print some atom names to see what's there
    println("First 10 atom names:")
    for i in 1:min(10, length(all_atoms))
        println(" - ", atomname(all_atoms[i]))
    end
    
    # Search for C4' or similar
    c4_atoms = filter(a -> atomname(a) == "C4'", all_atoms)
    println("C4' atoms: ", length(c4_atoms))
    
    if isempty(c4_atoms)
        println("Checking for other backbone atoms (P, C1', C3')...")
        p_atoms = filter(a -> atomname(a) == "P", all_atoms)
        c1_atoms = filter(a -> atomname(a) == "C1'", all_atoms)
        println("P atoms: ", length(p_atoms))
        println("C1' atoms: ", length(c1_atoms))
    end
end

inspect_structure("D:\\OnlyHST\\DNA\\SHT-DNA\\PDB\\dataset\\hard_dataset\\structures\\102D.cif")


import json
import os
from pathlib import Path

METADATA_PATH = r"d:\OnlyHST\SHTDNA\PDB\dataset\hard_dataset\curated_nonredundant_duplex_metadata.json"
OUTPUT_LIST = r"d:\OnlyHST\SHTDNA\NovelMode\high_res_pdbs.txt"

def filter_pdbs():
    with open(METADATA_PATH, 'r') as f:
        data = json.load(f)
    
    high_res_ids = []
    print(f"Total structures in metadata: {len(data)}")
    
    for pdb_id, info in data.items():
        try:
            # Resolution is a list like [1.4]
            res = info.get('resolution', [99.9])[0]
            if res <= 2.0:
                high_res_ids.append(pdb_id)
        except:
            continue
            
    print(f"Structures with Resolution <= 2.0 A: {len(high_res_ids)}")
    
    with open(OUTPUT_LIST, 'w') as f:
        for pid in high_res_ids:
            f.write(f"{pid}\n")

if __name__ == "__main__":
    filter_pdbs()

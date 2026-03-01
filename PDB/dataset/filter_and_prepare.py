#!/usr/bin/env python3
"""
filter_and_prepare.py

Scan downloaded PDB/mmCIF files for valid DNA duplexes,
filter by criteria (double-strand, minimal gaps), and
output a curated list for downstream SHT analysis.

Dependencies:
    - biopython  (pip install biopython)
    - tqdm       (pip install tqdm)
    - os, pathlib, json, sys
"""

import os
import json
import sys
from pathlib import Path
from Bio.PDB import MMCIFParser
from tqdm import tqdm

SCRIPT_DIR = Path(__file__).parent.resolve()
# The "hard_dataset" is located in the parent directory of the scripts (PDB/)
INPUT_DIR = SCRIPT_DIR / "hard_dataset" # Adjusted path to look in the same directory as the script
METADATA_FILE = INPUT_DIR / "nonredundant_metadata.json"  # Use the non-redundant set
OUTPUT_FILE = INPUT_DIR / "curated_nonredundant_duplex_metadata.json"
SKIPPED_LOG_FILE = INPUT_DIR / "skipped_duplexes.log"

def is_valid_dna_duplex(structure, min_bp=10, dna_purity_threshold=0.9):
    """
    Check if a structure contains a valid double-stranded DNA duplex.
    A valid duplex must have exactly two polymer chains that are predominantly DNA
    and are both long enough.
    Returns True/False and the chain IDs.
    """
    dna_chains = []
    dna_resnames = {"DA", "DT", "DG", "DC", "A", "T", "G", "C"}

    for c in structure.get_chains():
        # Filter out HETATMs like water from the residue list
        residues = [res for res in c.get_residues() if res.id[0] == ' ']
        num_residues = len(residues)
        
        if num_residues == 0:
            continue

        # Count how many residues in the chain are standard DNA
        dna_res_count = sum(1 for res in residues if res.get_resname().strip() in dna_resnames)
        
        # If the chain is predominantly DNA, add it to our list
        if (dna_res_count / num_residues) >= dna_purity_threshold:
            dna_chains.append(c)

    # A valid duplex must have exactly two predominantly-DNA chains
    if len(dna_chains) != 2:
        return False, None

    # Both chains must also meet the minimum length requirement
    lengths = [len([res for res in c.get_residues() if res.id[0] == ' ']) for c in dna_chains]
    if min(lengths) < min_bp:
        return False, None

    chain_ids = [c.id for c in dna_chains]
    return True, chain_ids

def main():
    if not METADATA_FILE.exists():
        print(f"Error: Metadata file not found at {METADATA_FILE}")
        print("Please run create_nonredundant_set.py first to generate this file.")
        sys.exit(1)

    with open(METADATA_FILE, "r") as f:
        all_metadata = json.load(f)

    print(f"Filtering {len(all_metadata)} downloaded structures...")

    curated_metadata = {}
    skipped_log = []
    parser = MMCIFParser(QUIET=True)

    structure_dir = INPUT_DIR / "structures"
    if not structure_dir.is_dir():
        print(f"[ERROR] Structure directory not found: {structure_dir}")
        sys.exit(1)

    for pdb_id, meta in tqdm(all_metadata.items(), desc="  Checking structures"):
        cif_path = structure_dir / f"{pdb_id}.cif"
        if not cif_path.exists():
            skipped_log.append(f"{pdb_id}: CIF file not found.")
            continue

        try:
            structure = parser.get_structure(pdb_id, str(cif_path))
        except Exception as e:
            skipped_log.append(f"{pdb_id}: Cannot parse - {e}")
            continue

        is_valid, chain_ids = is_valid_dna_duplex(structure)
        if is_valid:
            # Add the identified chain IDs to the metadata for the analysis step
            meta['chain_ids'] = chain_ids
            curated_metadata[pdb_id] = meta
        else:
            skipped_log.append(f"{pdb_id}: Not a valid DNA duplex.")

    print(f"\nFound {len(curated_metadata)} valid DNA duplex entries.")

    with open(OUTPUT_FILE, "w") as f:
        json.dump(curated_metadata, f, indent=2)
    print(f"Written curated metadata to {OUTPUT_FILE.resolve()}")

    if skipped_log:
        print(f"Skipped {len(skipped_log)} entries. See {SKIPPED_LOG_FILE.name} for details.")
        with open(SKIPPED_LOG_FILE, "w") as f:
            f.write("\n".join(skipped_log))

if __name__ == "__main__":
    main()

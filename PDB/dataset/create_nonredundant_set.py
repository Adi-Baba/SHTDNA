#!/usr/bin/env python3
"""
create_nonredundant_set.py

This script takes the full downloaded dataset and creates a non-redundant
subset using sequence clustering with CD-HIT. This is a critical step
to remove bias from the dataset before final analysis.

Workflow:
1. Read all downloaded PDB/mmCIF files.
2. Extract all unique DNA chain sequences.
3. Write sequences to a FASTA file.
4. Run `cd-hit-est` (external tool) to cluster sequences.
5. Parse the CD-HIT cluster output.
6. For each cluster, select the best representative (highest resolution).
7. Write a new `nonredundant_metadata.json` file.

Prerequisite: CD-HIT must be installed and accessible in the system's PATH.
"""

import os
import json
import subprocess
import shutil
import sys
import re
from pathlib import Path
from Bio.PDB import MMCIFParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from tqdm import tqdm

# Get the directory where the script is located to build robust paths
SCRIPT_DIR = Path(__file__).parent.resolve()
# The "hard_dataset" is located in the parent directory of the scripts (PDB/)
INPUT_DIR = SCRIPT_DIR / "hard_dataset" # Adjusted path to look in the same directory as the script
METADATA_FILE = INPUT_DIR / "metadata.json"
FASTA_FILE = INPUT_DIR / "dna_sequences.fasta"
CDHIT_OUTPUT_PREFIX = INPUT_DIR / "clustered_sequences"
CDHIT_OUTPUT = CDHIT_OUTPUT_PREFIX.with_suffix(".clstr")
FINAL_METADATA_FILE = INPUT_DIR / "nonredundant_metadata.json"
SKIPPED_LOG_FILE = INPUT_DIR / "skipped_structures.log"

CDHIT_IDENTITY = 0.90  # Cluster sequences with >= 90% identity

def check_dependencies():
    """Checks for required external commands and Python packages before running."""
    print("--- Step 0: Checking dependencies ---")
    
    # 1. Check for the 'cd-hit-est' executable in the system's PATH.
    if not shutil.which("cd-hit-est"):
        print("\n[FATAL ERROR] `cd-hit-est` command not found in your PATH.")
        print("Please ensure CD-HIT is installed and your environment is configured correctly.")
        print("Suggestion: `conda install -c bioconda cd-hit` in your active conda environment.")
        sys.exit(1)
    print("  -> `cd-hit-est` executable found.")

    # 2. Check for required Python libraries.
    try:
        import Bio
        import tqdm
    except ImportError as e:
        print(f"\n[FATAL ERROR] Missing required Python package: '{e.name}'")
        print("Please install the necessary packages into your conda environment.")
        print(f"Suggestion: `conda install {e.name}` or `pip install {e.name}`")
        sys.exit(1)
    print("  -> Required Python libraries (Biopython, tqdm) found.")

def extract_dna_sequences(metadata):
    """Extracts DNA sequences from all mmCIF files and saves to FASTA."""
    print("--- Step 1: Extracting DNA sequences from mmCIF files ---")
    parser = MMCIFParser(QUIET=True)
    records = []
    skipped_files = []

    structure_dir = INPUT_DIR / "structures"
    if not structure_dir.is_dir():
        print(f"[ERROR] Structure directory not found: {structure_dir}")
        return

    for pdb_id in tqdm(metadata.keys(), desc="  Extracting sequences"):
        cif_path = structure_dir / f"{pdb_id}.cif"
        if not cif_path.exists():
            continue

        try:
            structure = parser.get_structure(pdb_id, str(cif_path))
        except Exception as e:
            # This handles common Biopython errors with mmCIF files missing optional fields or other parsing issues.
            print(f"  [WARN] Skipping {pdb_id} due to parsing error: {e}", file=sys.stderr)
            skipped_files.append(f"{pdb_id}: {e}")
            continue

        for chain in structure.get_chains():
            # Heuristic to identify DNA chains by checking for the presence of canonical DNA residue names.
            # This includes standard 3-letter codes (e.g., 'DA') and 1-letter codes.
            dna_resnames = {"DA", "DT", "DG", "DC", "A", "T", "G", "C"}
            resnames = {res.get_resname().strip() for res in chain.get_residues()}
            
            # If there's no intersection between residues in the chain and our DNA set, skip it.
            if not resnames.intersection(dna_resnames):
                continue

            # Build the sequence string by taking the last character of each DNA residue's name (e.g., 'A' from 'DA').
            sequence_chars = [res.get_resname().strip()[-1] for res in chain if res.get_resname().strip() in dna_resnames]
            sequence = "".join(sequence_chars)
            if sequence:
                # Create a unique ID for each chain: PDBID_ChainID
                record_id = f"{pdb_id}_{chain.id}"
                records.append(SeqRecord(Seq(sequence), id=record_id, description=""))

    SeqIO.write(records, FASTA_FILE, "fasta")
    print(f"  -> Wrote {len(records)} DNA chain sequences to {FASTA_FILE}")

    if skipped_files:
        print(f"  -> Skipped {len(skipped_files)} structures due to parsing errors. See {SKIPPED_LOG_FILE.name} for details.")
        with open(SKIPPED_LOG_FILE, "w") as log_f:
            log_f.write("PDB_ID: Error\n")
            log_f.write("\n".join(skipped_files))

def run_cd_hit():
    """Runs the cd-hit-est command."""
    print("\n--- Step 2: Running CD-HIT to cluster sequences ---")
    command = [
        'cd-hit-est',
        '-i', str(FASTA_FILE),
        '-o', str(CDHIT_OUTPUT_PREFIX), # Output prefix
        '-c', str(CDHIT_IDENTITY),
        '-n', '8' # Word size
    ]
    try:
        print(f"  -> Running command: {' '.join(command)}")
        result = subprocess.run(command, check=True, capture_output=True, encoding='utf-8')
        print(f"  -> CD-HIT clustering complete. Output: {CDHIT_OUTPUT}")
        if result.stdout:
            print("  CD-HIT Output:\n", result.stdout)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print("\n[ERROR] CD-HIT execution failed.")
        print("Please ensure 'cd-hit-est' is installed and in your system's PATH.")
        if isinstance(e, subprocess.CalledProcessError):
            print("  CD-HIT STDERR:\n", e.stderr)
        raise

def select_representatives(metadata):
    """Parses CD-HIT clusters and selects the best PDB from each."""
    print("\n--- Step 3: Selecting best representative from each cluster ---")
    representatives = set()
    with open(CDHIT_OUTPUT, 'r') as f:
        clusters = f.read().split('>Cluster')

    for cluster_text in clusters:
        if not cluster_text.strip():
            continue
        
        lines = cluster_text.strip().split('\n')
        cluster_pdb_ids = []
        # Robustly parse PDB IDs from cluster member lines (e.g., "0  123nt, >1ABC_A... *")
        for line in lines:
            match = re.search(r'>([a-zA-Z0-9]{4})_', line)
            if match:
                # FIX 1: Convert PDB ID to uppercase to match metadata.json keys.
                cluster_pdb_ids.append(match.group(1).upper())
        
        # Find the PDB in this cluster with the best resolution
        best_pdb = None
        best_res = float('inf')
        for pdb_id in set(cluster_pdb_ids):
            res_list = metadata.get(pdb_id, {}).get('resolution')
            
            # FIX 2: Extract the numeric value from the list (e.g., [2.5]) for comparison.
            if isinstance(res_list, list) and res_list and isinstance(res_list[0], (float, int)):
                if res_list[0] < best_res:
                    best_res = res_list[0]
                    best_pdb = pdb_id
        
        if best_pdb:
            representatives.add(best_pdb)
            
    print(f"  -> Selected {len(representatives)} non-redundant PDB structures.")
    return {pdb_id: metadata[pdb_id] for pdb_id in representatives}

def main():
    check_dependencies()

    with open(METADATA_FILE, "r") as f:
        all_metadata = json.load(f)
    
    extract_dna_sequences(all_metadata)
    run_cd_hit()
    nonredundant_metadata = select_representatives(all_metadata)

    with open(FINAL_METADATA_FILE, "w") as f:
        json.dump(nonredundant_metadata, f, indent=2)
    print(f"\nSuccess! Non-redundant metadata saved to {FINAL_METADATA_FILE.resolve()}")
    print("You should now run filter_and_prepare.py on this new file.")

if __name__ == "__main__":
    main()
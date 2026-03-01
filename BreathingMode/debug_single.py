
import sys
import argparse
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / 'PDB' / 'code'))
import sht_dna_analysis


def test_single(file_path, num_slices=100):
    file_path = Path(file_path)
    print(f"Testing on {file_path}...")

    try:
        print("Calling compute_sigma_sht_from_pdb...")
        results = sht_dna_analysis.compute_sigma_sht_from_pdb(str(file_path), num_slices=num_slices)
        print("Success!")
        if 'beta_h' in results:
            print(f"Beta head: {results['beta_h'][:5]}")
        if 'A_h' in results:
            print(f"Area head: {results['A_h'][:5]}")
        if 'h' in results:
            print(f"h head: {results['h'][:5]}")

    except Exception as e:
        print(f"Failed: {e}")
        import traceback
        traceback.print_exc()


def main():
    parser = argparse.ArgumentParser(description='Debug run of compute_sigma_sht_from_pdb on one PDB/CIF')
    parser.add_argument('pdb', nargs='?', default=str(PROJECT_ROOT / 'PDB' / 'dataset' / 'hard_dataset' / 'structures' / '102D.cif'),
                        help='Path to PDB/CIF file')
    parser.add_argument('--slices', type=int, default=100, help='Number of slices for SHT')
    args = parser.parse_args()
    test_single(args.pdb, num_slices=args.slices)


if __name__ == '__main__':
    main()

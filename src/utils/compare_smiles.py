#!/usr/bin/env python3
"""
check_smiles.py

Check SMILES molecules for duplicates across multiple CSV files,
or compare two specific files for differences.

Usage:
1. Scan all smiles.csv files recursively in a directory (default: current dir):
    python check_smiles.py

2. Scan a specific directory for duplicates:
    python check_smiles.py --dir path/to/directory

3. Compare two specific CSV files (pairwise comparison):
    python check_smiles.py --file1 file1.csv --file2 file2.csv

Notes:
- Requires pandas and rdkit.
- SMILES columns are detected automatically if not specified.
- Reports duplicates or differences in a readable format.
"""

import os
import pandas as pd
import argparse
from rdkit import Chem

def canonical_smiles(smiles):
    """Return canonical SMILES string or None if invalid."""
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol else None

def read_smiles_from_file(filepath, column=None):
    """Read SMILES from a CSV file, returns list of canonical SMILES."""
    try:
        df = pd.read_csv(filepath)
    except Exception as e:
        print(f"Failed to read {filepath}: {e}")
        return []

    # If column not specified, pick the first column that contains 'smiles' (case-insensitive)
    if column is None:
        for i, col_name in enumerate(df.columns):
            if 'smiles' in col_name.lower():
                column = col_name
                break
        else:
            print(f"No SMILES column found in {filepath}")
            return []

    smiles_list = []
    for s in df[column]:
        if isinstance(s, str):
            can = canonical_smiles(s)
            if can:
                smiles_list.append(can)
    return smiles_list

def find_duplicates_recursively(start_dir='.'):
    """Scan all smiles.csv files recursively and report duplicates across files."""
    seen = {}
    for root, dirs, files in os.walk(start_dir):
        for file in files:
            if file.lower() == 'smiles.csv':
                filepath = os.path.join(root, file)
                smiles_list = read_smiles_from_file(filepath)
                for s in smiles_list:
                    seen.setdefault(s, set()).add(filepath)

    duplicates_found = False
    print("\n=== Duplicate Molecules Across Files ===")
    for smi, sources in seen.items():
        if len(sources) > 1:
            duplicates_found = True
            print(f"Molecule: {smi}")
            for src in sorted(sources):
                print(f"  - {src}")
            print()
    if not duplicates_found:
        print("No duplicates found across files.")

def compare_two_files(file1, file2):
    """Compare SMILES in two files and report differences."""
    smiles1 = read_smiles_from_file(file1)
    smiles2 = read_smiles_from_file(file2)

    set1 = set(smiles1)
    set2 = set(smiles2)

    common = set1 & set2
    only_in_file1 = set1 - set2
    only_in_file2 = set2 - set1

    print(f"\n=== Comparing {file1} and {file2} ===")
    print(f"Total molecules in {file1}: {len(smiles1)}")
    print(f"Total molecules in {file2}: {len(smiles2)}")
    print(f"Molecules in both files: {len(common)}")
    print(f"Molecules only in {file1}: {len(only_in_file1)}")
    print(f"Molecules only in {file2}: {len(only_in_file2)}")

    if only_in_file1:
        print("\nOnly in", file1)
        for s in list(only_in_file1)[:20]:
            print("  ", s)
        if len(only_in_file1) > 20:
            print(f"  ...and {len(only_in_file1)-20} more")
    if only_in_file2:
        print("\nOnly in", file2)
        for s in list(only_in_file2)[:20]:
            print("  ", s)
        if len(only_in_file2) > 20:
            print(f"  ...and {len(only_in_file2)-20} more")

def main():
    parser = argparse.ArgumentParser(description="Check SMILES for duplicates or compare two files.")
    parser.add_argument("--file1", type=str, help="First file for pairwise comparison")
    parser.add_argument("--file2", type=str, help="Second file for pairwise comparison")
    parser.add_argument("--dir", type=str, default='.', help="Directory to scan for duplicates (default current dir)")
    args = parser.parse_args()

    if args.file1 and args.file2:
        compare_two_files(args.file1, args.file2)
    else:
        find_duplicates_recursively(args.dir)

if __name__ == "__main__":
    main()


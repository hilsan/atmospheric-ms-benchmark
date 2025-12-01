#!/usr/bin/env python3
"""
remove_duplicate_SMILES_entries.py

Remove duplicate molecules from a CSV file based on canonical SMILES.
Keeps the first occurrence of each molecule and reindexes the dataset.
Saves a log of duplicates including their original SMILES and canonical form.

Usage:
    python remove_duplicate_SMILES_entries.py -i input.csv -o output.csv [--smiles_col SMILES] [--log_file duplicates.csv]
"""

import argparse
import pandas as pd
from rdkit import Chem

def canonical_smiles(smi):
    """Return canonical SMILES with stereochemistry preserved, or None if invalid."""
    mol = Chem.MolFromSmiles(smi)

    if not mol:
        return None
    # Preserve stereochemistry
    try:
        Chem.Kekulize(mol)
        can_smi = Chem.MolToSmiles(mol, kekuleSmiles=True, isomericSmiles=True)
    except:
        print("kukule error")
        can_smi = Chem.MolToSmiles(mol, kekuleSmiles=True,  isomericSmiles=True)
    return can_smi # Optional: uppercase for readability

def find_duplicates(df, smiles_col='SMILES'):
    """
    Find duplicate molecules based on canonical SMILES.
    Returns:
        to_drop: list of indices to drop
        duplicate_info: list of dicts with original index, kept index, original SMILES, canonical SMILES
        canonical_series: pandas Series of canonical SMILES
    """
    canonical_series = df[smiles_col].apply(lambda s: canonical_smiles(s) if pd.notnull(s) else None)
    
    smiles_to_indices = {}
    for idx, can_smi in canonical_series.items():
        if can_smi is None:
            continue
        smiles_to_indices.setdefault(can_smi, []).append(idx)

    to_drop = set()
    duplicate_info = []
    for indices in smiles_to_indices.values():
        if len(indices) > 1:
            kept = indices[0]
            dropped = indices[1:]
            to_drop.update(dropped)
            for d in dropped:
                duplicate_info.append({
                    'original_duplicate_index': d,
                    'original_kept_index': kept,
                    'original_SMILES': df.at[d, smiles_col],
                    'canonical_SMILES': canonical_series.at[d]
                })

    return sorted(to_drop), duplicate_info, canonical_series

def main():
    parser = argparse.ArgumentParser(description="Remove duplicate SMILES entries in a CSV file.")
    parser.add_argument("-i", "--input_file", required=True, help="Input CSV file with SMILES column")
    parser.add_argument("-o", "--output_file", required=True, help="Output CSV file for unique SMILES")
    parser.add_argument("--smiles_col", default='SMILES', help="Name of the SMILES column (default: SMILES)")
    parser.add_argument("--log_file", default=None, help="Optional CSV file to save duplicate mapping")
    args = parser.parse_args()

    # Read input CSV
    df = pd.read_csv(args.input_file)

    # Find duplicates and canonicalize
    to_drop, duplicate_info, canonical_series = find_duplicates(df, smiles_col=args.smiles_col)
    
    # Replace SMILES with canonical
    df[args.smiles_col] = canonical_series

    # Drop duplicates and reindex
    df_clean = df.drop(index=to_drop).reset_index(drop=True)

    # Save cleaned CSV
    df_clean.to_csv(args.output_file, index=False)
    print(f"Saved {len(df_clean)} unique molecules to {args.output_file}")

    # Save duplicate log if requested
    if args.log_file:
        pd.DataFrame(duplicate_info).to_csv(args.log_file, index=False)
        print(f"Saved duplicate mapping to {args.log_file}")


if __name__ == "__main__":
    main()


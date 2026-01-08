#!/usr/bin/env python3
import argparse
import os, sys
import json
import pandas as pd
sys.path.append("..")
from src.utils import processing

"""


Generate molecule directories and SDF files from a list of TMS-derivatized SMILES.

For each molecule in the input file:
    - Create a directory named by its index or a unique identifier
    - Save the derivatized SMILES as a text file
    - Generate an initial 3D SDF file

Usage (from shell or notebook):

    python make_directories_and_sdfs.py \
        --smiles-file path/to/toy_dataset_TMS.csv \
        --outdir path/to/\
        [--max-iters 200]

Arguments:
    --smiles-file   Path to input CSV file containing derivatized SMILES. 
                    Must contain a column "Modified_SMILES".
    --outdir        Directory where molecule subdirectories will be created.
    --max-iters     Optional. Maximum iterations for 3D geometry optimization (default: 200).

Notes:
    - If a molecule directory already exists, the script will not overwrite it; it only adds the SMILES and SDF files.
    - Designed for use in pipelines or notebooks. You can run directly from a notebook cell using the `!python ...` syntax.
"""


def ensure_dir(path):
    """Create directory if it does not exist."""
    os.makedirs(path, exist_ok=True)

def file_missing(path):
    """Return True if a file does not exist or is empty."""
    return not os.path.isfile(path) or os.path.getsize(path) == 0

def main():
    parser = argparse.ArgumentParser(description="Create molecule folders and SDFs from Modified_SMILES.")
    parser.add_argument("--input_csv", required=True,
                        help="CSV containing Modified_SMILES column.")
    parser.add_argument("--output_root", required=True,
                        help="Directory where molecule folders will be written.")
    parser.add_argument("--max_iters", type=int, default=200,
                        help="Max optimization steps for SDF generation.")
    args = parser.parse_args()

    df = pd.read_csv(args.input_csv)

    if "Modified_SMILES" not in df.columns:
        raise ValueError("Input CSV must contain 'Modified_SMILES' column.")

    ensure_dir(args.output_root)

    for idx, row in df.iterrows():
        smi = row["Modified_SMILES"]
        mol_name = f"{idx:04d}"
        mol_dir = os.path.join(args.output_root, mol_name)

        # Create folder if missing
        if not os.path.exists(mol_dir):
            print(f"[NEW] Creating folder: {mol_name}")
            ensure_dir(mol_dir)
        else:
            print(f"[EXISTING] Folder exists: {mol_name}")

        # File paths
        smi_file = os.path.join(mol_dir, "smiles.smi")
        sdf_file = os.path.join(mol_dir, "structure.sdf")
        metadata_file = os.path.join(mol_dir, "metadata.json")

        # --- Write SMILES file if missing ---
        if file_missing(smi_file):
            print(f"  → Writing SMILES file")
            with open(smi_file, "w") as f:
                f.write(smi)
        else:
            print(f"  → SMILES file exists, skipping")

        # --- Write / regenerate SDF if missing ---
        if file_missing(sdf_file):
            print(f"  → Generating SDF")
            try:
                processing.generate_sdf_from_smiles(
                    smiles_list=[smi],
                    output_file=sdf_file,
                    max_iters=args.max_iters
                )
            except Exception as e:
                print(f"[ERROR] SDF generation failed for {mol_name}: {e}")
        else:
            print(f"  → SDF exists, skipping")

        # --- Write metadata if missing ---
        if file_missing(metadata_file):
            metadata = {
                "index": idx,
                "original_smiles": row.get("Original_SMILES", None),
                "modified_smiles": smi,
                "total_replacements": row.get("Total_Replacements", None)
            }
            with open(metadata_file, "w") as f:
                json.dump(metadata, f, indent=4)
        else:
            print(f"  → metadata.json exists, skipping")

if __name__ == "__main__":
    main()


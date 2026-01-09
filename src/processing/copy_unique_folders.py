#!/usr/bin/env python3
"""
copy_unique_folders.py

Copy folders corresponding to unique molecules based on a duplicate mapping CSV.
This creates a new dataset folder with sequential naming (0000..00XX) for only
the non-duplicate entries.

Usage:
    python copy_unique_folders.py -m duplicate_mapping.csv -s franklin_tms -o franklin_tms_unique
"""

import argparse
import pandas as pd
from pathlib import Path
import shutil

def copy_unique_folders(mapping_file, source_dir, target_dir):
    """
    Copy folders from source_dir to target_dir based on old->new index mapping.

    Parameters:
        mapping_file (str or Path): CSV file with columns 'old_index' and 'new_index'.
        source_dir (str or Path): Folder containing original dataset folders.
        target_dir (str or Path): Folder where unique folders will be copied.
    """
    mapping_file = Path(mapping_file)
    source_dir = Path(source_dir)
    target_dir = Path(target_dir)

    # Create target folder if it does not exist
    target_dir.mkdir(parents=True, exist_ok=True)

    # Load mapping CSV
    mapping = pd.read_csv(mapping_file)
    if 'old_index' not in mapping.columns or 'new_index' not in mapping.columns:
        raise ValueError("Mapping CSV must contain 'old_index' and 'new_index' columns")

    # Copy folders
    for _, row in mapping.iterrows():
        old_idx = int(row['old_index'])
        new_idx = int(row['new_index'])
        old_folder = source_dir / f"{old_idx:04d}"
        new_folder = target_dir / f"{new_idx:04d}"

        if old_folder.exists():
            shutil.copytree(old_folder, new_folder)
            print(f"Copied {old_folder} -> {new_folder}")
        else:
            print(f"Warning: source folder {old_folder} does not exist!")

    print(f"\nFinished copying {len(mapping)} folders to {target_dir}")

def main():
    parser = argparse.ArgumentParser(description="Copy folders for unique molecules based on duplicate mapping.")
    parser.add_argument("-m", "--mapping_file", required=True, help="CSV file with old->new index mapping")
    parser.add_argument("-s", "--source_dir", required=True, help="Source directory containing original folders")
    parser.add_argument("-o", "--target_dir", required=True, help="Target directory for unique dataset folders")
    args = parser.parse_args()

    copy_unique_folders(args.mapping_file, args.source_dir, args.target_dir)

if __name__ == "__main__":
    main()

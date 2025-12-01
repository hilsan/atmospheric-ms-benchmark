import os
import pandas as pd
from rdkit import Chem

# Dictionary to store canonical smiles and their source files
seen_smiles = {}

# Walk through all directories and look for smiles.csv
for root, dirs, files in os.walk('.'):
    for file in files:
        if file == 'smiles.csv':
            filepath = os.path.join(root, file)
            try:
                df = pd.read_csv(filepath)
                for i, row in df.iterrows():
                    # Try to get the SMILES string
                    smiles = row[0] if isinstance(row[0], str) else None
                    if not smiles:
                        continue
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        canonical = Chem.MolToSmiles(mol)
                        if canonical in seen_smiles:
                            seen_smiles[canonical].append(filepath)
                        else:
                            seen_smiles[canonical] = [filepath]
            except Exception as e:
                print(f"Failed to process {filepath}: {e}")

# Report duplicates
print("\n=== Duplicate Molecules Across Files ===")
duplicates_found = False
for smiles, sources in seen_smiles.items():
    if len(set(sources)) > 1:
        duplicates_found = True
        print(f"Molecule: {smiles}")
        for src in sorted(set(sources)):
            print(f"  - {src}")
        print()

if not duplicates_found:
    print("No duplicates found across files.")


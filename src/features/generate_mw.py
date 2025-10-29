import os
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt
import argparse

def count_atoms(smiles, include_hydrogens=True):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    if include_hydrogens:
        mol = Chem.AddHs(mol)
    return mol.GetNumAtoms()

# Argument parsing
parser = argparse.ArgumentParser(description="Calculate molar weight and atom counts for molecules.")
parser.add_argument("-s", "--smiles", type=str, default='smiles.txt', help="Path to SMILES file")
parser.add_argument("-d", "--data_path", type=str, help="Path to read SMILES and write output files")
parser.add_argument("--include_hydrogens", action='store_true', help="Include hydrogens in atom count (default: False)")
args = parser.parse_args()

# Paths
smiles_path = os.path.join(args.data_path, args.smiles)
base_name = os.path.splitext(os.path.basename(smiles_path))[0]
output_file = os.path.join(args.data_path, f"{base_name}_molecular_basic.csv")

# Load data
try:
    df = pd.read_csv(smiles_path)
except FileNotFoundError:
    print(f"Error: The file {smiles_path} was not found.")
    sys.exit(1)

# Ensure SMILES and ID
if 'SMILES' not in df.columns:
    df.columns = ['SMILES']


# Calculate molecular properties
df['molar_mass'] = df['SMILES'].apply(lambda x: ExactMolWt(Chem.MolFromSmiles(x)))
df['atom_count'] = df['SMILES'].apply(lambda x: count_atoms(x, include_hydrogens=args.include_hydrogens))

# Save
df[['SMILES', 'molar_mass', 'atom_count']].to_csv(output_file, index=False)
print(f"Saved molecular mass and atom count to {output_file}")

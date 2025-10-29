#!/usr/bin/env python3
import os
import sys
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
from openbabel.pybel import readstring
from aprl_ssp.util import searchgroups
from aprl_ssp.substructure_search import count_groups
from aprl_ssp.userdef import *

# ---------------------------
# Argument parser
# ---------------------------
parser = ArgumentParser(
    description="Count SIMPOL groups and structural features from a SMILES column.",
    formatter_class=RawTextHelpFormatter
)

parser.add_argument('-g', '--groups', type=str,
                    default='../src/aprl_ssp/SMARTSpatterns/SIMPOLgroups_sane.csv',
                    help='CSV file with SIMPOL SMARTS patterns')
parser.add_argument('-s', '--smiles', type=str, required=True,
                    help='CSV file containing a column with SMILES')
parser.add_argument('-c', '--column', type=str, default='Original_SMILES',
                    help='Column name in the CSV containing SMILES')
parser.add_argument('-o', '--output', type=str, default=None,
                    help='Output CSV path')
parser.add_argument('-d', '--data_path', type=str, default='.',
                    help='Path to the folder containing the SMILES CSV')
args = parser.parse_args()

# ---------------------------
# Paths
# ---------------------------
smiles_path = os.path.join(args.data_path, args.smiles)
smiles_filename = os.path.basename(smiles_path).split('.')[0]
output_file = args.output or os.path.join(args.data_path, f"{smiles_filename}_SIMPOL.csv")

# ---------------------------
# Load data
# ---------------------------
try:
    groups = pd.read_csv(args.groups).set_index('substructure')
    data = pd.read_csv(smiles_path)
except FileNotFoundError as e:
    print(f"File not found: {e.filename}")
    sys.exit(1)

if args.column not in data.columns:
    print(f"Error: Column '{args.column}' not found in {smiles_path}")
    sys.exit(1)

data['SMILES'] = data[args.column]
data['compound'] = data['SMILES']

# ---------------------------
# Count SIMPOL groups
# ---------------------------
search = searchgroups(groups.pattern)
group_counts = count_groups(data.set_index('compound'), search).reset_index(drop=True)

# ---------------------------
# Add metadata
# ---------------------------
group_counts['SMILES'] = data['SMILES']
group_counts['oxygen_count'] = data['SMILES'].str.count('O')
data['mols'] = data['SMILES'].apply(lambda x: readstring('smiles', x))
group_counts['aromatic_ring'] = data['mols'].map(count_aromatic_rings)
group_counts['non_aromatic_ring'] = data['mols'].map(count_nonaromatic_rings)
group_counts['nitrophenol'] = data['mols'].map(
    lambda x: count_nitrophenols(
        x,
        '[OX2H][cX3]:[c]',
        '[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]'
    )
)

# ---------------------------
# Save output
# ---------------------------
group_counts.to_csv(output_file, index=False)
print(f"SIMPOL group counts saved to {output_file}")

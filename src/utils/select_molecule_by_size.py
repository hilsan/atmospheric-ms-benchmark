import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import seaborn as sns
import matplotlib.pyplot as plt

def collect_molecular_weights(root_dir):
    mw_data = []

    for subdir, _, files in os.walk(root_dir):
        for file in files:
            if file == 'SMILES.csv':
                path = os.path.join(subdir, file)
                try:
                    df = pd.read_csv(path)
                    for i, row in df.iterrows():
                        smiles = row[0]  # assuming SMILES is in the first column
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            mw = Descriptors.MolWt(mol)
                            mw_data.append({'subdir': subdir, 'index': i, 'mw': mw, 'smiles': smiles})
                except Exception as e:
                    print(f"Failed to read {path}: {e}")

    return pd.DataFrame(mw_data)

def plot_histogram(mw_df):
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.histplot(mw_df['mw'], kde=True, bins=30, color="skyblue", edgecolor="black")
    plt.title('Molecular Weight Distribution')
    plt.xlabel('Molecular Weight (g/mol)')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.show()

def find_extreme_molecules(mw_df):
    sorted_df = mw_df.sort_values(by='mw').reset_index(drop=True)
    small = sorted_df.iloc[0]
    medium = sorted_df.iloc[len(sorted_df)//2]
    large = sorted_df.iloc[-1]
    return small, medium, large

if __name__ == "__main__":
    root_dir = "."  # change if needed
    mw_df = collect_molecular_weights(root_dir)

    if mw_df.empty:
        print("No molecular weights found.")
    else:
        print(f"Total molecules found: {len(mw_df)}")

        small, medium, large = find_extreme_molecules(mw_df)
        print("\nSelected Molecules:")
        print(f"Smallest: index {small['index']}, MW = {small['mw']:.2f}, SMILES = {small['smiles']}")
        print(f"Median : index {medium['index']}, MW = {medium['mw']:.2f}, SMILES = {medium['smiles']}")
        print(f"Largest: index {large['index']}, MW = {large['mw']:.2f}, SMILES = {large['smiles']}")

        plot_histogram(mw_df)

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
import argparse
import os

def derivatize_molecule(smiles, substitute_secondary, substitute_primary):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)  # Add explicit hydrogens

    replacements = {
        'OH': 0,
        'SH': 0,
        'Secondary Amine': 0,
        'Primary Amine': 0,
        'Imine': 0,
        'OOH': 0,
        'COOH': 0
    }

    def add_tms_group(pattern, replacement_pattern, mol, group_name):
        matches_before = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        num_replacements = len(matches_before)

        if num_replacements > 0:
            mod_mol = Chem.ReplaceSubstructs(mol, 
                                             Chem.MolFromSmarts(pattern), 
                                             Chem.MolFromSmarts(replacement_pattern), 
                                             replaceAll=True)
            mol = mod_mol[0]  # Replace mol with the new modified molecule
            replacements[group_name] += num_replacements
        
        return mol

    mol = add_tms_group('CO-O[H]', 'CO-O[Si](C)(C)(C)', mol, 'OOH')          # Peroxide (OOH)
    mol = add_tms_group('C(=O)O[H]', 'C(=O)(O[Si](C)(C)(C))', mol, 'COOH')  # Carboxylic acid (COOH)
    mol = add_tms_group('CO[H]', 'CO[Si](C)(C)(C)', mol, 'OH')
    mol = add_tms_group('CS[H]', 'CS[Si](C)(C)(C)', mol, 'SH')               # Thiol (SH)
    mol = add_tms_group('C=N([H])', 'C=N[Si](C)(C)(C)', mol, 'Imine')        # Imine (C=N)

    if substitute_primary:
        if not mol.HasSubstructMatch(Chem.MolFromSmarts('aN([H])([H])')):
            mol = add_tms_group('N([H])([H])', 'N([Si](C)(C)(C))([H])', mol, 'Primary Amine')

    if substitute_secondary:
        mol = add_tms_group('N([H])', 'N([Si](C)(C)(C))', mol, 'Secondary Amine')

    mol = Chem.DeleteSubstructs(mol, Chem.MolFromSmarts("[#1]"))

    has_multiple_fragments = len(Chem.GetMolFrags(mol)) > 1

    new_smiles = Chem.MolToSmiles(mol)
    return new_smiles, replacements, has_multiple_fragments

# Define the color scheme for the functional groups (7 colors)
COLORS = {
    'OH': '#9DC858',  # Greenish for OH
    'SH': '#3D6B3F',  # Light blue for SH
    'Secondary Amine': '#E05AB1',  # Pink for Secondary Amine
    'Primary Amine': '#F2A85A',  # Orange for Primary Amine
    'Imine': '#6C5B7B',  # Purple for Imines
    'OOH':'#E05AB1',  # Yellow for OOH
    'COOH': '#8BD2F8',  # Dark Green for COOH
    'black': '#000000',  # Black for borders
    'None': '#7f7f7f'  # Gray for no replacements
}


def plot_number_of_molecules_with_substitutions(df):
    # Add a 'None' column to represent molecules with no substitutions
    df['None'] = (df['Total_Replacements'] == 0).astype(int)

    # Select only the functional group columns
    group_columns = ['OH', 'SH', 'Secondary Amine', 'Primary Amine', 'Imine', 'OOH', 'COOH', 'None']

    # Prepare a DataFrame for the contributions of functional groups to each bar
    plot_data = []

    for _, row in df.iterrows():
        total_replacements = row['Total_Replacements']

        if total_replacements == 0:
            # Molecules with no substitutions contribute to the 'None' group
            plot_data.append({
                'Total_Replacements': 0,
                'Functional Group': 'None',
                'Molecule Contribution': 1.0
            })
        else:
            # Molecules with substitutions: distribute contribution proportionally
            total_group_count = sum(row[group] for group in group_columns if group != 'None')

            for group in group_columns:
                if group != 'None' and row[group] > 0:
                    contribution = row[group] / total_group_count
                    plot_data.append({
                        'Total_Replacements': total_replacements,
                        'Functional Group': group,
                        'Molecule Contribution': contribution
                    })

    # Create a DataFrame for plotting
    plot_df = pd.DataFrame(plot_data)

    # Group and pivot to prepare for plotting
    plot_df_grouped = plot_df.groupby(['Total_Replacements', 'Functional Group'])['Molecule Contribution'].sum().unstack(fill_value=0)

    # Plot the stacked bar chart
    plt.figure(figsize=(12, 8))
    plot_df_grouped.plot(
        kind='bar',
        stacked=True,
        color=[COLORS[group] for group in plot_df_grouped.columns],
        width=0.9
    )

    # Add labels, title, and legend
    plt.xlabel('Number of Substitutions per Molecule', fontsize=14)
    plt.ylabel('Number of Molecules', fontsize=14)
    plt.xticks(rotation=0, fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(title='Functional Groups', fontsize=12, bbox_to_anchor=(1.05, 1), loc='upper left')

    # Adjust layout and save the plot
    plt.tight_layout()
    plt.savefig('count_subs.png', dpi=340)
    plt.show()




def process_smiles_file(input_file, output_file, substitute_secondary, substitute_primary):
    df = pd.read_csv(input_file)

    # Add columns for functional group replacements
    total_replacements = {
        'OH': 0,
        'SH': 0,
        'Secondary Amine': 0,
        'Primary Amine': 0,
        'Imine': 0,
        'OOH': 0,
        'COOH': 0
    }

    processed_data = []

    with open(output_file, 'w') as f:
        f.write("Original_SMILES,Modified_SMILES,Total_Replacements,OH,SH,Primary_Amine,Secondary_Amine,Imine,OOH,COOH\n")

        for original_smiles in df['SMILES']:
            try:
                der_smiles, replacements, has_multiple_fragments = derivatize_molecule(original_smiles, substitute_secondary, substitute_primary)

                total_replacements_count = sum(replacements.values())

                # Append data for plotting
                processed_data.append({
                    'Original_SMILES': original_smiles,
                    'Modified_SMILES': der_smiles,
                    'Total_Replacements': total_replacements_count,
                    **replacements
                })

                f.write(f"{original_smiles},{der_smiles},{total_replacements_count},"
                        f"{replacements['OH']},{replacements['SH']},"
                        f"{replacements['Primary Amine']},{replacements['Secondary Amine']},"
                        f"{replacements['Imine']},{replacements['OOH']},{replacements['COOH']}\n")

            except Exception as e:
                print(f"Error processing SMILES '{original_smiles}': {e}")
                f.write(f"{original_smiles},ERROR\n")

    plot_df = pd.DataFrame(processed_data)

    # Plot the fraction of substitutions per functional group
    plot_number_of_molecules_with_substitutions(plot_df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process SMILES for TMS derivatization.")
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input CSV file containing SMILES.')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to the output CSV file for storing results.')
    parser.add_argument('--substitute_secondary', action='store_false', help='Substitutions of secondary amines.')
    parser.add_argument('--both_on_primary', action='store_true', help='Substitute two hydrogens on primary amines.')

    args = parser.parse_args()

    process_smiles_file(args.input_file, args.output_file, not args.substitute_secondary, not args.both_on_primary)

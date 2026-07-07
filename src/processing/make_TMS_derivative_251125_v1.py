#!/usr/bin/env python3
"""
make_TMS_derivative_v3.py

Derivatize molecules in a CSV file using TMS on functional groups,
plot substitution histogram, and optionally compare with reference TMS counts.

Usage:
    python make_TMS_derivative_251125_v1.py -i input.csv -o output.csv [--compare_ref ref.csv]
    python make_TMS_derivative_251125_v1.py -i input.csv -o output.csv --compare_ref \
        --paper_dir reports/paper/franklin_tms --dataset_name franklin_tms
"""

import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem

# ----------------------------------------
# Plot style
# ----------------------------------------
sns.set_theme(style="whitegrid")
plt.rcParams.update({
    "font.size": 24,
    "axes.labelsize": 24,
    "xtick.labelsize": 24,
    "ytick.labelsize": 24,
    "legend.fontsize": 24
})

COLORS = {
    'OH': '#9DC858',
    'SH': '#3D6B3F',
    'Secondary Amine': '#E05AB1',
    'Primary Amine': '#F2A85A',
    'Imine': '#6C5B7B',
    'OOH': '#E05AB1',
    'COOH': '#8BD2F8',
    'None': '#7f7f7f'
}

# ----------------------------------------
# Derivatization logic
# ----------------------------------------
def derivatize_molecule(smiles, substitute_secondary, substitute_primary):
    mol = Chem.MolFromSmiles(smiles)
    Chem.Kekulize(mol)

    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    
    mol = Chem.AddHs(mol)

    
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
        smarts = Chem.MolFromSmarts(pattern)
        matches = mol.GetSubstructMatches(smarts)
        num_matches = len(matches)

        if num_matches > 0:
            mod_mol = Chem.ReplaceSubstructs(mol, smarts, Chem.MolFromSmarts(replacement_pattern), replaceAll=True)
            mol = mod_mol[0]
            replacements[group_name] += num_matches
        return mol

    # Functional group replacements
    mol = add_tms_group('CO-O[H]', 'CO-O[Si](C)(C)(C)', mol, 'OOH')
    mol = add_tms_group('C(=O)O[H]', 'C(=O)(O[Si](C)(C)(C))', mol, 'COOH')
    mol = add_tms_group('O[H]', 'O[Si](C)(C)(C)', mol, 'OH')
    mol = add_tms_group('CS[H]', 'CS[Si](C)(C)(C)', mol, 'SH')
    mol = add_tms_group('C=N([H])', 'C=N[Si](C)(C)(C)', mol, 'Imine')

    # Amines
    if substitute_primary:
        mol = add_tms_group('N([H])([H])', 'N([Si](C)(C)(C))([H])', mol, 'Primary Amine')
    if substitute_secondary:
        mol = add_tms_group('N([H])', 'N([Si](C)(C)(C))', mol, 'Secondary Amine')
    
    # Remove explicit hydrogens
   # mol = Chem.DeleteSubstructs(mol, Chem.MolFromSmarts("[#1]"))
   
    # SMILES generation
    try:
        Chem.Kekulize(mol)
        mol = Chem.RemoveHs(mol)
        new_smiles = Chem.MolToSmiles(mol, kekuleSmiles=True, allHsExplicit=False, isomericSmiles=True)
    except:
        print('in the except')
        new_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)

    return new_smiles, replacements

# ----------------------------------------
# Plotting function
# ----------------------------------------
def plot_substitutions(processed_data):
    df = pd.DataFrame(processed_data)
    df['None'] = (df['Total_Replacements'] == 0).astype(int)
    group_columns = ['OH', 'SH', 'Secondary Amine', 'Primary Amine', 'Imine', 'OOH', 'COOH', 'None']
    
    plot_data = []
    for _, row in df.iterrows():
        total_replacements = row['Total_Replacements']
        if total_replacements == 0:
            plot_data.append({'Total_Replacements': 0, 'Functional Group': 'None', 'Molecule Contribution': 1.0})
        else:
            total_group_count = sum(row[group] for group in group_columns if group != 'None')
            for group in group_columns:
                if group != 'None' and row[group] > 0:
                    contribution = row[group] / total_group_count
                    plot_data.append({'Total_Replacements': total_replacements, 'Functional Group': group, 'Molecule Contribution': contribution})
    
    plot_df = pd.DataFrame(plot_data)
    plot_df_grouped = plot_df.groupby(['Total_Replacements', 'Functional Group'])['Molecule Contribution'].sum().unstack(fill_value=0)
    
    plt.figure(figsize=(20, 15))
    ax = plot_df_grouped.plot(kind='bar', stacked=True, color=[COLORS[g] for g in plot_df_grouped.columns], edgecolor='black', linewidth=1)
    plt.xlabel('Number of Substitutions per Molecule', labelpad=15)
    plt.ylabel('Fraction of Molecules', labelpad=15)
    plt.xticks(rotation=0)
    plt.legend(title='', loc='upper right', frameon=True, facecolor='white', edgecolor='black')
    plt.tight_layout()
    plt.show()

# ----------------------------------------
# Reference comparison function
# ----------------------------------------
def compare_with_reference(derived_df, reference_df):
    derived_df = derived_df.reset_index(drop=True)
    reference_df = reference_df.reset_index(drop=True)
    
    if 'TMS' not in reference_df.columns:
        print("Warning: Reference file has no 'TMS' column. Skipping comparison.")
        return
    
    s1 = derived_df['Total_Replacements'].reset_index(drop=True)
    s2 = reference_df['TMS'].reset_index(drop=True)

    comparison = s1 == s2
    percentage_equal = comparison.mean() * 100

    print(f"\nPercentage of equality: {percentage_equal:.2f}%")

    mismatches = s1[~comparison]
    mismatched_s2 = s2[~comparison]

    if not mismatches.empty:
        print("\nMismatched rows:")
        for reset_index in mismatches.index:
            smiles1 = derived_df.loc[reset_index, 'Modified_SMILES']
            smiles2 = reference_df.loc[reset_index, 'SMILES']
            print(f"Index {reset_index}: s1 = {mismatches.loc[reset_index]}, s2 = {mismatched_s2.loc[reset_index]}, "
                  f"SMILES s1 = {smiles1}, SMILES s2 = {smiles2}")
    else:
        print("\nNo mismatches found.")

# ── SI LaTeX table ────────────────────────────────────────────────────────────
def write_franklin_si_table(input_df, tms_rows, paper_dir, dataset_name, nist_ids=None):
    """Write a longtable-style SI compound table for Franklin-style datasets.

    Columns: #, Compound name, SMILES, TMS (ref), NIST ID, TMS SMILES, TMS NIST ID.
    Bold rows indicate a mismatch between the reference TMS count and the derived count.
    TMS=0/derived=0 rows show 'not derivatized' spanning the two TMS columns.
    nist_ids: optional dict {index_str: (nist_id_str, tms_nist_id_str)}.
    """
    def esc(t):
        return str(t).replace("_", r"\_").replace("&", r"\&")

    def smi_esc(s):
        # Escape SMILES backslash stereo bonds (\C, \N, etc.) before %, # so
        # they become \textbackslash{}C rather than LaTeX command tokens inside \smi{}.
        return (str(s).replace("\\", r"\textbackslash{}")
                      .replace("#", r"\#")
                      .replace("%", r"\%"))

    nist_ids = nist_ids or {}
    tms_by_smi = {r["Original_SMILES"]: r for r in tms_rows}

    col_spec = r"C{0.6cm} L{3.8cm} L{4.2cm} C{0.6cm} C{1.2cm} L{4.2cm} C{1.2cm}"
    hdr = (r"\textbf{\#} & \textbf{Compound name} & \textbf{SMILES} & "
           r"\textbf{TMS} & \textbf{NIST ID} & \textbf{TMS SMILES} & \textbf{TMS NIST ID} \\")
    label    = f"tab:compounds_{dataset_name}"
    ds_label = esc(dataset_name)

    rows_tex = []
    for i, (_, row) in enumerate(input_df.iterrows()):
        bidx         = f"{i:04d}"
        name         = esc(row["Coumpound.Name.Code"])
        p_smi        = smi_esc(str(row["SMILES"]))
        tms_ref      = int(float(row.get("TMS", 0) or 0))
        tms_r        = tms_by_smi.get(str(row["SMILES"]), {})
        tms_der      = int(float(tms_r.get("Total_Replacements", tms_ref) or tms_ref))
        t_smi        = smi_esc(str(tms_r.get("Modified_SMILES", "")))
        bad          = tms_ref != tms_der
        nist_id, tms_nist_id = nist_ids.get(bidx, ("---", "---"))

        if tms_ref == 0 and tms_der == 0:
            line = (f"{bidx} & {name} & \\smi{{{p_smi}}} & 0 & {nist_id} & "
                    r"\multicolumn{2}{c}{\textit{not derivatized}}")
        elif bad:
            tms_smi_col = f"\\smi{{{t_smi}}}" if t_smi else "---"
            line = (f"\\textbf{{{bidx}}} & \\textbf{{{name}}} & \\smi{{{p_smi}}} & "
                    f"\\textbf{{{tms_ref}}} & \\textbf{{{nist_id}}} & "
                    f"{tms_smi_col} & \\textbf{{{tms_nist_id}}}")
        else:
            tms_smi_col = f"\\smi{{{t_smi}}}" if t_smi else "---"
            line = (f"{bidx} & {name} & \\smi{{{p_smi}}} & "
                    f"{tms_ref} & {nist_id} & {tms_smi_col} & {tms_nist_id}")

        rows_tex.append(line + r" \\[2pt]")

    tex = "\n".join([
        r"% Requires: \usepackage{booktabs,longtable}",
        r"% Uses \smi{} command and C{}/L{} column types from preamble",
        f"\\begin{{longtable}}{{{col_spec}}}",
        f"\\caption{{Dataset compounds ({ds_label}). "
        r"Index corresponds to the numerical identifier used throughout this work. "
        r"TMS gives the number of trimethylsilyl groups as reported in the experimental "
        r"source; the derived TMS count was computed from automated SMILES derivatisation. "
        r"\textbf{Bold rows} indicate a mismatch between reported and derived TMS count. "
        r"For compounds where no derivatisation was applied (TMS\,=\,0, derived\,=\,0), "
        r"TMS columns are omitted. NIST IDs are given where a spectrum was found in the "
        r"NIST/EPA/NIH Mass Spectral Library; --- indicates no spectrum was found.}",
        f"\\label{{{label}}} \\\\",
        r"\toprule",
        hdr,
        r"\midrule",
        r"\endfirsthead",
        f"\\multicolumn{{7}}{{c}}{{\\textit{{Table \\ref{{{label}}} "
        f"continued from previous page}}}} \\\\[4pt]",
        r"\toprule",
        hdr,
        r"\midrule",
        r"\endhead",
        r"\midrule",
        r"\multicolumn{7}{r}{\textit{Continued on next page}} \\",
        r"\endfoot",
        r"\bottomrule",
        r"\endlastfoot",
        "",
        *rows_tex,
        "",
        r"\end{longtable}",
    ])

    out = Path(paper_dir) / f"si_compounds_{dataset_name}.tex"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(tex)
    print(f"Saved SI table: {out}")


# ----------------------------------------
# Main function
# ----------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Process SMILES for TMS derivatization.")
    parser.add_argument('-i', '--input_file', type=str, required=True)
    parser.add_argument('-o', '--output_file', type=str, required=True)
    parser.add_argument('--substitute_secondary', action='store_true', default=False, help="Substitute secondary amines")
    parser.add_argument('--both_on_primary', action='store_true', default=False, help="Substitute both hydrogens on primary amines")
    parser.add_argument('--compare_ref', action='store_true', default=False,
                    help="Compare Total_Replacements with reference TMS column")
    parser.add_argument('--ref_file', type=str, default=None,
                    help="Optional CSV file to use as reference for comparison (default: input file)")
    parser.add_argument('--paper_dir', type=str, default=None,
                    help="Paper output dir for SI table (e.g. reports/paper/franklin_tms)")
    parser.add_argument('--dataset_name', type=str, default=None,
                    help="Dataset name used in SI table label and filename")
    parser.add_argument('--nist_id_file', type=str, default=None,
                    help="CSV with columns Index,NIST_ID,TMS_NIST_ID for SI table")
    args = parser.parse_args()

    df = pd.read_csv(args.input_file)
    processed_data = []

    with open(args.output_file, 'w') as f:
        f.write("Original_SMILES,Modified_SMILES,Total_Replacements,OH,SH,Primary_Amine,Secondary_Amine,Imine,OOH,COOH\n")
        for smi in df['SMILES']:

            try:
                new_smiles, replacements = derivatize_molecule(
                    smi,
                    substitute_secondary=args.substitute_secondary,
                    substitute_primary=args.both_on_primary
                )

                total = sum(replacements.values())
                processed_data.append({'Original_SMILES': smi, 'Modified_SMILES': new_smiles, 'Total_Replacements': total, **replacements})
                f.write(f"{smi},{new_smiles},{total},{replacements['OH']},{replacements['SH']},{replacements['Primary Amine']},{replacements['Secondary Amine']},{replacements['Imine']},{replacements['OOH']},{replacements['COOH']}\n")
            except Exception as e:
                print(f"Error processing SMILES '{smi}': {e}")
                f.write(f"{smi},ERROR\n")

    # Plot substitutions
    plot_substitutions(processed_data)

    # Reference comparison
    if args.compare_ref:
        ref_file = args.ref_file if args.ref_file else args.input_file
        ref_df = pd.read_csv(ref_file)
        compare_with_reference(pd.DataFrame(processed_data), ref_df)

    # SI table for paper
    if args.paper_dir and args.dataset_name:
        nist_ids = {}
        if args.nist_id_file:
            ndf = pd.read_csv(args.nist_id_file, dtype=str).fillna("")
            for _, r in ndf.iterrows():
                idx = f"{int(r['Index']):04d}" if r['Index'].strip().isdigit() else r['Index'].strip()
                nid     = r.get("NIST_ID", "").strip()     or "---"
                tms_nid = r.get("TMS_NIST_ID", "").strip() or "---"
                nist_ids[idx] = (nid, tms_nid)
        write_franklin_si_table(df, processed_data, args.paper_dir, args.dataset_name,
                                nist_ids=nist_ids)

if __name__ == "__main__":
    main()


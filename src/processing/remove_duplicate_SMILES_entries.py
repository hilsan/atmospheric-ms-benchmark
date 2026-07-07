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
from pathlib import Path
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

def write_franklin_si_table(df, paper_dir, dataset_name, nist_ids=None):
    """Write a parent-only SI compound longtable (no TMS derivatisation done).

    Used by the non-TMS Franklin notebook where DERIVATIZE=False.
    TMS SMILES column shows --- for TMS>0 compounds (derivatisation not performed).
    nist_ids: optional dict {index_str: (nist_id_str, tms_nist_id_str)}.
    """
    def esc(t):
        return str(t).replace("_", r"\_").replace("&", r"\&")

    def smi_esc(s):
        return (str(s).replace("\\", r"\textbackslash{}")
                      .replace("#", r"\#")
                      .replace("%", r"\%"))

    nist_ids = nist_ids or {}

    col_spec = r"C{0.6cm} L{3.8cm} L{4.2cm} C{0.6cm} C{1.2cm} L{4.2cm} C{1.2cm}"
    hdr = (r"\textbf{\#} & \textbf{Compound name} & \textbf{SMILES} & "
           r"\textbf{TMS} & \textbf{NIST ID} & \textbf{TMS SMILES} & \textbf{TMS NIST ID} \\")
    label    = f"tab:compounds_{dataset_name}"
    ds_label = esc(dataset_name)

    rows_tex = []
    for i, (_, row) in enumerate(df.iterrows()):
        bidx    = f"{i:04d}"
        name    = esc(row["Coumpound.Name.Code"])
        p_smi   = smi_esc(str(row["SMILES"]))
        tms_ref = int(float(row.get("TMS", 0) or 0))
        nist_id, tms_nist_id = nist_ids.get(bidx, ("---", "---"))

        if tms_ref == 0:
            line = (f"{bidx} & {name} & \\smi{{{p_smi}}} & 0 & {nist_id} & "
                    r"\multicolumn{2}{c}{\textit{not derivatized}}")
        else:
            line = (f"{bidx} & {name} & \\smi{{{p_smi}}} & "
                    f"{tms_ref} & {nist_id} & --- & {tms_nist_id}")

        rows_tex.append(line + r" \\[2pt]")

    tex = "\n".join([
        r"% Requires: \usepackage{booktabs,longtable}",
        r"% Uses \smi{} command and C{}/L{} column types from preamble",
        f"\\begin{{longtable}}{{{col_spec}}}",
        f"\\caption{{Dataset compounds ({ds_label}). "
        r"Index corresponds to the numerical identifier used throughout this work. "
        r"TMS gives the number of trimethylsilyl groups as reported in the experimental "
        r"source. This dataset was analysed without TMS derivatisation (DERIVATIZE\,=\,False); "
        r"TMS SMILES columns are therefore not applicable. "
        r"NIST IDs are given where a spectrum was found in the "
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


def main():
    parser = argparse.ArgumentParser(description="Remove duplicate SMILES entries in a CSV file.")
    parser.add_argument("-i", "--input_file", required=True, help="Input CSV file with SMILES column")
    parser.add_argument("-o", "--output_file", required=True, help="Output CSV file for unique SMILES")
    parser.add_argument("--smiles_col", default='SMILES', help="Name of the SMILES column (default: SMILES)")
    parser.add_argument("--log_file", default=None, help="Optional CSV file to save duplicate mapping")
    parser.add_argument("--index_map_file", default=None, help="Optional CSV file to save old->new index mapping")
    parser.add_argument("--paper_dir", type=str, default=None,
                        help="Paper output dir for SI table (e.g. reports/paper/franklin)")
    parser.add_argument("--dataset_name", type=str, default=None,
                        help="Dataset name used in SI table label and filename")
    parser.add_argument("--nist_id_file", type=str, default=None,
                        help="CSV with columns Index,NIST_ID,TMS_NIST_ID for SI table")
    args = parser.parse_args()

    # Read input CSV
    df = pd.read_csv(args.input_file)

    # Find duplicates and canonicalize
    to_drop, duplicate_info, canonical_series = find_duplicates(df, smiles_col=args.smiles_col)
    
    # Replace SMILES with canonical
    df[args.smiles_col] = canonical_series

    # Drop duplicates
    df_clean = df.drop(index=to_drop).reset_index(drop=True)

    # Build old->new index mapping
    old_to_new_index = {old_idx: new_idx for new_idx, old_idx in enumerate(df_clean.index)}
    # Actually, df_clean.index is already reset, we need mapping from original df indices
    kept_indices = sorted(set(df.index) - set(to_drop))
    old_to_new_index = {old_idx: new_idx for new_idx, old_idx in enumerate(kept_indices)}

    # Save cleaned CSV
    df_clean.to_csv(args.output_file, index=False)
    print(f"Saved {len(df_clean)} unique molecules to {args.output_file}")

    # Save duplicate log if requested
    if args.log_file:
        pd.DataFrame(duplicate_info).to_csv(args.log_file, index=False)
        print(f"Saved duplicate mapping to {args.log_file}")

    # Save old->new index mapping if requested
    if args.index_map_file:
        pd.DataFrame(list(old_to_new_index.items()), columns=['old_index', 'new_index']).to_csv(args.index_map_file, index=False)
        print(f"Saved old->new index mapping to {args.index_map_file}")

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
        write_franklin_si_table(df_clean, args.paper_dir, args.dataset_name, nist_ids=nist_ids)


if __name__ == "__main__":
    main()


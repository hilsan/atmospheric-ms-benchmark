#!/usr/bin/env python3

import argparse
import os
import pandas as pd
from rdkit import Chem
import pubchempy as pcp
import re

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def file_missing(path):
    return not os.path.isfile(path) or os.path.getsize(path) == 0

def smiles_to_inchikey(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return Chem.MolToInchiKey(mol)

def fetch_cas_from_pubchem(identifier: str, id_type='inchikey') -> str:
    """
    Fetch first CAS number from PubChem using InChIKey or SMILES.
    Looks through synonyms using regex to extract CAS numbers.
    """
    try:
        compounds = pcp.get_compounds(identifier, id_type)
        if not compounds:
            return ""

        cid = compounds[0].cid
        # Get all synonyms
        synonym_dicts = pcp.get_synonyms(cid)
        all_synonyms = []
        for d in synonym_dicts:
            all_synonyms.extend(d.get('Synonym', []))

        # Regex to find CAS numbers (format: XX-XX-X up to 7 digits in first block)
        cas_pattern = re.compile(r"\b\d{2,7}-\d{2}-\d\b")
        for syn in all_synonyms:
            match = cas_pattern.search(syn)
            if match:
                return match.group()

        return ""
    except Exception as e:
        print(f"[ERROR] PubChem lookup failed for {identifier}: {e}")
        return ""

def main():
    parser = argparse.ArgumentParser(description="Create EXP directories from a CSV of SMILES.")
    parser.add_argument("--input_csv", required=True, help="CSV file containing SMILES")
    parser.add_argument("--output_root", required=True, help="Directory where EXP folders (0000, 0001, ...) will be created")
    parser.add_argument("--smiles_column", default="SMILES", help="Column containing SMILES (default: SMILES)")
    args = parser.parse_args()

    df = pd.read_csv(args.input_csv)
    if args.smiles_column not in df.columns:
        raise ValueError(f"CSV must contain column '{args.smiles_column}'. Available columns: {list(df.columns)}")

    ensure_dir(args.output_root)

    for idx, row in df.iterrows():
        smi = str(row[args.smiles_column]).strip()
        mol_dir = os.path.join(args.output_root, f"{idx:04d}")
        ensure_dir(mol_dir)

        smiles_file = os.path.join(mol_dir, "smiles.txt")
        inchikey_file = os.path.join(mol_dir, "inchikey.txt")
        cas_file = os.path.join(mol_dir, "cas.txt")

        # Write SMILES
        if file_missing(smiles_file):
            with open(smiles_file, "w") as f:
                f.write(smi + "\n")

        # Generate InChIKey
        if file_missing(inchikey_file):
            try:
                ik = smiles_to_inchikey(smi)
                with open(inchikey_file, "w") as f:
                    f.write(ik + "\n")
            except Exception as e:
                manual_ik_file = os.path.join(mol_dir, "manual_inchikey.txt")
                if os.path.isfile(manual_ik_file):
                    print(f"[INFO] Using manual InChIKey for {mol_dir}")
                    with open(manual_ik_file) as mf, open(inchikey_file, "w") as f:
                        f.write(mf.read().strip() + "\n")
                else:
                    print(f"[WARNING] Could not generate InChIKey for {mol_dir}: {e}")
                    continue

        # Generate CAS number
        if file_missing(cas_file):
            if os.path.isfile(inchikey_file):
                with open(inchikey_file) as f:
                    ik = f.read().strip()
                cas_number = fetch_cas_from_pubchem(ik, 'inchikey')
            else:
                cas_number = fetch_cas_from_pubchem(smi, 'smiles')
            with open(cas_file, "w") as f:
                f.write(cas_number + "\n")

if __name__ == "__main__":
    main()

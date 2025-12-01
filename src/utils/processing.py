"""
Reusable utility functions for spectrum processing and SMILES handling.
All functions here are generic and can be imported by multiple scripts.
"""

import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter

# -----------------------------------------
# Spectrum loading and processing
# -----------------------------------------
def parse_msp(filepath):
    m_z_values, intensities = [], []
    with open(filepath, 'r') as file:
        in_ions_section = False
        for line in file:
            line = line.strip()
            if line.startswith('Begin Ions'):
                in_ions_section = True
                continue
            if line.startswith('End Ions'):
                in_ions_section = False
                continue
            if in_ions_section and line:
                try:
                    mz, intensity = map(float, line.split())
                    m_z_values.append(mz)
                    intensities.append(intensity)
                except ValueError:
                    continue
    return np.array(m_z_values), np.array(intensities)


def parse_csv(filepath):
    m_z_values, intensities = [], []
    with open(filepath, 'r') as file:
        for line in file:
            try:
                mz, intensity = map(float, line.strip().split(','))
                m_z_values.append(mz)
                intensities.append(intensity)
            except ValueError:
                continue
    return np.array(m_z_values), np.array(intensities)


def parse_sdf(filepath):
    m_z_values, intensities = [], []
    with open(filepath, 'r') as file:
        for line in file:
            if "V2000" in line or line.startswith("M  END"):
                continue
            try:
                mz, intensity = map(float, line.strip().split())
                m_z_values.append(mz)
                intensities.append(intensity)
            except ValueError:
                continue
    return np.array(m_z_values), np.array(intensities)


def scale_intensities(intensities, target_max=999):
    """Scale intensities so that the maximum intensity is target_max."""
    max_intensity = max(intensities) if len(intensities) > 0 else 0
    return [intensity / max_intensity * target_max if max_intensity > 0 else intensity for intensity in intensities]


def write_spectrum(filepath, m_z_values, intensities, fmt="msp"):
    """Write spectrum data to a file in MSP, CSV, or SDF format."""
    with open(filepath, 'w') as file:
        if fmt == "msp":
            file.write("Begin Ions\n")
            for mz, intensity in zip(m_z_values, intensities):
                file.write(f"{mz} {intensity}\n")
            file.write("End Ions\n")
        elif fmt == "csv":
            for mz, intensity in zip(m_z_values, intensities):
                file.write(f"{mz},{intensity}\n")
        elif fmt == "sdf":
            for mz, intensity in zip(m_z_values, intensities):
                file.write(f"{mz} {intensity}\n")
        else:
            raise ValueError(f"Unsupported format: {fmt}")


def get_file_format(filepath):
    _, ext = os.path.splitext(filepath)
    ext = ext.lower()
    if ext == ".csv":
        return "csv"
    elif ext == ".msp":
        return "msp"
    elif ext == ".sdf":
        return "sdf"
    else:
        raise ValueError(f"Unsupported file format: {ext}")


# -----------------------------------------
# SMILES / SDF handling
# -----------------------------------------
def generate_sdf_from_smiles(input_file, output_file, error_file=None, max_iters=200):
    """
    Generate an SDF from a SMILES input file. Writes errors to error_file if provided.
    """
    error_file = error_file or os.path.join(os.path.dirname(input_file), "errors.log")
    writer = SDWriter(output_file)
    
    with open(input_file, 'r') as file:
        smiles_list = [line.strip() for line in file]

    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Invalid SMILES: {smiles}")
            mol = Chem.AddHs(mol)
            Chem.Kekulize(mol, clearAromaticFlags=True)
            
            # Embed molecule in 3D
            params = AllChem.ETKDGv3()
            params.useExpTorsionAnglePrefs = False
            params.useBasicKnowledge = False
            if AllChem.EmbedMolecule(mol, params) != 0:
                if AllChem.EmbedMolecule(mol, useRandomCoords=True) != 0:
                    raise ValueError("Embedding failed")
            
            # Optimize
            AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
            
            writer.write(mol)
        except Exception as e:
            with open(error_file, 'a') as ef:
                ef.write(f"Error processing SMILES '{smiles}': {e}\n")
    
    writer.close()
    return output_file


def derivatize_molecule(smiles, substitute_secondary=True, substitute_primary=True):
    """Derivatize functional groups with TMS or similar."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    
    mol = Chem.AddHs(mol)

    # Example: Add your SMARTS replacements here if needed
    # ... (reuse derivatization logic from your previous scripts)
    
    # Remove explicit Hs before returning
    mol = Chem.DeleteSubstructs(mol, Chem.MolFromSmarts("[#1]"))
    
    try:
        new_smiles = Chem.MolToSmiles(mol, kekuleSmiles=True)
    except:
        new_smiles = Chem.MolToSmiles(mol)
    
    return new_smiles


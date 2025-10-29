import os
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SDWriter

def generate_sdf_from_smiles(input_file, output_file, error_file, max_iters=200):
    # Read the SMILES strings from the input file
    with open(input_file, 'r') as file:
        smiles_list = [line.strip() for line in file.readlines()]
    
    # Open the SDF writer
    writer = SDWriter(output_file)
    
    for smiles in smiles_list:
        try:
            # Create a molecule object from the SMILES
            mol = Chem.MolFromSmiles(smiles)
            
            # Check if the molecule was created successfully
            if mol is None:
                raise ValueError(f"Invalid SMILES: {smiles}")
            
            # Add hydrogens to the molecule while preserving aromaticity
            mol_with_h = Chem.AddHs(mol)
            
            # Ensure the molecule has its aromaticity preserved (sometimes necessary after adding hydrogens)
            Chem.Kekulize(mol_with_h, clearAromaticFlags=True)
            mol = mol_with_h 
            # Attempt 3D embedding
            try:
                params = AllChem.ETKDGv3()
                params.useExpTorsionAnglePrefs = False
                params.useBasicKnowledge = False
                if AllChem.EmbedMolecule(mol, params) != 0:
                    print(f"Falling back to random coordinates for {smiles}")
                    if AllChem.EmbedMolecule(mol, useRandomCoords=True) != 0:
                        raise ValueError("Embedding failed")
            except Exception as e:
                raise ValueError(f"3D embedding failed for SMILES: {smiles} with error {str(e)}")
            
            # Optimize the structure (log warnings if it doesn't converge)
            try:
                result = AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
                if result != 0:
                    print(f"Warning: UFF optimization did not converge for {smiles} (result code: {result}).")
            except Exception as e:
                print(f"UFF optimization failed for {smiles} with error: {str(e)}")
            
            # Write the molecule to the SDF file regardless of optimization success
            writer.write(mol)
        
        except Exception as e:
            # Write the error to the error file
            with open(error_file, 'a') as ef:
                ef.write(f"Error processing SMILES '{smiles}': {str(e)}\n")
    
    # Close the writer
    writer.close()
    print(f"SDF file '{output_file}' has been generated.")

def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Generate SDF files from SMILES strings.")
    parser.add_argument("--input_file", required=True, help="Path to the input file containing SMILES strings (one per line).")
    parser.add_argument(
        "--output_file", 
        help="Path to the output SDF file. Defaults to input file name with .sdf extension.",
        default=None
    )
    parser.add_argument(
        "--max_iters", 
        type=int, 
        help="Maximum number of iterations for MMFF optimization. Default is 200.",
        default=200
    )
    
    args = parser.parse_args()
    
    # Determine the output file name
    input_file = args.input_file
    output_file = args.output_file or os.path.splitext(input_file)[0] + ".sdf"
    
    # Get the directory of the input file for error log
    input_dir = os.path.dirname(input_file)
    error_file = os.path.join(input_dir, "errors.log")
    
    # Generate the SDF file and handle errors
    generate_sdf_from_smiles(input_file, output_file, error_file, max_iters=args.max_iters)

if __name__ == "__main__":
    main()

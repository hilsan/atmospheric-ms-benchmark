import argparse
from src.utils import processing

def main():
    parser = argparse.ArgumentParser(description="Generate SDF files from SMILES strings using utils.")
    parser.add_argument("--input_file", required=True, help="Path to input file with SMILES strings (one per line).")
    parser.add_argument(
        "--output_file",
        default=None,
        help="Path to output SDF file. Defaults to input file name with .sdf extension."
    )
    parser.add_argument(
        "--max_iters",
        type=int,
        default=200,
        help="Maximum number of iterations for MMFF/UFF optimization. Default is 200."
    )
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file or processing.default_sdf_name(input_file)

    # Read SMILES using utils
    smiles_list = processing.load_smiles_file(input_file)

    # Process and write SDF
    processing.generate_sdf(
        smiles_list=smiles_list,
        output_file=output_file,
        max_iters=args.max_iters
    )

if __name__ == "__main__":
    main()


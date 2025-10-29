import re
import argparse

def parse_file1(file_path):
    """
    Parse File 1 to extract m/z and intensity values.
    """
    mz_values = []
    intensities = []
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        in_ions_section = False
        
        for line in lines:
            line = line.strip()
            
            # Look for "m/z Values and Intensities" to start extracting peak data
            if "m/z Values and Intensities" in line:
                in_ions_section = True
                continue
            
            if in_ions_section:
                # Extract m/z and intensity pairs (with potential '|' separators)
                peak_data = re.findall(r"(\d+)\s+(\d+)", line)
                for mz, intensity in peak_data:
                    mz_values.append(int(mz))
                    intensities.append(int(intensity))
                    
    return mz_values, intensities

def save_as_file2(mz_values, intensities, output_file):
    """
    Save m/z and intensity values as File 2 format, without extra metadata.
    """
    with open(output_file, 'w') as file:
        # Begin ions section
        file.write("Begin Ions\n")
        
        for mz, intensity in zip(mz_values, intensities):
            file.write(f"{mz} {intensity}\n")
        
        # End ions section
        file.write("End Ions\n")
        
    print(f"File saved to {output_file}")

def main(input_file, output_file):
    # Parse file1
    mz_values, intensities = parse_file1(input_file)
    
    # Save the parsed m/z and intensities into the format of file2
    save_as_file2(mz_values, intensities, output_file)

if __name__ == "__main__":
    # Set up argparse to handle input and output file paths
    parser = argparse.ArgumentParser(description="Parse File 1 and convert to File 2 format.")
    parser.add_argument('-i', '--input', type=str, required=True, help="Path to the input File 1")
    parser.add_argument('-o', '--output', type=str, default='exp.msp', help="Path to output File 2 (default: 'exp.msp')")
    
    args = parser.parse_args()
    
    # Call main with the input and output files from argparse
    main(args.input, args.output)

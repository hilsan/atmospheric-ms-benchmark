import argparse
import os

def scale_intensities(intensities, target_max=999):
    """Scale intensities so that the maximum intensity is target_max."""
    max_intensity = max(intensities) if intensities else 0
    return [intensity / max_intensity * target_max if max_intensity > 0 else intensity for intensity in intensities]

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
    return m_z_values, intensities

def parse_csv(filepath):
    m_z_values, intensities = [], []
    with open(filepath, 'r') as file:
        for line in file:
            try:
                mz, intensity = map(float, line.strip().split(','))
                m_z_values.append(mz)
                intensities.append(intensity)
            except ValueError:
                print(f"Skipping invalid line in CSV: {line.strip()}")
    return m_z_values, intensities

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
                print(f"Skipping invalid line in SDF: {line.strip()}")
    return m_z_values, intensities

def write_spectrum(filepath, m_z_values, intensities, format):
    """Writes scaled spectrum data to a new file."""
    with open(filepath, 'w') as file:
        if format == "msp":
            file.write("Begin Ions\n")
            for mz, intensity in zip(m_z_values, intensities):
                file.write(f"{mz} {intensity}\n")
            file.write("End Ions\n")
        elif format == "csv":
            for mz, intensity in zip(m_z_values, intensities):
                file.write(f"{mz},{intensity}\n")
        elif format == "sdf":
            for mz, intensity in zip(m_z_values, intensities):
                file.write(f"{mz} {intensity}\n")

def main():
    parser = argparse.ArgumentParser(description="Scale intensities in a spectrum file.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input spectrum file")
    parser.add_argument("-o", "--output", required=True, help="Path to save the scaled spectrum")
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    _, ext = os.path.splitext(input_file)
    
    if ext.lower() == ".msp":
        mz_values, intensities = parse_msp(input_file)
        format_type = "msp"
    elif ext.lower() == ".csv":
        mz_values, intensities = parse_csv(input_file)
        format_type = "csv"
    elif ext.lower() == ".sdf":
        mz_values, intensities = parse_sdf(input_file)
        format_type = "sdf"
    else:
        raise ValueError(f"Unsupported file format: {ext}")

    scaled_intensities = scale_intensities(intensities, target_max=999)
    write_spectrum(output_file, mz_values, scaled_intensities, format_type)
    print(f"Scaled spectrum saved to {output_file}")

if __name__ == "__main__":
    main()

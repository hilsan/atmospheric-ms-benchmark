import argparse
import numpy as np
import matplotlib.pyplot as plt
import os

# Color mapping for spectrum types
SPECTRUM_COLORS = {
    "Experimental": "#4D4D4D",  # Dark Grey
    "QCxMS": "#1F77B4",         # Blue
    "QCxMS2-xtb": "#2CA02C",    # Green
    "NEIMS": "#E377C2",         # Pink
    "QCxMS2-DFT": "#FF7F0E",    # Orange
    "Default": "#1F77B4"        # Fallback color
}

def load_spectrum_csv(filepath):
    m_z_values, intensities = [], []
    with open(filepath, 'r') as file:
        for line in file:
            try:
                mz, intensity = map(float, line.strip().split(','))
                m_z_values.append(mz)
                intensities.append(intensity)
            except ValueError:
                print(f"Skipping invalid line in CSV: {line.strip()}")
    return np.array(m_z_values), np.array(intensities)

def load_spectrum_msp(filepath):
    m_z_values, intensities = [], []
    with open(filepath, 'r') as file:
        for line in file:
            if line.startswith("Name") or line.startswith("Comment") or "Ions" in line:
                continue
            try:
                mz, intensity = map(float, line.strip().split())
                m_z_values.append(mz)
                intensities.append(intensity)
            except ValueError:
                print(f"Skipping invalid line in MSP: {line.strip()}")
    return np.array(m_z_values), np.array(intensities)

def parse_sdf(filepath):
    m_z_values, intensities = [], []
    reading_spectrum = False

    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            if "PREDICTED SPECTRUM" in line:
                reading_spectrum = True
                continue
            if reading_spectrum and line == "":
                break
            if reading_spectrum:
                try:
                    mz, intensity = map(float, line.split())
                    m_z_values.append(mz)
                    intensities.append(intensity)
                except ValueError:
                    print(f"Skipping invalid line in SDF: {line}")
    return np.array(m_z_values), np.array(intensities)

def filter_peaks(mz_values, intensities, threshold=1):
    mask = intensities >= threshold
    return mz_values[mask], intensities[mask]

def bin_to_integer(mz_values, intensities):
    # Round m/z to nearest integer
    mz_rounded = np.round(mz_values).astype(int)
    unique_mz = np.unique(mz_rounded)
    binned_intensities = np.zeros_like(unique_mz, dtype=float)

    for i, mz_int in enumerate(unique_mz):
        binned_intensities[i] = intensities[mz_rounded == mz_int].sum()

    return unique_mz, binned_intensities

def plot_stacked_spectra(spectrum_data, output_file):
    num_spectra = len(spectrum_data)
    figsize_width = 12
    figsize_height = 6 * num_spectra * 0.8

    fig, axes = plt.subplots(num_spectra, 1, figsize=(figsize_width, figsize_height), sharex=True)

    if num_spectra == 1:
        axes = [axes]

    for i, (mz_values, intensities, spectrum_type) in enumerate(spectrum_data):
        color = SPECTRUM_COLORS.get(spectrum_type, SPECTRUM_COLORS["Default"])
        axes[i].bar(mz_values, intensities, width=1.0, color=color, alpha=0.7, label=spectrum_type)
        axes[i].set_ylabel("Intensity", fontsize=30)
        axes[i].tick_params(axis='x', labelsize=30)
        axes[i].tick_params(axis='y', labelsize=30)
        axes[i].legend(fontsize=24)
        axes[i].set_xlim(min(mz_values) - 1, max(mz_values) + 1)
        if i == num_spectra - 1:
            axes[i].set_xlabel("m/z", fontsize=30)
            axes[i].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight', transparent=True)
    plt.close()
    print(f"Stacked plot saved to {output_file}")

def get_file_format(filepath):
    _, extension = os.path.splitext(filepath)
    if extension.lower() == ".csv":
        return "csv"
    elif extension.lower() == ".msp":
        return "msp"
    elif extension.lower() == ".sdf":
        return "sdf"
    else:
        raise ValueError(f"Unsupported file format: {extension}")

def main():
    parser = argparse.ArgumentParser(description="Filter, bin, normalize, and plot mass spectra.")
    parser.add_argument("--experimental", required=True, help="Path to the Experimental spectrum file")
    parser.add_argument("--qcxms", required=True, help="Path to the QCxMS spectrum file")
    parser.add_argument("--qcxms2_xtb", required=True, help="Path to the QCxMS2-xtb spectrum file")
    parser.add_argument("--qcxms2_dft", required=True, help="Path to the QCxMS2-DFT spectrum file")
    parser.add_argument("--neims", required=True, help="Path to the NEIMS spectrum file")
    parser.add_argument("--output", required=True, help="Path to save the stacked spectrum plot")

    args = parser.parse_args()

    spectrum_data = []
    file_paths_and_types = [
        (args.experimental, "Experimental"),
        (args.qcxms, "QCxMS"),
        (args.qcxms2_xtb, "QCxMS2-xtb"),
        (args.qcxms2_dft, "QCxMS2-DFT"),
        (args.neims, "NEIMS")
    ]

    for file_path, spectrum_type in file_paths_and_types:
        file_format = get_file_format(file_path)
        if file_format == "csv":
            load_func = load_spectrum_csv
        elif file_format == "msp":
            load_func = load_spectrum_msp
        elif file_format == "sdf":
            load_func = parse_sdf

        mz_values, intensities = load_func(file_path)
        intensities = intensities.astype(float)

        # Filter peaks below threshold
        mz_values, intensities = filter_peaks(mz_values, intensities)

        # Bin to nearest integer m/z
        mz_values, intensities = bin_to_integer(mz_values, intensities)

        # Normalize intensities to max 999
        max_intensity = intensities.max()
        print(f"{spectrum_type} (binned & filtered): max = {max_intensity}, first 10 = {intensities[:10]}")
        if max_intensity > 0:
            intensities = intensities / max_intensity * 999

        spectrum_data.append((mz_values, intensities, spectrum_type))

    plot_stacked_spectra(spectrum_data, args.output)

if __name__ == "__main__":
    main()

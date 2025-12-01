#!/usr/bin/env python3
"""
Plot stacked mass spectra from multiple methods (Experimental, QCxMS, NEIMS, etc.)
Usage:
python plot_stacked_MS.py --experimental exp.csv --qcxms qcxms.csv --neims neims.sdf --output stacked.png
"""

import argparse
import matplotlib.pyplot as plt
from utils.plotting import load_spectrum_csv, load_spectrum_msp, load_spectrum_sdf, bin_and_normalize, get_file_format
import numpy as np
import os

# Color mapping for spectrum types
SPECTRUM_COLORS = {
    "Experimental": "#4D4D4D",
    "QCxMS": "#1F77B4",
    "QCxMS2-xtb": "#2CA02C",
    "NEIMS": "#E377C2",
    "QCxMS2-DFT": "#FF7F0E",
    "Default": "#1F77B4"
}

def plot_stacked_spectra(spectrum_data, output_file):
    num_spectra = len(spectrum_data)
    fig, axes = plt.subplots(num_spectra, 1, figsize=(12, 6*num_spectra*0.8), sharex=True)
    if num_spectra == 1:
        axes = [axes]

    for i, (mz_values, intensities, label) in enumerate(spectrum_data):
        color = SPECTRUM_COLORS.get(label, SPECTRUM_COLORS["Default"])
        axes[i].bar(mz_values, intensities, width=1.0, color=color, alpha=0.7, label=label)
        axes[i].set_ylabel("Intensity", fontsize=18)
        axes[i].legend(fontsize=14)
        axes[i].set_xlim(min(mz_values)-1, max(mz_values)+1)

        if i == num_spectra-1:
            axes[i].set_xlabel("m/z", fontsize=18)

    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight', transparent=True)
    plt.close()
    print(f"Stacked plot saved to {output_file}")

def load_and_process_spectrum(filepath, threshold=1.0, top_n=20):
    fmt = get_file_format(filepath)
    loader = {"csv": load_spectrum_csv, "msp": load_spectrum_msp, "sdf": load_spectrum_sdf}[fmt]
    mz, inten = loader(filepath)
    mz, inten = bin_and_normalize(mz, inten, threshold, top_n)
    return mz, inten

def main():
    parser = argparse.ArgumentParser(description="Plot stacked mass spectra.")
    parser.add_argument("--experimental", required=True)
    parser.add_argument("--qcxms", required=False)
    parser.add_argument("--qcxms2_xtb", required=False)
    parser.add_argument("--qcxms2_dft", required=False)
    parser.add_argument("--neims", required=False)
    parser.add_argument("--threshold", type=float, default=1.0)
    parser.add_argument("--top_n", type=int, default=20)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    spectrum_data = []

    for label, filepath in [("Experimental", args.experimental),
                            ("QCxMS", args.qcxms),
                            ("QCxMS2-xtb", args.qcxms2_xtb),
                            ("QCxMS2-DFT", args.qcxms2_dft),
                            ("NEIMS", args.neims)]:
        if filepath:
            mz, inten = load_and_process_spectrum(filepath, args.threshold, args.top_n)
            spectrum_data.append((mz, inten, label))

    plot_stacked_spectra(spectrum_data, args.output)

if __name__ == "__main__":
    main()


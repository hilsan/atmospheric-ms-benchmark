#!/usr/bin/env python3
"""
Plot mirrored mass spectra:
- Experimental spectrum inverted
- One or more simulated spectra upright
Usage:
python plot_mirrored_spectra.py --exp exp_file.msp --simulated sim1.sdf sim2.csv --output_dir ./plots
"""

import argparse
import os
import matplotlib.pyplot as plt
from utils.plotting import load_spectrum_csv, load_spectrum_msp, load_spectrum_sdf, bin_and_normalize, get_file_format

DISPLAY_NAME_MAP = {
    'annotated.sdf': 'NEIMS',
    'qcxms.csv': 'QCxMS'
}

def plot_mirrored(exp_mz, exp_int, sim_mz, sim_int, sim_label='Simulated', exp_label='Experimental', sim_color="#BAE5FC", exp_color="#E79ECD", output_file="mirrored.png"):
    fig, ax = plt.subplots(figsize=(8,6))
    ax.vlines(sim_mz, 0, sim_int, color=sim_color, linewidth=1.5, label=sim_label)
    ax.vlines(exp_mz, 0, -exp_int, color=exp_color, linewidth=1.5, label=exp_label)
    ax.axhline(0, color='black', linewidth=1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend()
    plt.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Plot mirrored spectra.")
    parser.add_argument("--exp", required=True, help="Experimental spectrum file")
    parser.add_argument("--simulated", required=True, nargs='+', help="Simulated spectra files")
    parser.add_argument("--threshold", type=float, default=1.0)
    parser.add_argument("--top_n", type=int, default=20)
    parser.add_argument("--output_dir", default=".")
    args = parser.parse_args()

    exp_format = get_file_format(args.exp)
    load_func = {"csv": load_spectrum_csv, "msp": load_spectrum_msp, "sdf": load_spectrum_sdf}[exp_format]
    exp_mz, exp_int = bin_and_normalize(*load_func(args.exp), threshold=args.threshold, top_n=args.top_n)

    for sim_file in args.simulated:
        sim_format = get_file_format(sim_file)
        load_func = {"csv": load_spectrum_csv, "msp": load_spectrum_msp, "sdf": load_spectrum_sdf}[sim_format]
        sim_mz, sim_int = bin_and_normalize(*load_func(sim_file), threshold=args.threshold, top_n=args.top_n)
        sim_label = DISPLAY_NAME_MAP.get(os.path.basename(sim_file), os.path.basename(sim_file))
        output_file = os.path.join(args.output_dir, f"{os.path.splitext(os.path.basename(sim_file))[0]}_mirrored.png")
        plot_mirrored(exp_mz, exp_int, sim_mz, sim_int, sim_label=sim_label, output_file=output_file)

if __name__ == "__main__":
    main()


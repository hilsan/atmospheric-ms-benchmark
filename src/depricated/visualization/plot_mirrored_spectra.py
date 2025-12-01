#!/usr/bin/env python3
import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import os
import matplotlib as mpl

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],  # or Helvetica/DejaVu Sans
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 1.2,
    "axes.grid": True,
    "grid.color": "#E0E0E0",
    "grid.linestyle": "--",
    "grid.linewidth": 0.8,
    "legend.frameon": True,
    "legend.fancybox": True,
    "legend.framealpha": 0.8,
    "legend.facecolor": "white",
    "legend.edgecolor": "#E0E0E0",
    "figure.dpi": 340
})



display_name_map = {
    'annotated.sdf': 'NEIMS',
    'qcxms.csv': 'QCxMS'
}

def bin_and_normalize(filepath, threshold=1.0, top_n=20):
    ext = os.path.splitext(filepath)[1].lower()
    mz, intensity = [], []
    with open(filepath, 'r') as f:
        for line in f:
            try:
                if ext == '.csv':
                    x, y = map(float, line.strip().split(','))
                else:
                    x, y = map(float, line.strip().split())
                mz.append(x)
                intensity.append(y)
            except:
                continue

    # Bin to nearest integer
    binned = defaultdict(float)
    for m, inten in zip(mz, intensity):
        binned[round(m)] += inten

    if not binned:
        return [], []

    max_int = max(binned.values())
    norm_binned = {mz_: inten / max_int * 999 for mz_, inten in binned.items() if inten / max_int * 999 >= threshold}

    if top_n > 0:
        norm_binned = dict(sorted(norm_binned.items(), key=lambda x: -x[1])[:top_n])

    sorted_mz = sorted(norm_binned)
    sorted_int = [norm_binned[m] for m in sorted_mz]
    return sorted_mz, sorted_int



def plot_mirrored(exp_mz, exp_int, sim_mz, sim_int, sim_filename='simulated', exp_label='Experimental',
                  sim_color="#BAE5FC", exp_color="#E79ECD", output_file="mirrored_spectrum.png"):
    # Map the filename to a display name
    sim_label = display_name_map.get(sim_filename, sim_filename)

    # Set general font size
    plt.rcParams.update({'font.size': 18})

    fig, ax = plt.subplots(figsize=(8,6))  # less wide, taller aspect ratio

    # Simulated upright
    ax.vlines(sim_mz, 0, sim_int, color=sim_color, linewidth=1.5, label=sim_label)
    # Experimental inverted
    ax.vlines(exp_mz, 0, [-i for i in exp_int], color=exp_color, linewidth=1.5, label=exp_label)

    # Draw black reference lines at x=0 and y=0
    ax.axhline(0, color='black', linewidth=1)
    ax.axvline(0, color='black', linewidth=1)

    # Clean style: hide ticks
    ax.set_xticks([])
    ax.set_yticks([])

    # Set limits
    ax.set_xlim(0, max(sim_mz + exp_mz)+1)

    # Remove the default bottom and left spine lines
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Legend
    ax.legend(frameon=True, facecolor="white", framealpha=0.8)

    plt.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Plot mirrored spectra: simulated upright, experimental inverted.")
    parser.add_argument("--exp", required=True, help="Experimental spectrum file (.msp)")
    parser.add_argument("--simulated", required=True, nargs='+', help="Simulated spectrum files (e.g., annotated.sdf)")
    parser.add_argument("--threshold", type=float, default=1.0, help="Minimum normalized intensity (default=1.0)")
    parser.add_argument("--top_n", type=int, default=20, help="Top N peaks to keep (0 = all, default=20)")
    parser.add_argument("--output_dir", default=".", help="Directory to save plots")
    args = parser.parse_args()

    # Parse experimental spectrum
    exp_mz, exp_int = bin_and_normalize(args.exp, args.threshold, args.top_n)

    for sim_file in args.simulated:
        sim_mz, sim_int = bin_and_normalize(sim_file, args.threshold, args.top_n)
        base_name = os.path.basename(sim_file)
        output_file = os.path.join(args.output_dir, f"{os.path.splitext(base_name)[0]}_mirrored.png")
        plot_mirrored(
            exp_mz, exp_int,
            sim_mz, sim_int,
            sim_filename=base_name,  # <--- use this to map display name
            output_file=output_file
        )


if __name__ == "__main__":
    main()

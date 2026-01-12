#!/usr/bin/env python3
"""
postprocessing_manager.py

Automates post-processing of QCxMS/NEIMS spectra:
- Scaling intensities
- Comparing simulated vs experimental spectra
- Plotting similarity histograms
"""

import argparse
import subprocess
import sys
from pathlib import Path

# Import feature scripts directly if possible
FEATURES_PATH = Path(__file__).resolve().parent.parent / "features"
sys.path.insert(0, str(FEATURES_PATH))

# Optional: import functions if available as modules
# from scale_intensities import scale_intensities_main
# from compare_spectra import compare_spectra_main
# from plot_similarity_histograms import plot_similarity_histograms_main

def run_command(cmd, verbose=True):
    """Run a shell command and optionally print output."""
    if verbose:
        print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"⚠ Command failed: {cmd}")
        print(result.stderr)
    return result

def scale_files(directory: Path):
    """Scale QCxMS and QCxMS2 data."""
    # Define input/output mapping
    scale_commands = [
        (directory.glob("*2/*2/peaks.csv"), directory / "scaled_qcxms2-xtb.csv"),
        (directory.glob("*2/*w*/peaks.csv"), directory / "scaled_qcxms2-dft.csv"),
        ([directory / "QCxMS/MS-run/result.csv"], directory / "scaled_qcxms.csv")
    ]

    for inputs, output in scale_commands:
        for infile in inputs:
            run_command(f"python {FEATURES_PATH}/scale_intensities.py -i {infile} -o {output}")

def compare_spectra(directory: Path, nomax=False):
    """Compare simulated spectra to experimental spectra."""
    cmp_script = "compare_spectra_nomax.py" if nomax else "compare_spectra.py"
    outputs = {
        "neims": directory / f"exp2neims{'_nomax' if nomax else ''}.png",
        "qcxms": directory / f"exp2qcxms{'_nomax' if nomax else ''}.png",
        "qcxms_xtb": directory / f"exp2qcxms-xtb{'_nomax' if nomax else ''}.png",
        "qcxms_dft": directory / f"exp2qcxms-dft{'_nomax' if nomax else ''}.png"
    }

    comparisons = [
        (directory / "annotated.sdf", outputs["neims"]),
        (directory / "scaled_qcxms.csv", outputs["qcxms"]),
        (directory / "scaled_qcxms2-xtb.csv", outputs["qcxms_xtb"]),
        (directory / "scaled_qcxms2-dft.csv", outputs["qcxms_dft"])
    ]

    for sim_file, out_file in comparisons:
        run_command(f"python {FEATURES_PATH}/{cmp_script} --simulated {sim_file} "
                    f"--reference {directory / 'exp.msp'} --output {out_file}")

def plot_histograms(base_dir: Path):
    """Plot similarity histograms for all directories."""
    # All peaks
    allpeaks_dir = base_dir / "res-allpeaks"
    allpeaks_dir.mkdir(exist_ok=True)
    run_command(f"python {FEATURES_PATH}/plot_similarity_histograms.py --dir {base_dir} "
                f"--output {allpeaks_dir} --qcxms_flag --qcxms2_flag --qcxms2_dft_flag")

    # No max
    nomax_dir = base_dir / "res-nomax"
    nomax_dir.mkdir(exist_ok=True)
    run_command(f"python {FEATURES_PATH}/plot_similarity_histograms.py --dir {base_dir} "
                f"--output {nomax_dir} --qcxms_flag exp2qcxms_nomax_similarity_scores.csv "
                f"-qcxms2_flag exp2qcxms2_nomax_similarity_scores.csv "
                f"--qcxms2_dft_flag exp2qcxms2_dft_nomax_similarity_scores.csv "
                f"--neims exp2neims_nomax_similarity_scores.csv")

def main():
    parser = argparse.ArgumentParser(description="Post-processing manager for QCxMS/NEIMS spectra")
    parser.add_argument("--dir", type=Path, required=True,
                        help="Base directory containing simulation directories (0*/1*/etc.)")
    parser.add_argument("--nomax", action="store_true",
                        help="Use compare_spectra_nomax instead of compare_spectra")
    parser.add_argument("--skip-scaling", action="store_true",
                        help="Skip intensity scaling step")
    parser.add_argument("--skip-comparison", action="store_true",
                        help="Skip spectra comparison step")
    parser.add_argument("--skip-plots", action="store_true",
                        help="Skip plotting similarity histograms")
    args = parser.parse_args()

    base_dir = args.dir

    # Loop over all directories starting with 0*
    sim_dirs = [d for d in base_dir.iterdir() if d.is_dir() and d.name.startswith("0")]

    for sim_dir in sim_dirs:
        print(f"\n=== Processing directory: {sim_dir} ===")

        if not args.skip_scaling:
            scale_files(sim_dir)
        if not args.skip_comparison:
            compare_spectra(sim_dir, nomax=args.nomax)

    if not args.skip_plots:
        plot_histograms(base_dir)

    print("\n✅ Post-processing completed for all directories.")

if __name__ == "__main__":
    main()


"""
Run spectral comparisons between experimental and simulated spectra
for all molecules in a simulation results directory, then automatically
back-fill Entropy_Similarity using Python 3.8 (ms_entropy).
"""

import os
import subprocess
import sys

# System Python 3.8 on Puhti — has ms_entropy, numpy, pandas in user site-packages.
# ms_entropy requires Python >= 3.8; the notebook kernel runs 3.6.
_PYTHON38 = "/appl/opt/python/3.8.14-gnu850/bin/python3.8"


def run_comparison(sim_base_dir, methods):
    """Run compare_spectra.py then patch_entropy.py for all peak types and methods.

    Parameters
    ----------
    sim_base_dir : str
        Root directory containing one subfolder per method and an EXP subfolder.
    methods : list[str]
        Method names to include (must match folder names under sim_base_dir).
    """
    sim_base_dir = os.path.abspath(sim_base_dir)
    src_root = os.path.abspath(os.path.join(sim_base_dir, "../../../src"))
    compare_script = os.path.join(src_root, "analysis", "compare_spectra.py")
    entropy_script = os.path.join(src_root, "analysis", "patch_entropy.py")

    print(f"Simulation base directory: {sim_base_dir}")
    print(f"Comparison script:         {compare_script}")
    print(f"Methods:                   {methods}")

    cmd = [
        sys.executable, compare_script,
        "--base_dir", sim_base_dir,
        "--methods", *methods,
        "--include_all_peaks",
        "--include_top20",
        "--include_top10",
    ]

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    if result.returncode != 0:
        print("Error running compare_spectra.py:")
        print(result.stderr)
        return
    else:
        print("Comparison complete.")
        if result.stdout:
            print(result.stdout)

    # Back-fill Entropy_Similarity using Python 3.8 (ms_entropy not available in 3.6 kernel).
    if os.path.isfile(_PYTHON38):
        print(f"\nBack-filling Entropy_Similarity with {_PYTHON38} ...")
        ent_result = subprocess.run(
            [_PYTHON38, entropy_script, "--base_dir", sim_base_dir],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
        )
        if ent_result.returncode != 0:
            print("Warning: patch_entropy.py failed — entropy values will be NaN.")
            print(ent_result.stderr)
        else:
            if ent_result.stdout:
                print(ent_result.stdout)
    else:
        print(f"Warning: Python 3.8 not found at {_PYTHON38}.")
        print("Entropy_Similarity will be NaN. Run manually:\n"
              f"  {_PYTHON38} {entropy_script} --base_dir {sim_base_dir}")

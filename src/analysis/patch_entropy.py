#!/usr/bin/env python3
"""
patch_entropy.py

Back-fill the Entropy_Similarity column in all spectra_*_comparison.csv files
for a given simulation base directory, and optionally plot histograms of the
entropy similarity distribution per method.

Run with Python 3.8+ (where ms_entropy >= 1.0 is available):

    /appl/opt/python/3.8.14-gnu850/bin/python3.8 src/analysis/patch_entropy.py \
        --base_dir data/simulation_results/franklin_tms

The script reads each existing comparison CSV, computes entropy similarity from
the corresponding experimental and simulated spectra files, and writes the value
into the Entropy_Similarity column. Files that already have a non-NaN value are
skipped unless --force is passed.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import ms_entropy as _ms_entropy
except ImportError:
    raise SystemExit("ms_entropy not found — run with Python 3.8+:\n"
                     "  module load python-data/3.8-22.10")


# Map peak_type string (as stored in comparison CSVs) to:
#   - the comparison CSV stem
#   - the spectrum file name inside spectra/ subdirectory
PEAK_FILE_PREFIX = {
    "all":   ("spectra_all_comparison",   "spectra_all.csv"),
    "top20": ("spectra_top20_comparison", "spectra_top20.csv"),
    "10pct": ("spectra_10pct_comparison", "spectra_10pct.csv"),
}


def _read_spectrum_file(path):
    """Read mz,intensity CSV into {mz: intensity} dict."""
    if not Path(path).is_file():
        return None
    out = {}
    with open(path) as f:
        for line in f:
            try:
                x, y = map(float, line.strip().split(","))
                out[int(round(x))] = out.get(int(round(x)), 0.0) + y
            except Exception:
                pass
    return out or None


def _entropy_sim(spec_ref, spec_sim):
    if spec_ref is None or spec_sim is None:
        return np.nan
    arr_ref = np.array([[mz, i] for mz, i in spec_ref.items()], dtype=np.float32)
    arr_sim = np.array([[mz, i] for mz, i in spec_sim.items()], dtype=np.float32)
    if len(arr_ref) == 0 or len(arr_sim) == 0:
        return np.nan
    return float(_ms_entropy.calculate_entropy_similarity(
        arr_ref, arr_sim, ms2_tolerance_in_da=0.05, clean_spectra=True
    ))


def patch_base_dir(base_dir: str, force: bool = False):
    base = Path(base_dir).resolve()
    results_dir = base / "results"
    exp_base = base / "EXP"

    if not results_dir.exists():
        print(f"No results/ dir found in {base_dir}")
        return

    mol_dirs = sorted(d for d in results_dir.iterdir() if d.name.isdigit())
    print(f"Found {len(mol_dirs)} molecule result dirs in {results_dir}")

    total_updated = 0

    for mol_dir in mol_dirs:
        mol_idx = mol_dir.name

        for pt, (csv_stem, spectra_fname) in PEAK_FILE_PREFIX.items():
            csv_path = mol_dir / f"{csv_stem}.csv"
            if not csv_path.exists():
                continue

            df = pd.read_csv(csv_path)

            if "Entropy_Similarity" not in df.columns:
                df["Entropy_Similarity"] = np.nan

            rows_to_update = df.index if force else df.index[df["Entropy_Similarity"].isna()]
            if len(rows_to_update) == 0:
                continue

            # Load experimental spectrum for this peak type
            exp_file = exp_base / mol_idx / "spectra" / spectra_fname
            spec_ref = _read_spectrum_file(exp_file)

            updated = 0
            for idx in rows_to_update:
                method = df.at[idx, "Method"]
                sim_file = base / method / mol_idx / "spectra" / spectra_fname
                spec_sim = _read_spectrum_file(sim_file)
                val = _entropy_sim(spec_ref, spec_sim)
                df.at[idx, "Entropy_Similarity"] = val
                updated += 1

            df.to_csv(csv_path, index=False)
            total_updated += updated

    print(f"Done — updated {total_updated} Entropy_Similarity values in {base_dir}")


def main():
    parser = argparse.ArgumentParser(description="Back-fill Entropy_Similarity in comparison CSVs.")
    parser.add_argument("--base_dir", required=True,
                        help="Simulation base dir (e.g. data/simulation_results/franklin_tms)")
    parser.add_argument("--force", action="store_true",
                        help="Recompute even if Entropy_Similarity is already filled")
    args = parser.parse_args()
    patch_base_dir(args.base_dir, force=args.force)


if __name__ == "__main__":
    main()

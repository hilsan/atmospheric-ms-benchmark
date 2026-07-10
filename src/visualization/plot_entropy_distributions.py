#!/usr/bin/env python3
"""
Plot 3x3 grid of Shannon spectral entropy histograms.
Rows = datasets (Franklin TMS, UCB-GLOBES, Franklin underivatized)
Cols = methods (EXP, NEIMS, QCxMS)
"""
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
DATA_ROOT = Path(__file__).resolve().parents[2] / "data" / "simulation_results"

DATASETS = [
    {"label": "GoAmazon (TMS)",   "key": "franklin_tms",       "qcxms": "QCxMS_10_ps"},
    {"label": "External (UCB)",   "key": "ucb_globes_tracers", "qcxms": "QCxMS_10_ps"},
    {"label": "GoAmazon (underivatized)", "key": "franklin",   "qcxms": "QCxMS_10_ps"},
]

METHODS = ["EXP", "NEIMS", "QCxMS"]
METHOD_DIRS = {
    "EXP":   "EXP",
    "NEIMS": "NEIMS",
    "QCxMS": None,   # filled per-dataset from ds["qcxms"]
}

METHOD_COLORS = {
    "EXP":   "#444444",
    "NEIMS": "#D36EA5",
    "QCxMS": "#5A99D3",
}

GREY_TEXT = "#4C4C4C"


# ---------------------------------------------------------------------------
# Shannon entropy
# ---------------------------------------------------------------------------
def shannon_entropy(intensities: np.ndarray) -> float:
    total = intensities.sum()
    if total <= 0:
        return np.nan
    p = intensities / total
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))


def load_entropies(base_dir: Path, method_dir: str) -> np.ndarray:
    """Compute per-molecule Shannon entropy from spectra_all.csv files."""
    spectra_dir = base_dir / method_dir
    if not spectra_dir.exists():
        return np.array([])

    entropies = []
    for mol_dir in sorted(spectra_dir.iterdir()):
        spec_file = mol_dir / "spectra" / "spectra_all.csv"
        if not spec_file.exists():
            continue
        try:
            df = pd.read_csv(spec_file, header=None, names=["mz", "intensity"])
            h = shannon_entropy(df["intensity"].values)
            if not np.isnan(h):
                entropies.append(h)
        except Exception:
            continue
    return np.array(entropies)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def make_figure(output_path: Path):
    fig, axes = plt.subplots(
        len(DATASETS), len(METHODS),
        figsize=(7.08, 6.0),   # ACS 2-col width, ~6 in tall
        sharey=False, sharex=False,
        facecolor="none",
    )
    plt.rcParams.update({"font.size": 7})

    # collect all entropy arrays first (for shared x range per method)
    all_data = {}
    for ds in DATASETS:
        base = DATA_ROOT / ds["key"]
        for method in METHODS:
            mdir = ds["qcxms"] if method == "QCxMS" else METHOD_DIRS[method]
            h = load_entropies(base, mdir)
            all_data[(ds["key"], method)] = h

    for row_idx, ds in enumerate(DATASETS):
        for col_idx, method in enumerate(METHODS):
            ax = axes[row_idx][col_idx]
            h = all_data[(ds["key"], method)]

            if len(h) == 0:
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=6, color="grey")
                ax.set_visible(True)
                continue

            color = METHOD_COLORS[method]
            ax.hist(h, bins=10, color=color, edgecolor="white", linewidth=0.4,
                    alpha=0.85)

            # Annotate mean ± std
            mu, sd = h.mean(), h.std()
            ax.axvline(mu, color=color, linewidth=1.2, linestyle="--", alpha=0.9)
            ax.text(0.97, 0.93,
                    f"$\\bar{{H}}={mu:.2f}$\n$\\sigma={sd:.2f}$",
                    transform=ax.transAxes,
                    ha="right", va="top",
                    fontsize=6, color=GREY_TEXT)

            ax.set_facecolor("none")
            ax.tick_params(colors=GREY_TEXT, labelsize=6)
            for spine in ax.spines.values():
                spine.set_edgecolor("#cccccc")

            # Column header (top row only)
            if row_idx == 0:
                ax.set_title(method, fontsize=7, color=GREY_TEXT, pad=4)

            # Row label (left column only)
            if col_idx == 0:
                ax.set_ylabel(ds["label"], fontsize=6, color=GREY_TEXT, labelpad=4)

            # x-axis label (bottom row only)
            if row_idx == len(DATASETS) - 1:
                ax.set_xlabel("Shannon entropy $H$", fontsize=6, color=GREY_TEXT)

    fig.suptitle("Spectral Shannon entropy distributions", fontsize=8,
                 color=GREY_TEXT, y=1.01)
    plt.tight_layout(pad=0.5, h_pad=0.8, w_pad=0.6)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {output_path}")
    plt.close(fig)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", default=None,
                        help="Output PDF path (default: reports/combined/paper/entropy_distributions.pdf)")
    args = parser.parse_args()

    out = Path(args.output) if args.output else (
        Path(__file__).resolve().parents[2] / "reports" / "combined" / "paper" / "entropy_distributions.pdf"
    )
    make_figure(out)

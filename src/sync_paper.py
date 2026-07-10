#!/usr/bin/env python3
"""
Sync figures, tables, and numbers from benchmark reports/ to the paper repo,
then commit and push to GitHub.

Usage:
    python src/sync_paper.py [--dry-run] [--no-push]
"""

import argparse
import os
import shutil
import subprocess
from pathlib import Path
from datetime import datetime

BENCH = Path("/scratch/project_2006752/hsandstr/Project/atmospheric-ms-benchmark")
PAPER = Path("/scratch/project_2006752/hsandstr/Project/-benchmark-paper")
REPORTS = BENCH / "reports"

# ---------------------------------------------------------------------------
# File mapping: (source relative to REPORTS, dest relative to PAPER)
# ---------------------------------------------------------------------------
FIGURE_MAP = [
    # GoAmazon / TMS dataset figures
    ("franklin_tms/paper/franklin_tms_properties_figure.pdf",  "figures/tms/properties_figure.pdf"),
    ("franklin_tms/paper/franklin_tms_tms_distribution.pdf",   "figures/tms/properties.pdf"),
    ("franklin_tms/paper/franklin_tms_score_vs_properties.pdf","figures/tms/score_vs_properties.pdf"),
    ("franklin_tms/paper/franklin_tms_fg_correlation_heatmap.pdf", "figures/tms/fg_correlation_heatmap.pdf"),
    ("franklin_tms/paper/franklin_tms_fg_prevalence.pdf",      "figures/tms/fg_prevalence.pdf"),
    ("franklin_tms/paper/sensitivity_figure.png",              "figures/tms/sensitivity_figure.png"),
    ("franklin_tms/paper/sensitivity_figure.svg",              "figures/tms/sensitivity_figure.svg"),
    ("franklin_tms/paper/cosine_by_peak_type.pdf",             "figures/tms/cosine_by_peak_type.pdf"),
    # UCB / External dataset figures
    # NOTE: ucb properties_figure.pdf requires running notebook 6 with PAPER_MODE=True
    ("ucb_globes_tracers/paper/sensitivity_figure.png",        "figures/ucb_tracers/sensitivity_figure.png"),
    ("ucb_globes_tracers/paper/cosine_by_peak_type.pdf",       "figures/ucb_tracers/cosine_by_peak_type.pdf"),
    # Franklin underivatized (appendix)
    ("franklin/paper/franklin_properties_figure.pdf",          "figures/orig/properties_figure.pdf"),
    ("franklin/paper/franklin_fg_correlation_heatmap.pdf",     "figures/orig/fg_correlation_heatmap.pdf"),
    ("franklin/paper/franklin_score_vs_properties.pdf",        "figures/orig/score_vs_properties.pdf"),
    ("franklin/paper/sensitivity_figure.png",                  "figures/orig/sensitivity_figure.png"),
    # UMAP / combined figures (from presentation prep)
    ("combined/paper/umap_qcxms_neims.pdf",                    "figures/tms/umap_cosine_score.pdf"),
    ("combined/paper/umap_all_methods.pdf",                    "figures/combined/umap_all_methods.pdf"),
    ("combined/pres/neims_cosine_by_dataset.pdf",              "figures/combined/neims_cosine_by_dataset.pdf"),
    # Hyperparameter scan figures
    ("iee_scan/iee_scan_metrics.png",                          "figures/grid_hyperparam/iee_scan_metrics.png"),
    ("iee_scan/mirror_0039_iee03.png",                         "figures/grid_hyperparam/mirror_0039_iee03.png"),
    # Entropy distribution figure (appendix)
    ("combined/paper/entropy_distributions.pdf",               "figures/combined/entropy_distributions.pdf"),
]

TABLE_MAP = [
    # Main results table (TMS / GoAmazon, all peak-picking strategies)
    ("franklin_tms/paper/table_combined.tex",                  "tables/pres/table_combined.tex"),
    # TMS tables
    ("franklin_tms/paper/table_all.tex",                       "tables/tms/table_all.tex"),
    ("franklin_tms/paper/table_10pct.tex",                     "tables/tms/table_10pct.tex"),
    ("franklin_tms/paper/table_top20.tex",                     "tables/tms/table_top20.tex"),
    ("franklin_tms/paper/table_qcxms_variants.tex",            "tables/tms/table_qcxms_variants.tex"),
    ("franklin_tms/paper/table_strict_overlap_all.tex",        "tables/tms/table_strict_overlap_all.tex"),
    ("franklin_tms/paper/table_strict_overlap_10pct.tex",      "tables/tms/table_strict_overlap_10pct.tex"),
    ("franklin_tms/paper/table_strict_overlap_top20.tex",      "tables/tms/table_strict_overlap_top20.tex"),
    ("bin/paper/franklin_tms/si_compounds_franklin_tms.tex",    "tables/tms/si_compounds_franklin_tms.tex"),
    ("bin/paper/franklin/si_compounds_franklin.tex",            "tables/orig/si_compounds_franklin.tex"),
    ("franklin_tms/paper/diagnostics_coverage.tex",            "tables/tms/diagnostics_coverage.tex"),
    # Franklin underivatized tables
    ("franklin/paper/table_combined.tex",                      "tables/orig/table_combined.tex"),
    ("franklin/paper/table_all.tex",                           "tables/orig/table_all.tex"),
    ("franklin/paper/table_10pct.tex",                         "tables/orig/table_10pct.tex"),
    ("franklin/paper/table_top20.tex",                         "tables/orig/table_top20.tex"),
    ("franklin/paper/table_qcxms_variants.tex",                "tables/orig/table_qcxms_variants.tex"),
    ("franklin/paper/table_strict_overlap_all.tex",            "tables/orig/table_strict_overlap_all.tex"),
    ("franklin/paper/table_strict_overlap_10pct.tex",          "tables/orig/table_strict_overlap_10pct.tex"),
    ("franklin/paper/table_strict_overlap_top20.tex",          "tables/orig/table_strict_overlap_top20.tex"),
    ("franklin/paper/diagnostics_coverage.tex",                "tables/orig/diagnostics_coverage.tex"),
    ("franklin/paper/table_perclass_entropy.tex",              "tables/orig/table_perclass_entropy.tex"),
    # UCB / External dataset tables
    ("ucb_globes_tracers/paper/table_all.tex",                 "tables/tracers/table_all.tex"),
    ("ucb_globes_tracers/paper/table_combined.tex",            "tables/tracers/table_combined.tex"),
    ("ucb_globes_tracers/paper/table_10pct.tex",               "tables/tracers/table_10pct.tex"),
    ("ucb_globes_tracers/paper/table_top20.tex",               "tables/tracers/table_top20.tex"),
    ("ucb_globes_tracers/paper/table_strict_overlap_all.tex",  "tables/tracers/table_strict_overlap_all.tex"),
    # SI compound tables
    ("paper/ucb_globes_tracers/si_compounds_ucb_globes_tracers.tex",
                                                               "tables/tracers/si_compounds_ucb_globes_tracers.tex"),
]

# Files that need action before they exist (warn user)
MISSING_NOTES = [
    ("ucb_globes_tracers/paper/properties_figure.pdf",
     "Run notebook 6 with PAPER_MODE=True for ucb_globes_tracers → figures/ucb_tracers/properties_figure.pdf"),
]


def sync_file(src: Path, dst: Path, dry_run: bool) -> str:
    """Copy src → dst. Returns status string."""
    if not src.exists():
        return f"  MISSING  {src.relative_to(REPORTS)}"
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        src_mtime = src.stat().st_mtime
        dst_mtime = dst.stat().st_mtime
        if src_mtime <= dst_mtime:
            return f"  UP-TO-DATE  {dst.relative_to(PAPER)}"
        status = "UPDATED"
    else:
        status = "NEW"
    if not dry_run:
        shutil.copy2(src, dst)
    return f"  {status}  {dst.relative_to(PAPER)}"


def _load_strict_overlap(base_dir: Path, canonical: str = "QCxMS_10_ps") -> "pd.DataFrame":
    """Load results CSVs, rename canonical QCxMS, restrict to strict-overlap molecules."""
    import pandas as pd
    import numpy as np
    all_qcxms = {"QCxMS_10_ps", "QCxMS_25_ps", "QCxMS_10_ps_iee03"}
    dfs = []
    for f in sorted((base_dir / "results").glob("*/*.csv")):
        if f.parent.name.isdigit():
            df = pd.read_csv(f)
            df["Molecule"] = f.parent.name
            df["Peak_Type"] = f.stem.split("_")[1]
            dfs.append(df)
    if not dfs:
        return None
    data = pd.concat(dfs, ignore_index=True)
    data["Method"] = data["Method"].replace({canonical: "QCxMS"})
    data = data[~data["Method"].isin(all_qcxms - {canonical})]
    # Fill missing metric columns
    for col in ["Cosine", "Weighted_Dot", "Tanimoto", "%S_sim_in_ref", "%R_ref_in_sim", "Entropy_Similarity"]:
        if col not in data.columns:
            data[col] = float("nan")
    # Strict-overlap: molecules with results in all active methods
    methods = [m for m in ["QCxMS", "QCxMS2", "CFMID", "NEIMS"] if m in data["Method"].unique()]
    pivot = data.groupby(["Molecule", "Method"])["Cosine"].mean().unstack("Method")
    complete = pivot.dropna(subset=methods).index.tolist()
    return data[data["Molecule"].isin(complete)]


def _stat(data, method: str, metric: str, peak_type: str = "all") -> tuple:
    """Return (mean, std) rounded to 1 dp; (nan, nan) if no data."""
    import numpy as np
    sub = data[(data["Method"] == method) & (data["Peak_Type"] == peak_type)][metric].dropna()
    if len(sub) == 0:
        return float("nan"), float("nan")
    return round(sub.mean(), 1), round(sub.std(), 1)


def generate_numbers_tex(dry_run: bool):
    """Write numbers.tex by reading result CSVs — always up-to-date after sync."""
    import numpy as np

    DATA = BENCH / "data" / "simulation_results"
    tms  = _load_strict_overlap(DATA / "franklin_tms")
    ucb  = _load_strict_overlap(DATA / "ucb_globes_tracers")
    orig = _load_strict_overlap(DATA / "franklin")

    def nc(name, val):
        v = f"{val:.1f}" if not (isinstance(val, float) and np.isnan(val)) else "NA"
        return f"\\newcommand{{\\{name}}}{{{v}}}"

    lines = [
        "%% numbers.tex — auto-generated by src/sync_paper.py",
        "%% Do NOT edit manually. Re-run sync_paper.py to update.",
        "%% Include in manuscript preamble: \\input{numbers}",
        "%%",
        "%% GoAmazon (Franklin TMS) — all peaks, strict-overlap",
    ]
    if tms is not None:
        n_so = tms["Molecule"].nunique()
        lines += [
            f"%% N = {n_so} strict-overlap molecules",
            nc("neimsCosineTMS",        _stat(tms, "NEIMS",  "Cosine")[0]),
            nc("neimsCosineTMSstd",     _stat(tms, "NEIMS",  "Cosine")[1]),
            nc("neimsWeightedDotTMS",   _stat(tms, "NEIMS",  "Weighted_Dot")[0]),
            nc("neimsTanimotoTMS",      _stat(tms, "NEIMS",  "Tanimoto")[0]),
            nc("neimsRecallTMS",        _stat(tms, "NEIMS",  "%R_ref_in_sim")[0]),
            nc("neimsPrecisionTMS",     _stat(tms, "NEIMS",  "%S_sim_in_ref")[0]),
            nc("neimsCosineTopTwentyTMS", _stat(tms, "NEIMS", "Cosine", "top20")[0]),
            "%%",
            nc("qcxmsCosineTMS",        _stat(tms, "QCxMS",  "Cosine")[0]),
            nc("qcxmsWeightedDotTMS",   _stat(tms, "QCxMS",  "Weighted_Dot")[0]),
            nc("qcxmsTanimotoTMS",      _stat(tms, "QCxMS",  "Tanimoto")[0]),
            nc("qcxmsRecallTMS",        _stat(tms, "QCxMS",  "%R_ref_in_sim")[0]),
            nc("qcxmsPrecisionTMS",     _stat(tms, "QCxMS",  "%S_sim_in_ref")[0]),
            "%%",
            nc("qcxmsTwoCosineTMS",     _stat(tms, "QCxMS2", "Cosine")[0]),
            nc("qcxmsTwoWeightedDotTMS",_stat(tms, "QCxMS2", "Weighted_Dot")[0]),
            nc("qcxmsTwoTanimotoTMS",   _stat(tms, "QCxMS2", "Tanimoto")[0]),
            nc("qcxmsTwoRecallTMS",     _stat(tms, "QCxMS2", "%R_ref_in_sim")[0]),
            nc("qcxmsTwoPrecisionTMS",  _stat(tms, "QCxMS2", "%S_sim_in_ref")[0]),
            "%%",
            nc("cfmidCosineTMS",        _stat(tms, "CFMID",  "Cosine")[0]),
            nc("cfmidWeightedDotTMS",   _stat(tms, "CFMID",  "Weighted_Dot")[0]),
            nc("cfmidTanimotoTMS",      _stat(tms, "CFMID",  "Tanimoto")[0]),
            nc("cfmidRecallTMS",        _stat(tms, "CFMID",  "%R_ref_in_sim")[0]),
            nc("cfmidPrecisionTMS",     _stat(tms, "CFMID",  "%S_sim_in_ref")[0]),
        ]
    lines += ["%%", "%% External dataset (UCB-GLOBES) — all peaks, strict-overlap"]
    if ucb is not None:
        lines += [
            nc("neimsCosinUCB",         _stat(ucb, "NEIMS",  "Cosine")[0]),
            nc("neimsTanimotoUCB",      _stat(ucb, "NEIMS",  "Tanimoto")[0]),
            nc("neimsRecallUCB",        _stat(ucb, "NEIMS",  "%R_ref_in_sim")[0]),
            nc("neimsPrecisionUCB",     _stat(ucb, "NEIMS",  "%S_sim_in_ref")[0]),
            "%%",
            nc("qcxmsCosinUCB",         _stat(ucb, "QCxMS",  "Cosine")[0]),
            nc("qcxmsTanimotoUCB",      _stat(ucb, "QCxMS",  "Tanimoto")[0]),
            nc("qcxmsTwoCosineUCB",     _stat(ucb, "QCxMS2", "Cosine")[0]),
            nc("qcxmsTwoTanimotoUCB",   _stat(ucb, "QCxMS2", "Tanimoto")[0]),
            nc("qcxmsTwoRecallUCB",     _stat(ucb, "QCxMS2", "%R_ref_in_sim")[0]),
            nc("qcxmsTwoPrecisionUCB",  _stat(ucb, "QCxMS2", "%S_sim_in_ref")[0]),
        ]
    lines += ["%%", "%% Non-derivatized (Franklin) — all peaks, strict-overlap"]
    if orig is not None:
        lines += [
            nc("neimsCosineOrig",       _stat(orig, "NEIMS",  "Cosine")[0]),
            nc("qcxmsCosineOrig",       _stat(orig, "QCxMS",  "Cosine")[0]),
            nc("qcxmsTwoCosineOrig",    _stat(orig, "QCxMS2", "Cosine")[0]),
            nc("cfmidCosineOrig",       _stat(orig, "CFMID",  "Cosine")[0]),
        ]
    lines += [
        "%%",
        "%% Dataset sizes",
        "\\newcommand{\\nMolsFranklinTMS}{61}",
        "\\newcommand{\\nMolsUCB}{13}",
    ]

    content = "\n".join(lines) + "\n"
    dst = PAPER / "numbers.tex"
    if not dry_run:
        dst.write_text(content)
    print(f"  {'DRY-RUN' if dry_run else 'WRITTEN'}  numbers.tex")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--no-push", action="store_true")
    args = parser.parse_args()

    print(f"\n{'DRY RUN — ' if args.dry_run else ''}Paper sync {datetime.now():%Y-%m-%d %H:%M}\n")

    print("── Figures ──")
    for src_rel, dst_rel in FIGURE_MAP:
        print(sync_file(REPORTS / src_rel, PAPER / dst_rel, args.dry_run))

    print("\n── Tables ──")
    for src_rel, dst_rel in TABLE_MAP:
        print(sync_file(REPORTS / src_rel, PAPER / dst_rel, args.dry_run))

    print("\n── Numbers ──")
    generate_numbers_tex(args.dry_run)

    print("\n── Needs regeneration ──")
    for src_rel, note in MISSING_NOTES:
        if not (REPORTS / src_rel).exists():
            print(f"  TODO: {note}")

    if args.dry_run:
        print("\nDry run complete — no files changed.\n")
        return

    # Git commit and push
    print("\n── Git ──")
    result = subprocess.run(
        ["git", "status", "--short"],
        cwd=PAPER, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    changed = result.stdout.decode().strip()
    if not changed:
        print("  Nothing to commit.")
        return

    subprocess.run(["git", "add", "-A"], cwd=PAPER, check=True)
    msg = f"sync: update figures, tables, numbers ({datetime.now():%Y-%m-%d %H:%M})"
    subprocess.run(["git", "commit", "-m", msg], cwd=PAPER, check=True)

    if not args.no_push:
        subprocess.run(["git", "push"], cwd=PAPER, check=True)
        print("  Pushed to GitHub.")
    else:
        print("  Committed (skipped push).")


if __name__ == "__main__":
    main()

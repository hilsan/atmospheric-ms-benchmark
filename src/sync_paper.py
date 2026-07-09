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


def generate_numbers_tex(dry_run: bool):
    """Write numbers.tex from known table values."""
    content = """\
%% numbers.tex — auto-generated by src/sync_paper.py
%% Do NOT edit manually. Re-run sync_paper.py to update.
%% Include in manuscript preamble: \\input{numbers}
%%
%% GoAmazon (Franklin TMS) — all peaks
\\newcommand{\\neimsCosineTMS}{812.1}
\\newcommand{\\neimsCosineTMSstd}{70.8}
\\newcommand{\\neimsTanimotoTMS}{0.7}
\\newcommand{\\neimsRecallTMS}{85.9}
\\newcommand{\\neimsPrecisionTMS}{76.8}
\\newcommand{\\neimsCosineTopTwentyTMS}{840.0}
%%
\\newcommand{\\qcxmsCosineTMS}{566.4}
\\newcommand{\\qcxmsTanimotoTMS}{0.5}
\\newcommand{\\qcxmsRecallTMS}{92.1}
\\newcommand{\\qcxmsPrecisionTMS}{53.9}
%%
\\newcommand{\\qcxmsTwoCosineTMS}{486.4}
\\newcommand{\\qcxmsTwoTanimotoTMS}{0.4}
\\newcommand{\\qcxmsTwoRecallTMS}{77.8}
\\newcommand{\\qcxmsTwoPrecisionTMS}{47.9}
%%
\\newcommand{\\cfmidCosineTMS}{528.0}
\\newcommand{\\cfmidTanimotoTMS}{0.4}
\\newcommand{\\cfmidRecallTMS}{57.9}
\\newcommand{\\cfmidPrecisionTMS}{58.5}
%%
%% External dataset (UCB-GLOBES) — all peaks
\\newcommand{\\neimsCosinUCB}{685.0}
\\newcommand{\\neimsTanimotoUCB}{0.6}
\\newcommand{\\neimsRecallUCB}{72.0}
\\newcommand{\\neimsPrecisionUCB}{87.0}
%%
\\newcommand{\\qcxmsTwoCosineUCB}{358.0}
\\newcommand{\\qcxmsTwoTanimotoUCB}{0.5}
\\newcommand{\\qcxmsTwoRecallUCB}{77.1}
\\newcommand{\\qcxmsTwoPrecisionUCB}{61.5}
%%
%% Non-derivatized (Franklin) — all peaks, strict-overlap N=31
\\newcommand{\\neimsCosineOrig}{787.7}
\\newcommand{\\qcxmsCosineOrig}{557.7}
\\newcommand{\\qcxmsTwoCosineOrig}{445.3}
\\newcommand{\\cfmidCosineOrig}{516.3}
%%
%% Dataset sizes
\\newcommand{\\nMolsFranklinTMS}{61}
\\newcommand{\\nMolsUCB}{13}
"""
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

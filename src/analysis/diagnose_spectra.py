"""
Diagnostics for spectra completeness across methods/datasets.

Reports per-method coverage, lists missing molecules, and cross-method
overlap to distinguish problematic molecules from method-specific failures.
"""

from pathlib import Path

import numpy as np
import pandas as pd

EXPECTED_SPECTRA = ["spectra_all.csv", "spectra_10pct.csv", "spectra_top20.csv"]


def _spectra_complete(spectra_dir):
    return all((Path(spectra_dir) / f).exists() for f in EXPECTED_SPECTRA)


def diagnose_spectra(sim_base_dir, datasets, n_mols=None, plot=True, output_dir=None):
    """Check spectra completeness for all methods and report overlaps.

    Parameters
    ----------
    sim_base_dir : str
        Root directory containing one subfolder per method/dataset.
    datasets : list[str]
        Method/dataset folder names to check (e.g. ["QCxMS_25_ps", "NEIMS", "EXP"]).
    n_mols : int, optional
        Expected number of molecules. Auto-detected from the largest folder
        index found if not provided.
    plot : bool
        If True, show a presence/absence heatmap.

    Returns
    -------
    dict[str, dict]
        Per-dataset dicts with keys "complete" and "missing" (sets of folder names).
    """
    sim_base_dir = Path(sim_base_dir)

    # Auto-detect n_mols
    if n_mols is None:
        for dataset in datasets:
            dp = sim_base_dir / dataset
            if dp.exists():
                indices = [
                    int(d.name) for d in dp.iterdir()
                    if d.is_dir() and d.name.isdigit()
                ]
                if indices:
                    n_mols = max(indices) + 1
                    break
    if n_mols is None:
        raise ValueError("Could not auto-detect n_mols — no dataset directories found.")

    all_mols = [f"{i:04d}" for i in range(n_mols)]

    # ── Collect per-method status ────────────────────────────────────────────
    results = {}
    for dataset in datasets:
        dp = sim_base_dir / dataset
        complete, missing, no_dir = [], [], []
        for mol in all_mols:
            mol_path = dp / mol
            if not mol_path.exists():
                no_dir.append(mol)
                missing.append(mol)
            elif _spectra_complete(mol_path / "spectra"):
                complete.append(mol)
            else:
                missing.append(mol)
        results[dataset] = {
            "complete": set(complete),
            "missing":  set(missing),
            "no_dir":   set(no_dir),
        }

    # ── Summary table ────────────────────────────────────────────────────────
    print(f"\n{'═'*65}")
    print(f"  Spectra diagnostics — {sim_base_dir.name}")
    print(f"{'═'*65}")
    print(f"  {'Method':<28} {'Complete':>8} {'Missing':>8} {'No dir':>7}  {'%':>6}")
    print(f"  {'─'*28} {'─'*8} {'─'*8} {'─'*7}  {'─'*6}")
    for dataset in datasets:
        r    = results[dataset]
        n_ok = len(r["complete"])
        n_mi = len(r["missing"])
        n_nd = len(r["no_dir"])
        pct  = 100 * n_ok / n_mols if n_mols else 0
        flag = "  ✓" if n_mi == 0 else ""
        print(f"  {dataset:<28} {n_ok:>8} {n_mi:>8} {n_nd:>7}  {pct:>5.1f}%{flag}")
    print(f"  {'─'*28} {'─'*8} {'─'*8} {'─'*7}  {'─'*6}")
    print(f"  Total molecules expected: {n_mols}")

    # ── Per-method missing list ───────────────────────────────────────────────
    any_missing = any(results[d]["missing"] for d in datasets)
    if any_missing:
        print(f"\n{'─'*65}")
        print("  Missing molecules per method:")
        print(f"{'─'*65}")
        for dataset in datasets:
            missing = sorted(results[dataset]["missing"])
            if missing:
                print(f"  {dataset:<28}: {', '.join(missing)}")

    # ── Cross-method overlap ─────────────────────────────────────────────────
    sim_datasets = [d for d in datasets if d != "EXP"]
    if len(sim_datasets) > 1:
        # Count how many methods are missing each molecule
        miss_count = {
            mol: sum(1 for d in sim_datasets if mol in results[d]["missing"])
            for mol in all_mols
        }
        miss_count = {mol: n for mol, n in miss_count.items() if n > 0}

        if miss_count:
            print(f"\n{'─'*65}")
            print("  Cross-method overlap (simulation methods only):")
            print(f"{'─'*65}")

            threshold = max(2, len(sim_datasets) // 2)
            problematic    = sorted(m for m, n in miss_count.items() if n >= threshold)
            method_specific = sorted(m for m, n in miss_count.items() if n == 1)

            if problematic:
                print(f"\n  Missing in ≥{threshold}/{len(sim_datasets)} methods "
                      f"(likely problematic molecules):")
                for mol in problematic:
                    methods = [d for d in sim_datasets if mol in results[d]["missing"]]
                    print(f"    {mol}  ← {', '.join(methods)}")

            if method_specific:
                print(f"\n  Missing in exactly 1 method (method-specific failures):")
                for mol in method_specific:
                    method = next(d for d in sim_datasets if mol in results[d]["missing"])
                    print(f"    {mol}  ← {method}")

            in_between = sorted(
                m for m, n in miss_count.items()
                if 1 < n < threshold
            )
            if in_between:
                print(f"\n  Missing in 2–{threshold-1} methods:")
                for mol in in_between:
                    methods = [d for d in sim_datasets if mol in results[d]["missing"]]
                    print(f"    {mol}  ← {', '.join(methods)}")

    # ── Heatmap ──────────────────────────────────────────────────────────────
    if plot:
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns

            matrix = pd.DataFrame(
                {
                    d: [1 if mol in results[d]["complete"] else 0 for mol in all_mols]
                    for d in datasets
                },
                index=all_mols,
            )

            # Only show rows with at least one missing entry
            rows_with_missing = matrix[matrix.sum(axis=1) < len(datasets)]

            if rows_with_missing.empty:
                print("\n  All molecules complete across all methods — no heatmap needed.")
            else:
                fig, ax = plt.subplots(
                    figsize=(max(6, len(datasets) * 1.1), max(4, len(rows_with_missing) * 0.35 + 1))
                )
                sns.heatmap(
                    rows_with_missing,
                    cmap=["#D9534F", "#5A99D3"],
                    linewidths=0.4,
                    linecolor="white",
                    vmin=0, vmax=1,
                    ax=ax,
                    cbar=False,
                    yticklabels=True,
                )
                ax.set_title(
                    "Spectra completeness — molecules with ≥1 missing method\n"
                    "  blue = complete   red = missing",
                    fontsize=12, pad=10,
                )
                ax.set_xlabel("")
                ax.tick_params(axis="x", rotation=35, labelsize=10)
                ax.tick_params(axis="y", labelsize=9)
                plt.tight_layout()
                plt.show()
        except ImportError:
            print("\n  (install matplotlib + seaborn for heatmap visualisation)")

    # ── LaTeX output ─────────────────────────────────────────────────────────
    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        dataset_label = Path(sim_base_dir).name.replace("_", r"\_")

        # ── Coverage table ────────────────────────────────────────────────────
        cov_rows = []
        for dataset in datasets:
            r    = results[dataset]
            n_ok = len(r["complete"])
            n_mi = len(r["missing"])
            pct  = 100 * n_ok / n_mols if n_mols else 0
            method_tex = dataset.replace("_", r"\_")
            cov_rows.append(
                f"  {method_tex} & {n_ok} & {n_mi} & {pct:.1f}\\% \\\\"
            )

        cov_tex = "\n".join([
            "% Requires: \\usepackage{booktabs}",
            "\\begin{table}[ht]",
            "\\centering",
            f"\\caption{{Spectra completeness per method ({dataset_label}). "
            f"\\textit{{Complete}}: all three peak-type spectra present. "
            f"\\textit{{Missing}}: molecule absent or incomplete.}}",
            f"\\label{{tab:diagnostics_coverage_{Path(sim_base_dir).name}}}",
            "\\begin{tabular}{lccc}",
            "\\toprule",
            "Method & Complete & Missing & \\% complete \\\\",
            "\\midrule",
            *cov_rows,
            "\\midrule",
            f"  Total expected & {n_mols} & & \\\\",
            "\\bottomrule",
            "\\end{tabular}",
            "\\end{table}",
        ])
        cov_path = Path(output_dir) / "diagnostics_coverage.tex"
        cov_path.write_text(cov_tex)
        print(f"Saved: {cov_path}")

        # ── Missing-molecules table ───────────────────────────────────────────
        miss_rows = []
        for dataset in datasets:
            missing = sorted(results[dataset]["missing"])
            if missing:
                method_tex  = dataset.replace("_", r"\_")
                missing_tex = ", ".join(missing)
                miss_rows.append(f"  {method_tex} & {missing_tex} \\\\")

        if miss_rows:
            miss_tex = "\n".join([
                "% Requires: \\usepackage{booktabs}",
                "\\begin{table}[ht]",
                "\\centering",
                f"\\caption{{Missing molecules per method ({dataset_label}). "
                "Molecule IDs refer to the zero-padded four-digit index in the benchmark set.}",
                f"\\label{{tab:diagnostics_missing_{Path(sim_base_dir).name}}}",
                "\\begin{tabular}{lp{10cm}}",
                "\\toprule",
                "Method & Missing molecule IDs \\\\",
                "\\midrule",
                *miss_rows,
                "\\bottomrule",
                "\\end{tabular}",
                "\\end{table}",
            ])
            miss_path = Path(output_dir) / "diagnostics_missing.tex"
            miss_path.write_text(miss_tex)
            print(f"Saved: {miss_path}")

    return results

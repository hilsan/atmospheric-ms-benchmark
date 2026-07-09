"""
Visualization functions for the atmospheric-MS-benchmark results.

Provides histograms, heatmaps, and summary tables comparing multiple
simulation methods against experimental spectra.

Usage
-----
from src.visualization.plot_results import run_all, PALETTE_FULL, PALETTE_PRESENTATION

# Exploration (all QCxMS variants)
run_all("../data/simulation_results/franklin/", output_dir="plots/franklin",
        palette=PALETTE_FULL)

# Presentation (QCxMS_25_ps renamed to QCxMS, 10_ps dropped)
run_all("../data/simulation_results/franklin_tms/", output_dir="plots/franklin_tms_pres",
        palette=PALETTE_PRESENTATION, presentation=True)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

# ---------------------------------------------------------------------------
# Palettes
# ---------------------------------------------------------------------------
PALETTE_FULL = {
    "QCxMS":              "#5A99D3",   # canonical variant (renamed by load_results)
    "QCxMS_10_ps":        "#F0913A",   # kept for variant-comparison cell
    "QCxMS_10_ps_iee03":  "#E05AB1",
    "QCxMS_25_ps":        "#D36EA5",
    "NEIMS":              "#4C4C4C",
    "CFMID":              "#A6A6A6",
    "QCxMS2":             "#7B68EE",
    "QCxMS2_dft":         "#B088EE",
}

PALETTE_PRESENTATION = {
    "QCxMS":  "#5A99D3",
    "QCxMS2": "#4C4C4C",
    "QCxMS2_dft": "#888888",
    "CFMID":  "#A6A6A6",
    "NEIMS":  "#D36EA5",
}

# ---------------------------------------------------------------------------
# Shared config
# ---------------------------------------------------------------------------
METRICS = ["Cosine", "Weighted_Dot", "Tanimoto", "%S_sim_in_ref", "%R_ref_in_sim", "Entropy_Similarity"]

METRIC_LABELS = {
    "Cosine":             "Cosine",
    "Weighted_Dot":       "Weighted Dot",
    "Tanimoto":           "Tanimoto",
    "%S_sim_in_ref":      "% sim. peak in ref",
    "%R_ref_in_sim":      "% ref. peak in sim",
    "Entropy_Similarity": "Entropy similarity",
}

PEAK_LABELS = {
    "all":   "All",
    "top20": "Top 20",
    "10pct": "≥10%",
}

PEAK_ORDER = ["all", "top20", "10pct"]

CMAP = LinearSegmentedColormap.from_list(
    "custom", ["#A6A6A6", "#5A99D3", "#D36EA5", "#4C4C4C"], N=256
)

MIN_MOLECULES = 20


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_results(base_dir, presentation=False, qcxms_canonical=None, tms_mismatch=None):
    """Load all comparison CSVs from the results subfolder.

    Parameters
    ----------
    base_dir : str
        Root directory containing a ``results/`` subfolder.
    presentation : bool
        If True, keep only ``qcxms_canonical`` among QCxMS variants and
        rename it to "QCxMS".
    qcxms_canonical : str, optional
        The QCxMS variant to keep in presentation mode (e.g. "QCxMS_25_ps").
        Defaults to "QCxMS_25_ps" if not specified.
    tms_mismatch : dict, optional
        {"mol_indices": [...], "fallback_sim_dir": "..."} — for molecules where
        the automated TMS derivatisation assigned TMS groups that don't match the
        reference (TMS=0), substitute scores from an underivatised benchmark run.
        The raw TMS simulation results are kept on disk; only the in-memory data
        used for analysis is replaced.
    """
    _canonical = qcxms_canonical or "QCxMS_25_ps"
    _all_qcxms = {"QCxMS_10_ps", "QCxMS_25_ps", "QCxMS_10_ps_iee03"}

    def _read_results_dir(rdir):
        dfs = []
        for f in sorted(Path(rdir).glob("*/*.csv")):
            if f.parent.name.isdigit() and not f.stem.startswith("spectra_all_comparison_bin"):
                df = pd.read_csv(f)
                df["Molecule"]  = f.parent.name
                df["Peak_Type"] = f.stem.split("_")[1]
                dfs.append(df)
        return dfs

    dfs = _read_results_dir(Path(base_dir) / "results")
    data = pd.concat(dfs, ignore_index=True)

    if tms_mismatch:
        fb_base = Path(tms_mismatch["fallback_sim_dir"]).resolve() / "results"
        for mol_idx in tms_mismatch["mol_indices"]:
            fb_dfs = _read_results_dir(fb_base)
            fb_mol = [d for d in fb_dfs if d["Molecule"].iloc[0] == mol_idx] if fb_dfs else []
            if fb_mol:
                data = data[data["Molecule"] != mol_idx]
                data = pd.concat([data] + fb_mol, ignore_index=True)
                print(f"  TMS mismatch: substituted Franklin results for mol {mol_idx}")

    # Drop non-canonical QCxMS variants; rename canonical to "QCxMS" in all modes
    drop = _all_qcxms - {_canonical}
    data = data[~data["Method"].isin(drop)]
    data["Method"] = data["Method"].replace({_canonical: "QCxMS"})

    # Fill any metric columns missing from old CSVs (e.g. Entropy_Similarity) with NaN
    for col in METRICS:
        if col not in data.columns:
            data[col] = np.nan

    print(f"Loaded {data.shape[0]} rows from {len(dfs)} files — {base_dir}")
    return data


def get_active_methods(data, palette):
    """Return methods with data for at least MIN_MOLECULES molecules.

    Ordered by palette, with unknown methods appended at the end.
    """
    mol_counts = (
        data.dropna(subset=METRICS, how="all")
        .groupby("Method")["Molecule"]
        .nunique()
    )
    n_total = data["Molecule"].nunique()
    threshold = min(MIN_MOLECULES, max(1, n_total // 3))
    active = mol_counts[mol_counts >= threshold].index.tolist()
    ordered = [m for m in palette if m in active]
    ordered += [m for m in active if m not in ordered]

    excluded = mol_counts[mol_counts < threshold]
    if not excluded.empty:
        print(f"Methods excluded (< {threshold} molecules):")
        for m, n in excluded.items():
            print(f"  {m}: {n}")

    print(f"Active methods: {ordered}")
    return ordered


# ---------------------------------------------------------------------------
# Histogram plots
# ---------------------------------------------------------------------------
def plot_histograms(data, palette, output_dir=None, presentation=False, paper=False):
    """Bar-histogram of each similarity metric for each peak type."""
    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)

    method_order = get_active_methods(data, palette)
    palette_spaced = {m.replace("_", " "): v for m, v in palette.items() if m in method_order}
    method_order_spaced = [m.replace("_", " ") for m in method_order]
    plot_data = data.copy()
    plot_data["Method"] = plot_data["Method"].str.replace("_", " ")

    if presentation:
        sns.set_style("white")
        sns.set_context("paper")
        fig_w, fig_h, dpi = 12.0, 6.5, 300
        title_fs, label_fs, tick_fs, leg_fs, leg_title_fs = 14, 14, 13, 13, 14
    elif paper:
        sns.set_style("white")
        sns.set_context("paper")
        fig_w, fig_h, dpi = 7.08, 4.0, 300
        title_fs, label_fs, tick_fs, leg_fs, leg_title_fs = 9, 8, 8, 7, 8
    else:
        sns.set_style("white")
        sns.set_context("notebook", font_scale=1.5)
        fig_w, fig_h, dpi = 10, 6, 150
        title_fs, label_fs, tick_fs, leg_fs, leg_title_fs = 20, 16, 14, 12, 14

    BINS = 15
    n_methods = len(method_order_spaced)

    for pt in sorted(data["Peak_Type"].unique()):
        subset = plot_data[
            (plot_data["Peak_Type"] == pt)
            & (plot_data["Method"].isin(method_order_spaced))
        ].dropna(subset=METRICS)
        pt_label = PEAK_LABELS.get(pt, pt)

        for metric in METRICS:
            metric_label = METRIC_LABELS[metric]
            fig, ax = plt.subplots(figsize=(fig_w, fig_h))
            if presentation:
                fig.subplots_adjust(top=0.78)

            all_vals = subset[metric].dropna()
            bin_edges = np.linspace(all_vals.min(), all_vals.max(), BINS + 1)
            bin_width = bin_edges[1] - bin_edges[0]
            bar_width = bin_width / n_methods * 0.9

            for i, method in enumerate(method_order_spaced):
                vals = subset[subset["Method"] == method][metric].dropna()
                if vals.empty:
                    continue
                counts, _ = np.histogram(vals, bins=bin_edges, density=True)
                offsets = bin_edges[:-1] + i * bar_width
                ax.bar(offsets, counts, width=bar_width,
                       color=palette_spaced[method], alpha=0.85, label=method)

            ax.set_xlabel(metric_label, fontsize=label_fs)
            ax.set_ylabel("Density", fontsize=label_fs)
            ax.tick_params(labelsize=tick_fs)

            if presentation:
                leg = ax.legend(title="Method", fontsize=leg_fs,
                                title_fontsize=leg_title_fs, frameon=False,
                                ncol=2, loc="lower center",
                                bbox_to_anchor=(0.5, 1.01))
            elif paper:
                leg = ax.legend(title="Method", fontsize=leg_fs,
                                title_fontsize=leg_title_fs, frameon=False,
                                ncol=1, loc="upper left")
            else:
                leg = ax.legend(title="Method", fontsize=leg_fs,
                                title_fontsize=leg_title_fs, frameon=True,
                                framealpha=0.9)

            sns.despine(trim=False, ax=ax)
            plt.tight_layout()

            if output_dir:
                fname = Path(output_dir) / f"{metric_label} {pt_label}.png"
                plt.savefig(fname, dpi=dpi, bbox_inches="tight",
                            bbox_extra_artists=(leg,))
                print(f"Saved: {fname}")

            plt.show()
            plt.close()


def plot_peak_type_panels(data, palette, output_dir=None, presentation=False, paper=False):
    """Single figure with one panel per peak type showing cosine distributions.

    Panels: All peaks | ≥10% | Top 20 — one grouped histogram per panel.
    Designed for presentation mode; also works for paper/explore.
    """
    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)

    method_order = get_active_methods(data, palette)
    palette_spaced = {m.replace("_", " "): v for m, v in palette.items() if m in method_order}
    method_order_spaced = [m.replace("_", " ") for m in method_order]
    plot_data = data.copy()
    plot_data["Method"] = plot_data["Method"].str.replace("_", " ")

    if presentation:
        fig_w, fig_h, dpi = 13.0, 5.0, 300
        title_fs, label_fs, tick_fs, leg_fs, leg_title_fs = 14, 14, 13, 13, 14
    elif paper:
        fig_w, fig_h, dpi = 7.08, 2.5, 300
        title_fs, label_fs, tick_fs, leg_fs, leg_title_fs = 8, 8, 7, 7, 8
    else:
        fig_w, fig_h, dpi = 12.0, 4.0, 150
        title_fs, label_fs, tick_fs, leg_fs, leg_title_fs = 12, 11, 10, 10, 11

    BINS = 15
    n_methods = len(method_order_spaced)
    pt_order = [p for p in PEAK_ORDER if p in data["Peak_Type"].unique()]

    fig, axes = plt.subplots(1, len(pt_order), figsize=(fig_w, fig_h),
                              sharey=False, facecolor="none")
    if len(pt_order) == 1:
        axes = [axes]

    for ax, pt in zip(axes, pt_order):
        ax.set_facecolor("none")
        subset = plot_data[
            (plot_data["Peak_Type"] == pt)
            & (plot_data["Method"].isin(method_order_spaced))
        ].dropna(subset=["Cosine"])

        all_vals = subset["Cosine"].dropna()
        if all_vals.empty:
            ax.set_visible(False)
            continue

        bin_edges = np.linspace(all_vals.min(), all_vals.max(), BINS + 1)
        bin_width = bin_edges[1] - bin_edges[0]
        bar_width = bin_width / n_methods * 0.9

        for i, method in enumerate(method_order_spaced):
            vals = subset[subset["Method"] == method]["Cosine"].dropna()
            if vals.empty:
                continue
            counts, _ = np.histogram(vals, bins=bin_edges, density=True)
            offsets = bin_edges[:-1] + i * bar_width
            ax.bar(offsets, counts, width=bar_width,
                   color=palette_spaced[method], alpha=0.85, label=method)

        ax.set_title(PEAK_LABELS.get(pt, pt), fontsize=title_fs)
        ax.set_xlabel("Cosine score", fontsize=label_fs)
        ax.set_ylabel("Density" if pt == pt_order[0] else "", fontsize=label_fs)
        ax.tick_params(labelsize=tick_fs)
        sns.despine(ax=ax)

    # Single legend on the last panel, placed below the axes
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=palette_spaced[m], alpha=0.85)
        for m in method_order_spaced if m in palette_spaced
    ]
    fig.legend(handles, method_order_spaced,
               title="Method", fontsize=leg_fs, title_fontsize=leg_title_fs,
               frameon=False, ncol=n_methods,
               loc="upper center", bbox_to_anchor=(0.5, 0.0))

    plt.tight_layout()
    if output_dir:
        fname = Path(output_dir) / "cosine_by_peak_type.png"
        plt.savefig(fname, dpi=dpi, bbox_inches="tight")
        plt.savefig(Path(output_dir) / "cosine_by_peak_type.pdf", bbox_inches="tight")
        print(f"Saved: {fname}")
    plt.show()
    plt.close()


# ---------------------------------------------------------------------------
# Summary tables
# ---------------------------------------------------------------------------
def print_summary_tables(data, palette, output_dir=None, presentation=False):
    """Print mean ± std tables per peak type for active methods.

    If output_dir is given, also saves a .tex file per peak type.
    When presentation=True, Weighted Dot is omitted.
    """
    from IPython.display import display

    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)

    metrics = [m for m in METRICS if not (presentation and m == "Weighted_Dot")]
    method_order = get_active_methods(data, palette)

    for pt in sorted(data["Peak_Type"].unique()):
        subset = data[data["Peak_Type"] == pt].dropna(subset=metrics)
        rows, latex_rows = [], []
        for method in method_order:
            method_data = subset[subset["Method"] == method]
            row = {"Method": method}
            latex_row = {"Method": method.replace("_", r"\_")}
            for metric in metrics:
                vals = method_data[metric].dropna()
                if len(vals) > 0:
                    row[metric] = f"{vals.mean():.1f} ± {vals.std():.1f}"
                    latex_row[metric] = f"{vals.mean():.1f} $\\pm$ {vals.std():.1f}"
                else:
                    row[metric] = latex_row[metric] = "N/A"
            rows.append(row)
            latex_rows.append(latex_row)

        if not rows:
            print(f"\n=== Peak Type: {pt} === (no active methods)")
            continue
        df_table = pd.DataFrame(rows).set_index("Method")
        df_table.columns = [METRIC_LABELS[m] for m in df_table.columns]
        print(f"\n=== Peak Type: {pt} ===")
        display(df_table)

        if output_dir:
            df_latex = pd.DataFrame(latex_rows).set_index("Method")
            df_latex.columns = [METRIC_LABELS[m] for m in df_latex.columns]
            tex = df_latex.to_latex(escape=False)
            tex_path = Path(output_dir) / f"table_{pt}.tex"
            tex_path.write_text(tex)
            print(f"  Saved LaTeX: {tex_path}")


# ---------------------------------------------------------------------------
# Heatmaps
# ---------------------------------------------------------------------------
def plot_heatmaps(data, palette, output_dir=None, presentation=False, paper=False):
    """Heatmap of mean scores per method × peak type for each metric."""
    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)

    method_order = get_active_methods(data, palette)

    if presentation:
        sns.set_style("white")
        sns.set_context("paper")
        fig_w, fig_h, dpi = 12.0, 6.5, 300
        title_fs, label_fs, tick_fs, annot_fs, cbar_fs = 14, 14, 13, 13, 13
    elif paper:
        sns.set_style("white")
        sns.set_context("paper")
        fig_w, fig_h, dpi = 7.08, 3.5, 300
        title_fs, label_fs, tick_fs, annot_fs, cbar_fs = 9, 8, 8, 7, 7
    else:
        sns.set_style("white")
        sns.set_context("notebook", font_scale=1.5)
        fig_w, fig_h, dpi = 8, 5, 150
        title_fs, label_fs, tick_fs, annot_fs, cbar_fs = 16, 14, 14, 14, 12

    for metric, metric_label in METRIC_LABELS.items():
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))

        mean_matrix = pd.DataFrame(
            index=[m.replace("_", " ") for m in method_order],
            columns=[PEAK_LABELS[p] for p in PEAK_ORDER],
        )
        annot_matrix = mean_matrix.copy()

        for method in method_order:
            for pt in PEAK_ORDER:
                vals = data[
                    (data["Method"] == method) & (data["Peak_Type"] == pt)
                ][metric].dropna()
                mean_val = vals.mean() if len(vals) > 0 else np.nan
                mean_matrix.loc[method.replace("_", " "), PEAK_LABELS[pt]] = mean_val
                if presentation or paper:
                    annot_matrix.loc[method.replace("_", " "), PEAK_LABELS[pt]] = (
                        f"{mean_val:.1f}" if not np.isnan(mean_val) else "N/A"
                    )
                else:
                    annot_matrix.loc[method.replace("_", " "), PEAK_LABELS[pt]] = (
                        f"{mean_val:.1f}\n±{vals.std():.1f}"
                        if not np.isnan(mean_val) else "N/A"
                    )

        mean_matrix = mean_matrix.astype(float)

        sns.heatmap(
            mean_matrix,
            annot=annot_matrix,
            fmt="",
            cmap=CMAP,
            linewidths=0.5,
            linecolor="white",
            ax=ax,
            cbar_kws={"label": metric_label},
            annot_kws={"size": annot_fs, "color": "white"},
        )

        ax.set_xlabel("Peak Strategy", fontsize=label_fs)
        ax.set_ylabel("Method", fontsize=label_fs)
        ax.tick_params(labelsize=tick_fs)
        sns.despine(trim=False)

        cbar = ax.collections[0].colorbar
        cbar.ax.tick_params(labelsize=cbar_fs)
        cbar.set_label(metric_label, fontsize=cbar_fs)

        plt.tight_layout()

        if output_dir:
            fname = Path(output_dir) / f"Heatmap {metric_label}.png"
            plt.savefig(fname, dpi=dpi, bbox_inches="tight")
            print(f"Saved: {fname}")

        plt.show()
        plt.close()


# ---------------------------------------------------------------------------
# Master function
# ---------------------------------------------------------------------------
def run_all(base_dir, output_dir=None, palette=None, presentation=False,
            paper=False, qcxms_canonical=None, tms_mismatch=None):
    """Load results and generate all plots and tables.

    Parameters
    ----------
    base_dir : str
        Root directory containing a ``results/`` subfolder.
    output_dir : str, optional
        Directory to save plots. Plots are only shown if None.
    palette : dict, optional
        Colour map {method_name: hex_colour}. Defaults to PALETTE_FULL.
    presentation : bool
        Slide-optimised figures (7" wide, 16–24 pt fonts). Drops non-canonical
        QCxMS variants and renames the canonical one to "QCxMS".
    paper : bool
        Paper-optimised figures (ACS 1-col = 3.35", 8 pt fonts, 300 DPI).
        Like presentation, keeps only the canonical QCxMS variant.
    qcxms_canonical : str, optional
        QCxMS variant to use as canonical (default: "QCxMS_25_ps").
    tms_mismatch : dict, optional
        Passed to load_results — substitutes underivatised benchmark scores for
        molecules where automated TMS assignment disagrees with the reference.
    """
    paper_or_pres = presentation or paper
    if palette is None:
        palette = PALETTE_PRESENTATION if paper_or_pres else PALETTE_FULL

    data = load_results(base_dir, presentation=paper_or_pres,
                        qcxms_canonical=qcxms_canonical,
                        tms_mismatch=tms_mismatch)
    print_summary_tables(data, palette, output_dir=output_dir,
                         presentation=paper_or_pres)
    plot_histograms(data, palette, output_dir=output_dir,
                    presentation=presentation, paper=paper)
    plot_peak_type_panels(data, palette, output_dir=output_dir,
                          presentation=presentation, paper=paper)
    plot_heatmaps(data, palette, output_dir=output_dir,
                  presentation=presentation, paper=paper)
    return data


# ---------------------------------------------------------------------------
# Cross-dataset UMAP / t-SNE
# ---------------------------------------------------------------------------

_ALL_QCXMS = {"QCxMS_10_ps", "QCxMS_25_ps", "QCxMS_10_ps_iee03"}

# Colours used inside each UMAP panel to distinguish datasets (edge colour)
_DATASET_EDGE_COLORS = ["#222222", "#D36EA5", "#5A99D3", "#4C4C4C"]


def plot_cross_dataset_umap(
    datasets,
    output_dir,
    qcxms_canonical="QCxMS_10_ps",
    methods=None,
    figname="umap_cross_dataset",
    presentation=False,
    paper=False,
):
    """UMAP / t-SNE molecular-space plot pooling multiple datasets.

    Panels correspond to simulation methods; points are coloured by cosine
    score and shaped by dataset.

    Parameters
    ----------
    datasets : list of dict, each with keys:
        ``label``         display name shown in the legend
        ``base_dir``      path to simulation results root (contains results/)
        ``processed_csv`` path to processed molecule CSV
        ``smiles_col``    column name with SMILES used for Morgan fingerprints
        ``marker``        matplotlib marker string, e.g. "o", "^", "s"
    output_dir : str or Path
    qcxms_canonical : str
        QCxMS variant to rename to "QCxMS".
    methods : list[str] or None
        Methods to show as panels. None → ["QCxMS", "QCxMS2", "CFMID", "NEIMS"].
    figname : str
        Output file stem (no extension).
    presentation : bool
        Slide-optimised sizing (16 pt fonts, larger markers).
    paper : bool
        Publication sizing. Default width = ACS 2-column (7.08").
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D

    _PY38 = "/users/hsandstr/NEIMS/QCxMS2/qcxms2/bin/python3"

    def _umap_via_subprocess(X_arr, n_neighbors, min_dist, seed):
        """Run UMAP in Python 3.8 env (where umap-learn works) via subprocess."""
        import subprocess, tempfile, os
        with tempfile.TemporaryDirectory() as tmpdir:
            x_path = os.path.join(tmpdir, "X.npy")
            e_path = os.path.join(tmpdir, "emb.npy")
            np.save(x_path, X_arr)
            script = (
                f"import numpy as np; from umap import UMAP; "
                f"X=np.load('{x_path}'); "
                f"emb=UMAP(n_neighbors={n_neighbors},min_dist={min_dist},"
                f"random_state={seed}).fit_transform(X); "
                f"np.save('{e_path}',emb)"
            )
            r = subprocess.run([_PY38, "-c", script],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               timeout=300)
            if r.returncode != 0:
                raise RuntimeError(f"UMAP subprocess failed:\n{r.stderr.decode()[-2000:]}")
            return np.load(e_path)

    try:
        import umap as _umap_mod
        _reducer_name = "UMAP"
        _use_subprocess = False
    except Exception:
        try:
            import subprocess as _sp
            _sp.run([_PY38, "-c", "from umap import UMAP"],
                    check=True,
                    stdout=_sp.PIPE, stderr=_sp.PIPE,
                    timeout=120)
            _reducer_name = "UMAP"
            _use_subprocess = True
        except Exception:
            from sklearn.manifold import TSNE as _TSNE
            _reducer_name = "t-SNE"
            _use_subprocess = False

    _methods = methods or ["QCxMS", "QCxMS2", "CFMID", "NEIMS"]
    n = len(_methods)
    ncols = min(n, 2)
    nrows = (n + ncols - 1) // ncols

    # --- sizing presets ---
    if paper:
        col_w = 7.08 / ncols        # 2-col ACS width shared across columns
        panel_h = col_w * 0.95
        fs = dict(title=9, label=8, tick=7, cbar=7, legend=7)
        pt_size, lw = 35, 0.3
    elif presentation:
        col_w, panel_h = 5.5, 5.0
        fs = dict(title=18, label=15, tick=14, cbar=14, legend=14)
        pt_size, lw = 80, 0.5
    else:
        col_w, panel_h = 5.0, 4.5
        fs = dict(title=15, label=13, tick=11, cbar=11, legend=10)
        pt_size, lw = 60, 0.4

    # --- build fingerprint matrix ---
    all_fps, all_ds_idx, all_mol_ids = [], [], []
    for di, ds in enumerate(datasets):
        df = pd.read_csv(ds["processed_csv"]).reset_index(drop=True)
        scol = ds.get("smiles_col", "SMILES")
        for i, row in df.iterrows():
            smi = str(row.get(scol, "") or "")
            mol = Chem.MolFromSmiles(smi) if smi else None
            if mol is None:
                continue
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
            all_fps.append(list(fp))
            all_ds_idx.append(di)
            all_mol_ids.append((di, f"{i:04d}"))

    X = np.array(all_fps)
    all_ds_idx = np.array(all_ds_idx)
    print(f"Fingerprints: {X.shape} ({_reducer_name})")

    # --- dimensionality reduction ---
    if _reducer_name == "UMAP" and _use_subprocess:
        emb = _umap_via_subprocess(X, n_neighbors=min(10, len(X) - 1),
                                   min_dist=0.2, seed=42)
    elif _reducer_name == "UMAP":
        reducer = _umap_mod.UMAP(
            n_neighbors=min(10, len(X) - 1), min_dist=0.2, random_state=42
        )
        emb = reducer.fit_transform(X)
    else:
        reducer = _TSNE(
            n_components=2, perplexity=min(15, len(X) - 1), random_state=42
        )
        emb = reducer.fit_transform(X)

    # --- clip outliers in embedding space (IQR × 3) so one extreme mol doesn't
    #     collapse all others into a line ---
    for dim in range(emb.shape[1]):
        q1, q3 = np.percentile(emb[:, dim], [25, 75])
        iqr = q3 - q1
        lo, hi = q1 - 3 * iqr, q3 + 3 * iqr
        emb[:, dim] = np.clip(emb[:, dim], lo, hi)

    # --- load cosine scores ---
    scores = {m: {di: {} for di in range(len(datasets))} for m in _methods}
    for di, ds in enumerate(datasets):
        rdir = Path(ds["base_dir"]) / "results"
        if not rdir.exists():
            continue
        for mol_dir in sorted(rdir.iterdir()):
            if not mol_dir.name.isdigit():
                continue
            csv = mol_dir / "spectra_all_comparison.csv"
            if not csv.exists():
                continue
            df_r = pd.read_csv(csv).set_index("Method")
            if qcxms_canonical in df_r.index:
                df_r = df_r.rename(index={qcxms_canonical: "QCxMS"})
            df_r = df_r.drop(
                index=[m for m in df_r.index if m in _ALL_QCXMS],
                errors="ignore",
            )
            for m in _methods:
                if m in df_r.index:
                    val = df_r.loc[m, "Cosine"]
                    scores[m][di][mol_dir.name] = float(val) if pd.notna(val) else np.nan

    # --- plot ---
    norm = Normalize(vmin=0, vmax=1000)
    cmap = "viridis"

    fig, axes = plt.subplots(
        nrows, ncols, figsize=(col_w * ncols, panel_h * nrows), squeeze=False
    )

    for pidx, method in enumerate(_methods):
        ax = axes[pidx // ncols, pidx % ncols]
        c_all = np.array([
            scores[method][di].get(mid, np.nan) for di, mid in all_mol_ids
        ])
        valid = np.isfinite(c_all)

        for di, ds in enumerate(datasets):
            ds_mask = all_ds_idx == di
            miss = ds_mask & ~valid
            hit  = ds_mask & valid
            ec = _DATASET_EDGE_COLORS[di % len(_DATASET_EDGE_COLORS)]
            if miss.any():
                ax.scatter(
                    emb[miss, 0], emb[miss, 1],
                    c="lightgrey", s=pt_size * 0.7,
                    marker=ds["marker"],
                    edgecolors=ec, linewidths=lw * 0.5, zorder=1,
                )
            if hit.any():
                ax.scatter(
                    emb[hit, 0], emb[hit, 1],
                    c=c_all[hit], cmap=cmap, norm=norm,
                    s=pt_size, marker=ds["marker"],
                    edgecolors=ec, linewidths=lw, zorder=2,
                )

        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.85)
        cbar.set_label("Cosine score", fontsize=fs["cbar"])
        cbar.ax.tick_params(labelsize=fs["tick"])

        ax.set_title(method, fontsize=fs["title"], pad=6)
        ax.set_xlabel(f"{_reducer_name} 1", fontsize=fs["label"])
        ax.set_ylabel(f"{_reducer_name} 2", fontsize=fs["label"])
        ax.tick_params(labelsize=fs["tick"])

    # hide unused panels
    for pidx in range(n, nrows * ncols):
        axes[pidx // ncols, pidx % ncols].set_visible(False)

    # dataset legend on the first axis
    legend_handles = [
        Line2D(
            [0], [0], marker=ds["marker"], color="w",
            markerfacecolor="#999999",
            markeredgecolor=_DATASET_EDGE_COLORS[di % len(_DATASET_EDGE_COLORS)],
            markersize=7 if paper else 10,
            label=ds["label"],
        )
        for di, ds in enumerate(datasets)
    ]
    axes[0, 0].legend(
        handles=legend_handles, fontsize=fs["legend"],
        framealpha=0.85, loc="best",
    )

    plt.tight_layout()

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        plt.savefig(out / f"{figname}.{ext}", dpi=300, bbox_inches="tight")
    print(f"Saved: {out / figname}.png  (.pdf)")
    plt.show()
    plt.close()

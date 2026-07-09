"""
Molecular property analysis and visualisation for the atmospheric MS benchmark.

Main entry point:
    run_analysis(ds, paper=False, presentation=False)

where `ds` is one of the dataset dicts defined in notebook 6 (or equivalent).
"""

import os

import matplotlib.font_manager as fm
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors as RDDescriptors
    _RDKIT = True
except ImportError:
    _RDKIT = False

GREY_TEXT = '#4C4C4C'

METHOD_PALETTE = {
    'QCxMS':  '#5A99D3',
    'QCxMS2': '#4C4C4C',
    'CFMID':  '#A6A6A6',
    'NEIMS':  '#D36EA5',
}

QCXMS_CANONICAL = 'QCxMS_10_ps'

C1, C2, C3 = '#5A99D3', '#D36EA5', '#7BC8A4'

_inter_reg = False
for _fp in ['/tmp/inter_fonts/Inter-Regular.ttf', '/tmp/inter_fonts/Inter-Bold.ttf']:
    if Path(_fp).exists():
        fm.fontManager.addfont(_fp)
        _inter_reg = True
_PAPER_FONT = ['Inter', 'Nimbus Sans', 'DejaVu Sans'] if _inter_reg else ['Nimbus Sans', 'DejaVu Sans']


def run_analysis(ds, paper=False, presentation=False):
    DATASET          = ds['name']
    SIM_ROOT         = os.path.abspath(ds['sim_dir'])
    INCLUDE_TMS_PLOT = ds.get('include_tms_plot', False)
    tms_col          = ds.get('tms_col', 'TMS')
    _col_hist        = ds.get('hist_color', '#5A99D3')

    if presentation:
        OUT_DIR = Path(ds['pres_dir'])
    elif paper:
        OUT_DIR = Path(ds['paper_dir'])
    else:
        OUT_DIR = Path(ds['reports_dir'])
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"{DATASET}  |  mode: {'pres' if presentation else 'paper' if paper else 'explore'}")
    print(f"{'='*60}")

    # ── Load ─────────────────────────────────────────────────────────────────
    _simpol_path = Path(ds['simpol_csv'])
    if not _simpol_path.exists():
        raise FileNotFoundError(
            f"Underivatized SIMPOL CSV not found for {DATASET}: {_simpol_path}\n"
            f"Run cell 3 (SIMPOL generation) first."
        )
    simpol  = pd.read_csv(_simpol_path).reset_index(drop=True)
    dataset = pd.read_csv(ds['dataset_csv']).reset_index(drop=True)

    _simpol_tms_path = ds.get('simpol_tms_csv')
    if _simpol_tms_path:
        _simpol_tms_path = Path(_simpol_tms_path)
        if not _simpol_tms_path.exists():
            raise FileNotFoundError(
                f"TMS SIMPOL CSV not found for {DATASET}: {_simpol_tms_path}\n"
                f"Run cell 3 (SIMPOL generation) first."
            )
        simpol_tms = pd.read_csv(_simpol_tms_path).reset_index(drop=True)
    else:
        simpol_tms = simpol

    _deriv_path = ds.get('deriv_csv')
    deriv = pd.read_csv(_deriv_path).reset_index(drop=True) if _deriv_path else dataset

    idx_rows = []
    _qcxms_dir = Path(SIM_ROOT, QCXMS_CANONICAL)
    if _qcxms_dir.exists():
        for mol_dir in sorted(_qcxms_dir.iterdir()):
            if mol_dir.name.isdigit():
                idx_rows.append({'mol_idx': mol_dir.name})
    idx_df = pd.DataFrame(idx_rows).reset_index(drop=True)

    props = idx_df.copy()
    props['TMS'] = dataset[tms_col].values[:len(idx_df)] if tms_col in dataset.columns else 0
    for col in simpol.columns:
        if col != 'SMILES':
            props[col] = simpol[col].values[:len(idx_df)]
    print(f'Properties: {props.shape}')

    score_rows = []
    _results_dir = Path(SIM_ROOT, 'results')
    if _results_dir.exists():
        for mol_dir in sorted(_results_dir.iterdir()):
            if mol_dir.name.isdigit():
                csv = mol_dir / 'spectra_all_comparison.csv'
                if csv.exists():
                    df = pd.read_csv(csv)
                    df['mol_idx'] = mol_dir.name
                    score_rows.append(df)

    if score_rows:
        scores_long = pd.concat(score_rows, ignore_index=True)
        scores_long['Method'] = scores_long['Method'].replace({QCXMS_CANONICAL: 'QCxMS'})
        scores_long = scores_long[scores_long['Method'].isin(METHOD_PALETTE)]

        _tms_mm = ds.get('tms_mismatch')
        if _tms_mm:
            _fb_base = Path(os.path.abspath(_tms_mm['fallback_sim_dir'])) / 'results'
            for _mol_idx in _tms_mm['mol_indices']:
                _fb_csv = _fb_base / _mol_idx / 'spectra_all_comparison.csv'
                if _fb_csv.exists():
                    _fb_df = pd.read_csv(_fb_csv)
                    _fb_df['mol_idx'] = _mol_idx
                    _fb_df['Method'] = _fb_df['Method'].replace({QCXMS_CANONICAL: 'QCxMS'})
                    _fb_df = _fb_df[_fb_df['Method'].isin(METHOD_PALETTE)]
                    scores_long = scores_long[scores_long['mol_idx'] != _mol_idx]
                    scores_long = pd.concat([scores_long, _fb_df], ignore_index=True)
                    print(f'  TMS mismatch: substituted Franklin results for mol {_mol_idx}')
                else:
                    print(f'  TMS mismatch: fallback CSV not found for mol {_mol_idx} — keeping TMS results')

        scores_wide = scores_long.pivot(index='mol_idx', columns='Method', values='Cosine').reset_index()
        scores_wide.columns.name = None
        scores_wide = scores_wide.rename(columns={m: f'Cosine_{m}' for m in METHOD_PALETTE})
    else:
        scores_wide = pd.DataFrame(columns=['mol_idx'])

    master = props.merge(scores_wide, on='mol_idx', how='left')
    print(f'Master table: {master.shape}')

    _props_corr = idx_df.copy()
    for _c in simpol_tms.columns:
        if _c != 'SMILES':
            _props_corr[_c] = simpol_tms[_c].values[:len(idx_df)]
    master_corr = _props_corr.merge(scores_wide, on='mol_idx', how='left')

    # ── 60% completion filter ─────────────────────────────────────────────────
    n_total = len(props)
    active_palette = {}
    for m in METHOD_PALETTE:
        col = f'Cosine_{m}'
        pct = master[col].notna().sum() / n_total if col in master.columns else 0.0
        if pct >= 0.60:
            active_palette[m] = METHOD_PALETTE[m]
        else:
            print(f'  Dropping {m}: {pct:.0%} complete (< 60% threshold)')
    if not active_palette:
        print(f'  No method meets 60% threshold — skipping plots for {DATASET}')
        return

    # ── Score vs carbon / oxygen / TMS count ─────────────────────────────────
    n_rows  = 3 if INCLUDE_TMS_PLOT else 2
    methods = list(active_palette)
    if presentation:
        _fw, _lfs, _tfs, _pt_s = 12.0, 14, 13, 40
    elif paper:
        _fw, _lfs, _tfs, _pt_s = 7.08, 8, 8, 12
    else:
        _fw, _lfs, _tfs, _pt_s = 7.0, 8, 8, 12

    fig, axes_arr = plt.subplots(n_rows, len(methods),
                                  figsize=(_fw, n_rows * (3.5 if presentation else 2.0)),
                                  facecolor='none', squeeze=False)
    axes = [list(axes_arr[r]) for r in range(n_rows)]
    if n_rows == 2:
        axes.append([None] * len(methods))

    for col_i, method in enumerate(methods):
        score_col = f'Cosine_{method}'
        label = method.replace('_', ' ')

        ax = axes[0][col_i]
        ax.set_facecolor('none')
        md = master.dropna(subset=['carbon number', score_col])
        ax.scatter(md['carbon number'], md[score_col],
                   color=active_palette[method], alpha=0.75, s=_pt_s, linewidths=0)
        r_str = f"r={stats.spearmanr(md['carbon number'], md[score_col])[0]:.2f}" if len(md) > 1 else "n/a"
        ax.set_xlabel(f"{label}\n# Carbon ({r_str})", fontsize=_lfs, color=GREY_TEXT)
        ax.tick_params(labelsize=_tfs, colors=GREY_TEXT)
        for spine in ax.spines.values(): spine.set_edgecolor('#cccccc')
        if col_i == 0: ax.set_ylabel('Cosine score', fontsize=_lfs, color=GREY_TEXT)

        ax = axes[1][col_i]
        ax.set_facecolor('none')
        md = master.dropna(subset=['oxygen_count', score_col])
        ax.scatter(md['oxygen_count'], md[score_col],
                   color=active_palette[method], alpha=0.75, s=_pt_s, linewidths=0)
        r_str = f"r={stats.spearmanr(md['oxygen_count'], md[score_col])[0]:.2f}" if len(md) > 1 else "n/a"
        ax.set_xlabel(f"{label}\n# Oxygen ({r_str})", fontsize=_lfs, color=GREY_TEXT)
        ax.tick_params(labelsize=_tfs, colors=GREY_TEXT)
        for spine in ax.spines.values(): spine.set_edgecolor('#cccccc')
        if col_i == 0: ax.set_ylabel('Cosine score', fontsize=_lfs, color=GREY_TEXT)

        if INCLUDE_TMS_PLOT:
            ax = axes[2][col_i]
            ax.set_facecolor('none')
            tms_data = master.dropna(subset=['TMS', score_col]).copy()
            if tms_data.empty:
                ax.set_visible(False)
                continue
            tms_data['TMS_group'] = tms_data['TMS'].astype(int).astype(str)
            tms_order = sorted(tms_data['TMS_group'].unique(), key=int)
            sns.boxplot(data=tms_data, x='TMS_group', y=score_col, order=tms_order,
                        color=active_palette[method], ax=ax,
                        flierprops=dict(marker='o', markersize=3, alpha=0.5),
                        linewidth=0.8, width=0.5)
            ax.set_xlabel(f"{label}\nTMS groups", fontsize=_lfs, color=GREY_TEXT)
            ax.tick_params(labelsize=_tfs, colors=GREY_TEXT)
            for spine in ax.spines.values(): spine.set_edgecolor('#cccccc')
            if col_i == 0: ax.set_ylabel('Cosine score', fontsize=_lfs, color=GREY_TEXT)

    plt.tight_layout()
    plt.savefig(OUT_DIR / 'score_vs_properties.png', dpi=300, bbox_inches='tight')
    plt.savefig(OUT_DIR / 'score_vs_properties.pdf', bbox_inches='tight')
    plt.show(); plt.close()

    # ── Functional group correlation heatmap (uses TMS-derivatised SIMPOL) ────
    _FG_LIST = [
        'carbon number', 'oxygen_count', 'hydroxyl (alkyl)', 'carbonyl',
        'carboxylic acid', 'ester, all', 'ether', 'ketone', 'aldehyde',
        'nitrate', 'aromatic hydroxyl', 'aromatic_ring', 'non_aromatic_ring',
        'C=C (non-aromatic)', 'nitro', 'organosulfate',
    ]
    FG_COLS_CORR = [fg for fg in _FG_LIST if fg in master_corr.columns and master_corr[fg].sum() > 0]
    FG_COLS      = [fg for fg in _FG_LIST if fg in master.columns      and master[fg].sum()      > 0]
    corr_rows = []
    for method in active_palette:
        col = f'Cosine_{method}'
        md = master_corr.dropna(subset=FG_COLS_CORR + [col])
        for fg in FG_COLS_CORR:
            r = np.nan if (md[fg].std() == 0 or len(md) < 3) else stats.spearmanr(md[fg], md[col])[0]
            corr_rows.append({'Method': method, 'Feature': fg, 'r': r})
    corr_df = pd.DataFrame(corr_rows).pivot(index='Feature', columns='Method', values='r')
    corr_df.index = (corr_df.index
        .str.replace('oxygen_count', '# Oxygen')
        .str.replace('carbon number', '# Carbon')
        .str.replace('_', ' '))
    corr_df = corr_df[list(active_palette)]

    if presentation:
        _hfw, _hfh_min, _hfh_row, _htfs, _annot = 9.0, 6.0, 0.55, 14, 13
    else:
        _hfw, _hfh_min, _hfh_row, _htfs, _annot = 3.35, 4.5, 0.35, 7, 6

    fig, ax = plt.subplots(figsize=(_hfw, max(_hfh_min, _hfh_row * len(FG_COLS_CORR))),
                           facecolor='none')
    sns.heatmap(corr_df, annot=True, fmt='.2f', cmap='RdBu_r',
                center=0, vmin=-0.6, vmax=0.6,
                linewidths=0.4, linecolor='white',
                annot_kws={'size': _annot}, ax=ax,
                cbar_kws={'label': 'Spearman r', 'shrink': 0.6})
    ax.set_xlabel(''); ax.set_ylabel('')
    ax.tick_params(labelsize=_htfs, colors=GREY_TEXT)
    ax.collections[0].colorbar.ax.tick_params(labelsize=_htfs)
    ax.collections[0].colorbar.set_label('Spearman r', fontsize=_htfs)
    plt.tight_layout()
    plt.savefig(OUT_DIR / 'fg_correlation_heatmap.png', dpi=300, bbox_inches='tight')
    plt.savefig(OUT_DIR / 'fg_correlation_heatmap.pdf', bbox_inches='tight')
    plt.show(); plt.close()

    # ── Functional group prevalence bar chart ────────────────────────────────
    FG_DISPLAY = {
        'carbon number':      '# Carbon',
        'oxygen_count':       '# Oxygen',
        'hydroxyl (alkyl)':   'Hydroxyl (alkyl)',
        'carbonyl':           'Carbonyl',
        'carboxylic acid':    'Carboxylic acid',
        'ester, all':         'Ester',
        'ether':              'Ether',
        'ketone':             'Ketone',
        'aldehyde':           'Aldehyde',
        'nitrate':            'Nitrate',
        'aromatic hydroxyl':  'Aromatic hydroxyl',
        'aromatic_ring':      'Aromatic ring',
        'non_aromatic_ring':  'Non-aromatic ring',
        'C=C (non-aromatic)': 'C=C (non-aromatic)',
        'nitro':              'Nitro',
        'organosulfate':      'Organosulfate',
    }
    fg_means = {fg: master[fg].mean() for fg in FG_COLS if fg in master.columns}
    fg_sorted = sorted(fg_means, key=fg_means.get, reverse=True)
    fg_labels = [FG_DISPLAY.get(fg, fg.replace('_', ' ')) for fg in fg_sorted]
    fg_values = [fg_means[fg] for fg in fg_sorted]

    if presentation:
        _pfs, _pfw, _row_h = 14, 7.0, 0.55
    else:
        _pfs, _pfw, _row_h = 7, 3.35, 0.28
    fig, ax = plt.subplots(figsize=(_pfw, max(4.0 if presentation else 2.5, _row_h * len(fg_sorted))),
                           facecolor='none')
    ax.set_facecolor('none')
    ax.barh(range(len(fg_sorted)), fg_values, color='#5A99D3', alpha=0.85)
    ax.set_yticks(range(len(fg_sorted)))
    ax.set_yticklabels(fg_labels, fontsize=_pfs)
    ax.invert_yaxis()
    ax.set_xlabel('Mean count per molecule', fontsize=_pfs, color=GREY_TEXT)
    ax.tick_params(labelsize=_pfs, colors=GREY_TEXT)
    for spine in ax.spines.values():
        spine.set_edgecolor('#cccccc')
    plt.tight_layout()
    plt.savefig(OUT_DIR / 'fg_prevalence.png', dpi=300, bbox_inches='tight')
    plt.savefig(OUT_DIR / 'fg_prevalence.pdf', bbox_inches='tight')
    plt.show(); plt.close()

    # ── TMS substitution count distribution ──────────────────────────────────
    if INCLUDE_TMS_PLOT:
        tms_counts = master['TMS'].dropna().astype(int)
        tms_vals   = sorted(tms_counts.unique())
        tms_freq   = [int((tms_counts == v).sum()) for v in tms_vals]
        if presentation:
            _tfs2, _tfw, _tfh = 14, 6.0, 4.0
        else:
            _tfs2, _tfw, _tfh = 7, 3.35, 2.2

        fig, ax = plt.subplots(figsize=(_tfw, _tfh), facecolor='none')
        ax.set_facecolor('none')
        ax.bar(tms_vals, tms_freq, color='#D36EA5', alpha=0.85, width=0.6)
        ax.set_xticks(tms_vals)
        ax.set_xlabel('Number of TMS groups', fontsize=_tfs2, color=GREY_TEXT)
        ax.set_ylabel('Number of molecules', fontsize=_tfs2, color=GREY_TEXT)
        ax.tick_params(labelsize=_tfs2, colors=GREY_TEXT)
        for spine in ax.spines.values():
            spine.set_edgecolor('#cccccc')
        plt.tight_layout()
        plt.savefig(OUT_DIR / 'tms_distribution.png', dpi=300, bbox_inches='tight')
        plt.savefig(OUT_DIR / 'tms_distribution.pdf', bbox_inches='tight')
        plt.show(); plt.close()

    # ── Properties overview: 1×3 figure (MW / FG stacked / TMS stacked) ──────
    FG_OVERVIEW = {
        'non_aromatic_ring':  'Non-arom. ring',
        'carboxylic acid':    'Carboxylic acid',
        'C=C (non-aromatic)': 'C=C',
        'hydroxyl (alkyl)':   'Alkyl OH',
        'aromatic_ring':      'Aromatic ring',
        'ketone':             'Ketone',
        'ester, all':         'Ester',
        'aromatic hydroxyl':  'Aromatic OH',
        'aldehyde':           'Aldehyde',
        'ether':              'Ether',
    }

    _smiles_col = ds['smiles_col']
    if _RDKIT and _smiles_col in dataset.columns:
        def _mw(s):
            mol = Chem.MolFromSmiles(str(s)) if pd.notna(s) else None
            return RDDescriptors.ExactMolWt(mol) if mol else None
        mw_vals = dataset[_smiles_col].apply(_mw).dropna().values
    else:
        mw_vals = np.array([])

    _n = len(simpol)
    _fg_keys, _fg_labels, _fg_n1, _fg_n2, _fg_n3p = [], [], [], [], []
    for _k, _lbl in FG_OVERVIEW.items():
        if _k not in simpol.columns:
            continue
        _c = simpol[_k]
        if (_c > 0).sum() == 0:
            continue
        _fg_keys.append(_k); _fg_labels.append(_lbl)
        _fg_n1.append( 100 * (_c == 1).sum() / _n)
        _fg_n2.append( 100 * (_c == 2).sum() / _n)
        _fg_n3p.append(100 * (_c >= 3).sum() / _n)

    _C_OH    = C1; _C_COOH = C2; _C_BOTH = '#F0A040'; _C_OTHER = '#9AAFBA'; _C_NONE = '#D8DEE4'
    _tms_col_data = dataset[tms_col].reset_index(drop=True) if tms_col in dataset.columns else pd.Series(0, index=dataset.index)
    _tms_ints = sorted(_tms_col_data.dropna().astype(int).unique())
    _d_oh   = deriv['OH'].fillna(0).astype(int).reset_index(drop=True)   if 'OH'   in deriv.columns else pd.Series(0, index=range(len(deriv)))
    _d_cooh = deriv['COOH'].fillna(0).astype(int).reset_index(drop=True) if 'COOH' in deriv.columns else pd.Series(0, index=range(len(deriv)))
    _tms_n_oh, _tms_n_cooh, _tms_n_both, _tms_n_other, _tms_n_none = [], [], [], [], []
    for _t in _tms_ints:
        _mask = (_tms_col_data.fillna(-1).astype(int) == _t).values
        if _t == 0:
            _tms_n_oh.append(0); _tms_n_cooh.append(0); _tms_n_both.append(0)
            _tms_n_other.append(0); _tms_n_none.append(int(_mask.sum()))
        else:
            _oh_t = _d_oh.values[_mask]; _cooh_t = _d_cooh.values[_mask]
            _tms_n_none.append(0)
            _tms_n_oh.append(  int(((_oh_t > 0) & (_cooh_t == 0)).sum()))
            _tms_n_cooh.append(int(((_cooh_t > 0) & (_oh_t == 0)).sum()))
            _tms_n_both.append(int(((_oh_t > 0) & (_cooh_t > 0)).sum()))
            _tms_n_other.append(int(((_oh_t == 0) & (_cooh_t == 0)).sum()))

    _old_rc = {k: plt.rcParams[k] for k in
               ['font.family', 'font.sans-serif', 'axes.spines.top', 'axes.spines.right']}
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': _PAPER_FONT,
                         'axes.spines.top': False, 'axes.spines.right': False})

    if presentation:
        _fig_w, _fig_h, _ovfs = 13.0, 6.0, 14
    elif paper:
        _fig_w, _fig_h, _ovfs = 7.08, 3.0, 7
    else:
        _fig_w, _fig_h, _ovfs = 9.5, 3.8, 7

    _n_fg  = len(_fg_keys)
    _fig   = plt.figure(figsize=(_fig_w, _fig_h), constrained_layout=True)
    _gs    = gridspec.GridSpec(1, 3, figure=_fig, width_ratios=[1.0, 1.6, 1.0])
    _ax_a  = _fig.add_subplot(_gs[0])
    _ax_b  = _fig.add_subplot(_gs[1])
    _ax_c  = _fig.add_subplot(_gs[2])

    def _style(_ax, left=True):
        _ax.tick_params(colors=GREY_TEXT, length=2, labelsize=_ovfs)
        _ax.spines['left'].set_color('#cccccc' if left else 'none')
        _ax.spines['left'].set_visible(left)
        _ax.spines['bottom'].set_color('#cccccc')
        _ax.xaxis.label.set_color(GREY_TEXT)
        _ax.yaxis.label.set_color(GREY_TEXT)

    def _plabel(_ax, _t):
        _ax.text(-0.14, 1.05, _t, transform=_ax.transAxes,
                 fontsize=max(14, _ovfs + 3), fontweight='bold',
                 color=GREY_TEXT, va='bottom', ha='left')

    if len(mw_vals):
        _bins = np.linspace(50, 560, 22)
        _ax_a.hist(mw_vals, bins=_bins, color=_col_hist, alpha=0.80, linewidth=0)
    _ax_a.set_xlabel('Molecular weight (Da)', labelpad=3, fontsize=_ovfs)
    _ax_a.set_ylabel('Number of compounds', labelpad=3, fontsize=_ovfs)
    _ax_a.text(0.97, 0.95, f'n = {_n}', transform=_ax_a.transAxes,
               ha='right', va='top', fontsize=_ovfs, color=GREY_TEXT)
    _style(_ax_a)

    _y     = np.arange(_n_fg)
    _left2 = np.array(_fg_n1) + np.array(_fg_n2)
    _ax_b.barh(_y, _fg_n1,  height=0.5, color=C1, linewidth=0, label='1')
    _ax_b.barh(_y, _fg_n2,  height=0.5, color=C2, linewidth=0, left=_fg_n1, label='2')
    _ax_b.barh(_y, _fg_n3p, height=0.5, color=C3, linewidth=0, left=_left2, label='3+')
    _ax_b.set_yticks(_y)
    _ax_b.set_yticklabels(_fg_labels, fontsize=_ovfs, color=GREY_TEXT)
    _ax_b.invert_yaxis()
    _ax_b.set_xlabel('% of molecules', labelpad=2, fontsize=_ovfs)
    _ax_b.set_xlim(0, max(np.array(_fg_n1)+np.array(_fg_n2)+np.array(_fg_n3p)) * 1.15 if _fg_n1 else 10)
    _ax_b.tick_params(colors=GREY_TEXT, length=2, labelsize=_ovfs)
    _ax_b.tick_params(axis='y', length=0)
    _style(_ax_b, left=False)
    _ax_b.legend(title='Count', frameon=False, fontsize=_ovfs, title_fontsize=_ovfs,
                 labelcolor=GREY_TEXT, loc='upper center', bbox_to_anchor=(0.5, -0.18),
                 ncol=3, handlelength=1, handletextpad=0.4)

    _x  = np.array(_tms_ints)
    _b1 = np.array(_tms_n_oh); _b2 = _b1 + np.array(_tms_n_cooh)
    _b3 = _b2 + np.array(_tms_n_both); _b4 = _b3 + np.array(_tms_n_other)
    _ax_c.bar(_x, _tms_n_oh,    color=_C_OH,    linewidth=0, label='OH',    width=0.55)
    _ax_c.bar(_x, _tms_n_cooh,  color=_C_COOH,  linewidth=0, label='COOH',  width=0.55, bottom=_b1)
    _ax_c.bar(_x, _tms_n_both,  color=_C_BOTH,  linewidth=0, label='Both',  width=0.55, bottom=_b2)
    _ax_c.bar(_x, _tms_n_other, color=_C_OTHER, linewidth=0, label='Other', width=0.55, bottom=_b3)
    _ax_c.bar(_x, _tms_n_none,  color=_C_NONE,  linewidth=0, label='None',  width=0.55, bottom=_b4)
    _ax_c.set_xticks(_x)
    _ax_c.set_xlabel('Number of TMS groups', labelpad=2, fontsize=_ovfs)
    _ax_c.set_ylabel('Number of molecules', labelpad=3, fontsize=_ovfs)
    _style(_ax_c)
    _cat_totals = {'OH': sum(_tms_n_oh), 'COOH': sum(_tms_n_cooh), 'Both': sum(_tms_n_both),
                   'Other': sum(_tms_n_other), 'None': sum(_tms_n_none)}
    _hh, _ll = _ax_c.get_legend_handles_labels()
    _filtered = [(h, l) for h, l in zip(_hh, _ll) if _cat_totals.get(l, 0) > 0]
    if _filtered:
        _fh, _fl = zip(*_filtered)
        _ax_c.legend(_fh, _fl, title='Sub. type', frameon=False, fontsize=_ovfs, title_fontsize=_ovfs,
                     labelcolor=GREY_TEXT, loc='upper center', bbox_to_anchor=(0.5, -0.18),
                     ncol=5, handlelength=1, handletextpad=0.4)

    for _ax, _lbl in [(_ax_a,'a'), (_ax_b,'b'), (_ax_c,'c')]:
        _plabel(_ax, _lbl)

    plt.savefig(OUT_DIR / 'properties_figure.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.savefig(OUT_DIR / 'properties_figure.pdf',       bbox_inches='tight', pad_inches=0.05)
    plt.show(); plt.close()
    plt.rcParams.update(_old_rc)

    print(f'Done — outputs in {OUT_DIR}')

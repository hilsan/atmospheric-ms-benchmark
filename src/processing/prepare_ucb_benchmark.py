#!/usr/bin/env python3
"""
prepare_ucb_benchmark.py

Process a UCB-GLOBES tracer Excel file into a benchmark dataset:
  - De-duplicates entries by formula + pairwise cosine similarity
  - Collapses stereoisomer pairs (same formula, 0.90 ≤ cosine < 0.98) to one entry
  - Excludes ambiguous isomers (same formula, cosine < 0.90), mixtures, uncertain structures
  - Looks up parent SMILES via PubChem CIR and validates against Excel MW
  - Derivatizes Format-B compounds (common name + N TMS) via make_TMS_derivative
  - Writes EXP spectra folders (bin/scaled from Excel X/Y-Values)
  - Writes SI LaTeX table (reports/<dataset_name>/si_compounds_<dataset>.tex)
  - Writes benchmark CSV (data/processed/<dataset_name>/<dataset_name>.csv)

Usage:
    python prepare_ucb_benchmark.py \\
        -i data/raw/UCB_globes_tracers/UCB_GLOBES_ID_NotInNIST_toHilda_2026.04.16.xlsx \\
        -o data/simulation_results/ucb_globes_tracers \\
        --processed_dir data/processed/ucb_globes_tracers \\
        --reports_dir reports/ucb_globes_tracers \\
        --deriv_script src/processing/make_TMS_derivative_251125_v1.py
"""

import argparse
import re
import subprocess
import tempfile
import time
import urllib.parse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import requests

# ── Formula normalisation (Hill order) ───────────────────────────────────────
def _hill(formula):
    """Canonical Hill-order formula string, e.g. 'C18H32SiO3' → 'C18H32O3Si'."""
    import re
    counts = {}
    for elem, cnt in re.findall(r"([A-Z][a-z]?)(\d*)", str(formula)):
        if elem:
            counts[elem] = counts.get(elem, 0) + (int(cnt) if cnt else 1)
    order = (["C"] if "C" in counts else []) + (["H"] if "H" in counts else [])
    order += sorted(k for k in counts if k not in ("C", "H"))
    return "".join(f"{e}{counts[e]}" for e in order)


# ── Thresholds ────────────────────────────────────────────────────────────────
COSINE_DUPLICATE = 0.98   # same compound, different injection → keep max peaks
COSINE_STEREO    = 0.90   # likely stereoisomer → keep one entry
# cosine < COSINE_STEREO with same formula → ambiguous, exclude both

HARD_EXCLUDE = {
    "2- and 3-hydroxyglutaric acid mixture, 3 TMS",
    "isomer of 10-hydroxypinonic acid?",
}

# Compounds whose structure is unambiguously known even though a same-formula
# isomer exists in the database.  Overrides the ambiguous_isomer rule.
KNOWN_KEEP = {
    "MBTCA, 3 TMS",   # 3-methylbutane-1,2,3-tricarboxylic acid; isomer structure unknown
}

# ── Spectrum utilities ────────────────────────────────────────────────────────
def _parse_xy(row):
    x = np.array([float(v) for v in str(row["X-Values"]).strip("[]").split()])
    y = np.array([float(v) for v in str(row["Y-Values"]).strip("[]").split()])
    return x, y


def bin_scale(x, y):
    """Round m/z to int, sum collisions, normalise base peak to 999."""
    mz = np.round(x).astype(int)
    d = {}
    for m, i in zip(mz, y):
        d[m] = d.get(m, 0) + i
    mx = max(d.values()) if d else 1.0
    return {m: i * 999.0 / mx for m, i in d.items()}


def cosine_sim(s1, s2):
    if not s1 or not s2:
        return 0.0
    mzs = set(s1) | set(s2)
    v1 = np.array([s1.get(m, 0.0) for m in mzs])
    v2 = np.array([s2.get(m, 0.0) for m in mzs])
    d = np.linalg.norm(v1) * np.linalg.norm(v2)
    return float(np.dot(v1, v2) / d) if d > 0 else 0.0


def write_spectra(binned, out_dir):
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(sorted(binned.items()), columns=["mz", "intensity"])
    top20 = df.sort_values("intensity", ascending=False).head(20).sort_values("mz")
    p10   = df[df["intensity"] >= 0.1 * df["intensity"].max()]
    fmt   = "%.4f"
    df.to_csv(Path(out_dir) / "spectra_all.csv",   index=False, header=False, float_format=fmt)
    p10.to_csv(Path(out_dir) / "spectra_10pct.csv", index=False, header=False, float_format=fmt)
    top20.to_csv(Path(out_dir) / "spectra_top20.csv", index=False, header=False, float_format=fmt)


# ── Name parsing ──────────────────────────────────────────────────────────────
def parse_name(name):
    """Return (base_name, n_tms_expected, is_already_derivatized).

    Format A: IUPAC name starting with 'trimethylsilyl' → already derivatized.
    Format B: common name with ', N TMS' or ', TMS' suffix → look up parent.
    Format C: no TMS → underivatized compound.
    """
    name = name.strip()
    if re.match(r"trimethylsilyl|bis\(trimethylsilyl\)", name, re.IGNORECASE):
        return name, None, True

    base = name
    n_tms = 0

    m = re.search(r",\s*(\d+)\s*TMS\s*$", base, re.IGNORECASE)
    if m:
        n_tms = int(m.group(1))
        base  = base[: m.start()].strip()
    else:
        m = re.search(r",\s*TMS\s*$", base, re.IGNORECASE)
        if m:
            n_tms = 1
            base  = base[: m.start()].strip()

    # strip trailing isomer index  e.g. "_1", " 1", " 2"
    base = re.sub(r"[\s_]\d+\s*$", "", base).strip()

    return base, n_tms, False


# ── PubChem lookup ────────────────────────────────────────────────────────────
def pubchem_smiles(name, retries=2):
    """Return (isomeric_smiles, monoisotopic_mass) or (None, None)."""
    enc = urllib.parse.quote(name)
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{enc}"
        f"/property/IsomericSMILES,MonoisotopicMass/JSON"
    )
    for attempt in range(retries + 1):
        try:
            r = requests.get(url, timeout=15)
            if r.status_code == 200:
                props = r.json()["PropertyTable"]["Properties"][0]
                return props.get("IsomericSMILES"), props.get("MonoisotopicMass")
            if r.status_code == 404:
                return None, None
        except Exception:
            if attempt < retries:
                time.sleep(1)
    return None, None


# ── Curated SMILES (keyed by base_name.lower() after parse_name()) ────────────
# PubChem does not index informal atmospheric-chemistry names; SMILES verified
# against formula/MW from the UCB-GLOBES Excel file.
CURATED_SMILES = {
    # Format C — no TMS groups (aldehyde/ketone functional groups only)
    "4-((1s,2r)-3,3-dimethyl-2-(3-oxobutyl)cyclobutyl)pent-4-enal":
        "O=CCCC(=C)C1CC(C)(C)C1CCC(C)=O",   # beta-caryophyllene aldehyde; C15H24O2
    "4-((1s,2r)-3,3-dimethyl-2-(3-oxobutyl)cyclobutyl)-4-oxobutanal":
        "O=CCCC(=O)C1CC(C)(C)C1CCC(C)=O",   # beta-nocaryophyllone aldehyde; C15H24O2
    "pinonaldehyde":
        "O=CC1CC(C)(C)C1CCC(C)=O",           # C10H16O2; aldehyde not TMS-derivatised

    # Format B — parent SMILES (derivatize() adds TMS groups)
    "beta-caryophyllonic acid":
        "OC(=O)CCC(=C)C1CC(C)(C)C1CCC(C)=O",  # C15H24O3; 1 COOH → 1 TMS
    "beta-nocaryophyllonic acid":
        "OC(=O)CCC(=O)C1CC(C)(C)C1CCC(C)=O",  # C15H24O3; 1 COOH → 1 TMS
    "pinonic acid":
        "OC(=O)C1CC(C)(C)C1CCC(C)=O",          # C10H16O3; 1 COOH → 1 TMS
    "pinonic acid isomer":
        "OC(=O)C1CC(C)(C)C1CCC(C)=O",          # C10H16O3; structure uncertain — TMS will mismatch (expected 3)
    "pinic acid isomer":
        "OC(=O)C1CC(C)(C)C1C(=O)O",            # C9H14O4; 2 COOH → 2 TMS
    "tricarballylic acid":
        "OC(=O)CC(CC(=O)O)C(=O)O",             # C6H8O6; 3 COOH → 3 TMS
    "2-methyltetrol":
        "OCC(O)(C)C(O)CO",                     # C5H12O4; 4 OH → 4 TMS (isoprene SOA tracer)
    "3-hydroxy-4,4-dimethylglutaric acid":
        "OC(=O)CC(O)C(C)(C)C(=O)O",            # C7H12O5; 2 COOH + 1 OH → 3 TMS
    "2-methylglyceric acid":
        "OCC(O)(C)C(=O)O",                     # C4H8O4; 1 COOH + 2 OH → 3 TMS
    "mbo triol":
        "OCC(O)C(C)(C)O",                      # C5H12O3; 3 OH → 3 TMS (MBO SOA tracer)
    "mbtca":
        "OC(=O)CC(C)(C(=O)O)CC(=O)O",          # C8H12O6; 3 COOH → 3 TMS (alpha-pinene SOA tracer)
}


# ── TMS derivatization ────────────────────────────────────────────────────────
def derivatize(smiles, deriv_script):
    """Run make_TMS_derivative on a single SMILES; return (modified_smiles, replacements_dict)."""
    with tempfile.TemporaryDirectory() as tmp:
        inp = Path(tmp) / "in.csv"
        out = Path(tmp) / "out.csv"
        inp.write_text(f"SMILES\n{smiles}\n")
        result = subprocess.run(
            ["python", deriv_script, "-i", str(inp), "-o", str(out)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        if result.returncode != 0 or not out.exists():
            return None, {}
        df = pd.read_csv(out)
        if df.empty:
            return None, {}
        row = df.iloc[0]
        reps = {
            "OH": int(row.get("OH", 0)),
            "SH": int(row.get("SH", 0)),
            "Primary_Amine":   int(row.get("Primary_Amine", 0)),
            "Secondary_Amine": int(row.get("Secondary_Amine", 0)),
            "Imine": int(row.get("Imine", 0)),
            "OOH":   int(row.get("OOH", 0)),
            "COOH":  int(row.get("COOH", 0)),
        }
        return str(row["Modified_SMILES"]), reps


# ── De-duplication ────────────────────────────────────────────────────────────
def classify(df, spectra):
    """Return dict UID → status string."""
    status = {}

    # hard exclusions first
    for _, row in df.iterrows():
        uid = row["UID"]
        if row["Name"] in HARD_EXCLUDE:
            status[uid] = "excluded:mixture_or_uncertain"
        elif pd.isna(row.get("Formula", float("nan"))):
            status[uid] = "excluded:no_formula"
        else:
            status[uid] = "keep"

    # build peak-count lookup and normalise formulas
    peak_counts = {}
    for _, row in df.iterrows():
        try:
            peak_counts[row["UID"]] = int(row["Num Peaks"])
        except (TypeError, ValueError):
            peak_counts[row["UID"]] = 0

    df = df.copy()
    df["_formula_norm"] = df["Formula"].apply(
        lambda f: _hill(f) if not pd.isna(f) else "nan"
    )

    # within each same-formula group apply cosine-based rules
    for formula, grp in df.groupby("_formula_norm"):
        if formula == "nan" or len(grp) < 2:
            continue

        # sort by descending peak count so [0] is always the richest spectrum
        uids = sorted(grp["UID"].tolist(), key=lambda u: -peak_counts.get(u, 0))

        # pairwise cosines
        cs = {}
        for i, j in combinations(range(len(uids)), 2):
            u1, u2 = uids[i], uids[j]
            cs[(u1, u2)] = cosine_sim(spectra.get(u1, {}), spectra.get(u2, {}))

        def get_cs(u1, u2):
            return cs.get((u1, u2), cs.get((u2, u1), 0.0))

        # greedy merge: group UIDs with cosine ≥ COSINE_STEREO into same cluster
        clusters = []
        assigned = set()
        for u in uids:
            if u in assigned:
                continue
            cluster = [u]
            for u2 in uids:
                if u2 in assigned or u2 == u:
                    continue
                if get_cs(u, u2) >= COSINE_STEREO:
                    cluster.append(u2)
            clusters.append(cluster)
            assigned |= set(cluster)

        # classify within each cluster
        for cluster in clusters:
            keeper = cluster[0]  # already sorted by npeaks desc
            for u in cluster[1:]:
                c = get_cs(keeper, u)
                if status.get(keeper) != "keep":
                    continue
                if c >= COSINE_DUPLICATE:
                    status[u] = f"duplicate_of:{keeper}"
                else:
                    status[u] = f"stereo_of:{keeper}"

        # cross-cluster pairs with cosine < COSINE_STEREO → ambiguous isomers
        if len(clusters) > 1:
            for ci, cj in combinations(range(len(clusters)), 2):
                for u1 in clusters[ci]:
                    for u2 in clusters[cj]:
                        if get_cs(u1, u2) < COSINE_STEREO:
                            for u in (u1, u2):
                                if status.get(u) == "keep":
                                    status[u] = "ambiguous_isomer"

    # KNOWN_KEEP override: restore compounds whose structure is unambiguously known
    name_lookup = dict(zip(df["UID"], df["Name"]))
    for uid, s in status.items():
        if s == "ambiguous_isomer" and name_lookup.get(uid) in KNOWN_KEEP:
            status[uid] = "keep"

    return status


# ── Mirror plots ──────────────────────────────────────────────────────────────
def make_mirror_plots(df, status, spectra, out_dir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not available — skipping mirror plots")
        return

    Path(out_dir).mkdir(parents=True, exist_ok=True)
    seen = set()

    for uid, s in status.items():
        partner = None
        if s.startswith("stereo_of:"):
            partner = s.split(":", 1)[1]
        elif s == "ambiguous_isomer":
            row  = df[df["UID"] == uid].iloc[0]
            same = df[(df["Formula"] == row["Formula"]) & (df["UID"] != uid)]
            if not same.empty:
                partner = same.iloc[0]["UID"]

        if partner is None:
            continue
        key = tuple(sorted([uid, partner]))
        if key in seen:
            continue
        seen.add(key)

        s1 = spectra.get(uid, {})
        s2 = spectra.get(partner, {})
        cs = cosine_sim(s1, s2)

        all_mz = sorted(set(s1) | set(s2))
        name1  = df[df["UID"] == uid].iloc[0]["Name"][:50]
        name2  = df[df["UID"] == partner].iloc[0]["Name"][:50]

        kept1  = "KEPT"    if status.get(uid)     == "keep" else "EXCL"
        kept2  = "KEPT"    if status.get(partner) == "keep" else "EXCL"

        fig, ax = plt.subplots(figsize=(10, 4))
        for m in all_mz:
            if m in s1:
                ax.bar(m,  s1[m] / 999, width=0.8, color="#5A99D3", alpha=0.85)
            if m in s2:
                ax.bar(m, -s2[m] / 999, width=0.8, color="#D9534F", alpha=0.85)
        ax.axhline(0, color="black", lw=0.8)
        ax.set_xlabel("m/z")
        ax.set_ylabel("Rel. intensity")
        ax.set_title(f"Mirror plot  ·  cosine = {cs:.3f}", fontsize=11)
        ax.text(0.01, 0.97, f"↑ [{kept1}] {uid}  {name1}",
                transform=ax.transAxes, va="top", fontsize=8, color="#5A99D3")
        ax.text(0.01, 0.03, f"↓ [{kept2}] {partner}  {name2}",
                transform=ax.transAxes, va="bottom", fontsize=8, color="#D9534F")
        plt.tight_layout()
        fname = f"mirror_{uid}_{partner}.png"
        fig.savefig(Path(out_dir) / fname, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Mirror: {fname}  cosine={cs:.3f}")


# ── SI LaTeX table ────────────────────────────────────────────────────────────
def write_si_table(df, status, smiles_assigned, bench_extra_cols, output_dir, dataset_name):
    """Emit a longtable-style SI compound table (booktabs + \\smi{} macros).

    Excluded entries are italicised; footnote symbols explain each exclusion
    reason. Footnotes are placed after \\end{longtable} because threeparttable
    is incompatible with longtable page-breaking.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # assign benchmark indices to kept entries in original row order
    bench_idx = {}
    idx = 0
    for _, row in df.iterrows():
        if status.get(row["UID"]) == "keep":
            bench_idx[row["UID"]] = f"{idx:04d}"
            idx += 1

    MARKS = {
        "duplicate":        r"$^\dagger$",
        "stereo":           r"$^\ddagger$",
        "ambiguous_isomer": r"$^{*}$",
        "tms_mismatch":     r"$^{\P}$",
        "excluded":         r"$^\S$",
    }
    FOOTNOTES = {
        "duplicate":
            r"Duplicate entry (same compound; cosine $\geq$ 0.98). "
            r"Entry with most peaks retained.",
        "stereo":
            r"Stereoisomer pair (cosine $\geq$ 0.90). "
            r"Entry with most peaks retained; SMILES assignment identical.",
        "ambiguous_isomer":
            r"Excluded: isomer present in dataset with spectral cosine similarity $<$ 0.90 "
            r"to a paired entry. SMILES assignment ambiguous; omitted from benchmark.",
        "tms_mismatch":
            r"Excluded: automated TMS derivatisation count did not match the count stated "
            r"in the compound name, indicating the assigned parent SMILES is incorrect. "
            r"Omitted from benchmark.",
        "excluded":
            r"Excluded: mixture of compounds or uncertain structure.",
    }

    def mark_for(s):
        if s.startswith("duplicate"):        return MARKS["duplicate"]
        if s.startswith("stereo"):           return MARKS["stereo"]
        if s == "ambiguous_isomer":          return MARKS["ambiguous_isomer"]
        if s == "excluded:tms_mismatch":     return MARKS["tms_mismatch"]
        if s.startswith("excluded"):         return MARKS["excluded"]
        return ""

    import ast, re

    def escape(t):
        return str(t).replace("_", r"\_").replace("&", r"\&")

    def smi_esc(s):
        return (str(s).replace("\\", r"\textbackslash{}")
                      .replace("#", r"\#")
                      .replace("%", r"\%"))

    def parse_tms(name):
        m = re.search(r',?\s*(\d+)\s*TMS\b', name, re.IGNORECASE)
        if m:
            return int(m.group(1))
        if re.search(r'\bTMS\b', name, re.IGNORECASE):
            return 1
        # "trimethylsilyl" spelled out — count occurrences
        count = len(re.findall(r'\btrimethylsilyl\b', name, re.IGNORECASE))
        return count

    used = set()
    rows_tex = []
    for _, row in df.iterrows():
        uid     = row["UID"]
        s       = status.get(uid, "keep")
        mark    = mark_for(s)
        bidx    = bench_idx.get(uid, "---")
        name    = escape(row["Name"])
        dbidx   = str(row.get("Index", "---"))
        is_kept = (bidx != "---")

        extras  = bench_extra_cols.get(uid, {})
        tms_cnt = parse_tms(row["Name"])
        p_smi   = smi_esc(extras.get("Original_SMILES", "")) if is_kept else ""
        t_smi   = smi_esc(extras.get("Modified_SMILES", "")
                          or smiles_assigned.get(uid, "")) if is_kept else ""

        p_col = f"\\smi{{{p_smi}}}" if p_smi else "---"

        if is_kept:
            if tms_cnt == 0:
                rows_tex.append(
                    f"{bidx} & {name}{mark} & \\smi{{{p_smi}}} & 0 & {dbidx} & "
                    r"\multicolumn{2}{c}{\textit{not derivatized}} \\[2pt]"
                )
            else:
                t_col = f"\\smi{{{t_smi}}}" if t_smi else "---"
                rows_tex.append(
                    f"{bidx} & {name}{mark} & {p_col} & {tms_cnt} & {dbidx} & "
                    f"{t_col} \\\\[2pt]"
                )
        else:
            rows_tex.append(
                f"\\textit{{---}} & \\textit{{{name}{mark}}} & "
                f"\\textit{{---}} & \\textit{{{tms_cnt}}} & \\textit{{{dbidx}}} & "
                r"\textit{---} \\[2pt]"
            )

        key = ("duplicate"        if s.startswith("duplicate")       else
               "stereo"           if s.startswith("stereo")          else
               "ambiguous_isomer" if s == "ambiguous_isomer"         else
               "tms_mismatch"     if s == "excluded:tms_mismatch"    else
               "excluded"         if s.startswith("excluded")        else None)
        if key:
            used.add(key)

    fn_lines = []
    for k in ["duplicate", "stereo", "ambiguous_isomer", "tms_mismatch", "excluded"]:
        if k in used:
            fn_lines.append(f"{MARKS[k]}~{FOOTNOTES[k]}\\\\[2pt]")

    ds_label  = escape(dataset_name)
    label     = f"tab:si_ucb_{dataset_name}"
    col_spec  = r"C{0.6cm} L{3.8cm} L{4.2cm} C{0.6cm} C{1.2cm} L{4.2cm}"
    hdr       = (r"\textbf{\#} & \textbf{Name} & \textbf{SMILES} & \textbf{TMS} & "
                 r"\textbf{DB idx} & \textbf{TMS SMILES} \\")

    tex = "\n".join([
        r"% Requires: \usepackage{booktabs,longtable}",
        r"% Uses \smi{} command and C{}/L{} column types from preamble",
        f"\\begin{{longtable}}{{{col_spec}}}",
        f"\\caption{{UCB-GLOBES tracer compounds ({ds_label}). "
        r"Index~(\#) is the numerical identifier used throughout this work; "
        r"rows marked --- were excluded from the benchmark (see footnotes). "
        r"TMS gives the number of trimethylsilyl groups as stated in the compound name. "
        r"DB~idx is the integer UCB-GLOBES database index. "
        r"SMILES is the parent (underivatised) structure; TMS SMILES is the "
        r"derivatised structure used for simulation. "
        r"--- indicates no assignment was available. "
        r"\textit{Italic rows} denote excluded entries.}",
        f"\\label{{{label}}} \\\\",
        r"\toprule",
        hdr,
        r"\midrule",
        r"\endfirsthead",
        f"\\multicolumn{{6}}{{c}}{{\\textit{{Table \\ref{{{label}}} "
        r"continued from previous page}} \\[4pt]",
        r"\toprule",
        hdr,
        r"\midrule",
        r"\endhead",
        r"\midrule",
        r"\multicolumn{6}{r}{\textit{Continued on next page}} \\",
        r"\endfoot",
        r"\bottomrule",
        r"\endlastfoot",
        "",
        *rows_tex,
        "",
        r"\end{longtable}",
        r"{\footnotesize\noindent",
        "\n".join(fn_lines),
        r"}",
    ])
    out = Path(output_dir).parent / "paper" / dataset_name / f"si_compounds_{dataset_name}.tex"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(tex)
    print(f"Saved SI table: {out}")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description="Prepare UCB-GLOBES benchmark dataset.")
    ap.add_argument("-i",  "--input",         required=True)
    ap.add_argument("-o",  "--output_dir",    required=True,
                    help="Simulation results base dir (EXP/ written here)")
    ap.add_argument("--processed_dir",        required=True,
                    help="Dir for benchmark CSV (data/processed/<dataset>)")
    ap.add_argument("--reports_dir",          required=True,
                    help="Dir for SI table and mirror plots")
    ap.add_argument("--dataset_name",         default="ucb_globes_tracers")
    ap.add_argument("--deriv_script",         default=None,
                    help="Path to make_TMS_derivative_251125_v1.py")
    ap.add_argument("--skip_smiles",          action="store_true",
                    help="Skip PubChem lookup (for testing)")
    args = ap.parse_args()

    df = pd.read_excel(args.input, engine="openpyxl")
    print(f"Loaded {len(df)} entries")

    # build binned spectra cache
    spectra = {}
    for _, row in df.iterrows():
        try:
            x, y = _parse_xy(row)
            spectra[row["UID"]] = bin_scale(x, y)
        except Exception as e:
            spectra[row["UID"]] = {}
            print(f"  Warning: no spectrum for {row['UID']}: {e}")

    # classify
    print("\nClassifying entries ...")
    status = classify(df, spectra)
    n_keep  = sum(1 for s in status.values() if s == "keep")
    n_drop  = len(status) - n_keep
    print(f"  Kept: {n_keep}   Dropped/collapsed: {n_drop}")
    for uid, s in status.items():
        name = df[df["UID"] == uid]["Name"].values[0]
        flag = "" if s == "keep" else f"  [{s}]"
        print(f"  {uid:<14} {name[:50]}{flag}")

    # SMILES lookup + derivatization
    smiles_assigned  = {}   # UID → derivatized SMILES (for SI table)
    bench_extra_cols = {}   # UID → {Original_SMILES, Modified_SMILES, Total, ...}

    if not args.skip_smiles:
        print("\nLooking up SMILES via PubChem ...")
        for _, row in df.iterrows():
            uid = row["UID"]
            if status.get(uid) != "keep":
                continue

            name = str(row["Name"])
            base, n_tms_expected, is_derivatized = parse_name(name)
            mw_excel = (None if pd.isna(row.get("MW", float("nan")))
                        else float(row["MW"]))

            smiles = CURATED_SMILES.get(base.lower())
            pc_mw  = None
            if smiles is not None:
                print(f"  [{uid}] Curated SMILES for: {base!r}")
            else:
                smiles, pc_mw = pubchem_smiles(base)
                time.sleep(0.25)  # PubChem rate limit

            if smiles is None:
                print(f"  [{uid}] FAILED lookup: {base!r}")
                continue

            # MW validation (use average MW from PubChem monoisotopic as proxy)
            if mw_excel and pc_mw:
                diff = abs(float(pc_mw) - mw_excel)
                if diff > 2.0:
                    print(f"  [{uid}] MW mismatch Δ={diff:.1f} "
                          f"(PubChem={pc_mw:.2f}, Excel={mw_excel:.2f}) — {base!r}")

            if is_derivatized:
                # Format A: SMILES is already the TMS derivative
                smiles_assigned[uid] = smiles
                bench_extra_cols[uid] = {
                    "Original_SMILES": smiles,
                    "Modified_SMILES": smiles,
                    "Total_Replacements": n_tms_expected or 0,
                }
                print(f"  [{uid}] Format-A (already derivatised): {smiles[:60]}")
            elif n_tms_expected == 0 or n_tms_expected is None:
                # No TMS groups expected — use parent SMILES directly
                smiles_assigned[uid] = smiles
                bench_extra_cols[uid] = {
                    "Original_SMILES": smiles,
                    "Modified_SMILES": smiles,
                    "Total_Replacements": 0,
                }
                print(f"  [{uid}] Format-C (no TMS): {smiles[:60]}")
            else:
                # Format B: derivatize parent
                if args.deriv_script:
                    mod_smiles, reps = derivatize(smiles, args.deriv_script)
                    total = sum(reps.values())
                    if total != n_tms_expected:
                        print(f"  [{uid}] TMS count mismatch: got {total}, expected {n_tms_expected} — {name!r} — EXCLUDING")
                        status[uid] = "excluded:tms_mismatch"
                        # still record SMILES for SI table transparency
                        smiles_assigned[uid] = mod_smiles or smiles
                        bench_extra_cols[uid] = {
                            "Original_SMILES": smiles,
                            "Modified_SMILES": mod_smiles or smiles,
                            "Total_Replacements": total,
                            **reps,
                        }
                        continue
                    smiles_assigned[uid] = mod_smiles or smiles
                    bench_extra_cols[uid] = {
                        "Original_SMILES": smiles,
                        "Modified_SMILES": mod_smiles or smiles,
                        "Total_Replacements": total,
                        **reps,
                    }
                    print(f"  [{uid}] Format-B derivatised (TMS={total}/{n_tms_expected}): {(mod_smiles or smiles)[:60]}")
                else:
                    smiles_assigned[uid] = smiles
                    bench_extra_cols[uid] = {"Original_SMILES": smiles, "Modified_SMILES": smiles}
                    print(f"  [{uid}] Format-B parent (no deriv_script): {smiles[:60]}")

    # write EXP folders
    exp_dir = Path(args.output_dir) / "EXP"
    print(f"\nWriting EXP spectra → {exp_dir}")
    idx = 0
    bench_rows = []
    for _, row in df.iterrows():
        uid = row["UID"]
        if status.get(uid) != "keep":
            continue
        mol_idx  = f"{idx:04d}"
        out_path = exp_dir / mol_idx / "spectra"
        binned   = spectra.get(uid, {})
        if not binned:
            print(f"  [{uid}] No spectrum — skipping.")
            idx += 1
            continue
        write_spectra(binned, out_path)
        extras = bench_extra_cols.get(uid, {})
        bench_rows.append({
            "Index":            mol_idx,
            "UID":              uid,
            "Name":             row["Name"],
            "Contributor_ID":   row.get("Contributor_ID", ""),
            "DB_Index":         row.get("Index", ""),
            "Original_SMILES":  extras.get("Original_SMILES", ""),
            "Modified_SMILES":  extras.get("Modified_SMILES", ""),
            "Total_Replacements": extras.get("Total_Replacements", ""),
            "OH":               extras.get("OH", ""),
            "SH":               extras.get("SH", ""),
            "COOH":             extras.get("COOH", ""),
            "Primary_Amine":    extras.get("Primary_Amine", ""),
            "Secondary_Amine":  extras.get("Secondary_Amine", ""),
            "Imine":            extras.get("Imine", ""),
            "OOH":              extras.get("OOH", ""),
        })
        print(f"  [{mol_idx}] {uid} written")
        idx += 1

    # benchmark CSV
    Path(args.processed_dir).mkdir(parents=True, exist_ok=True)
    bench_csv = Path(args.processed_dir) / f"{args.dataset_name}.csv"
    pd.DataFrame(bench_rows).to_csv(bench_csv, index=False)
    print(f"\nBenchmark CSV ({len(bench_rows)} compounds): {bench_csv}")

    # mirror plots
    mirror_dir = Path(args.reports_dir) / "mirror_plots"
    print(f"\nMirror plots → {mirror_dir}")
    make_mirror_plots(df, status, spectra, mirror_dir)

    # SI table
    print(f"\nSI table → {args.reports_dir}")
    write_si_table(df, status, smiles_assigned, bench_extra_cols, args.reports_dir, args.dataset_name)


if __name__ == "__main__":
    main()

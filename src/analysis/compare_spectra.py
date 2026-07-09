#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import numpy as np

try:
    import ms_entropy as _ms_entropy
    _HAS_MS_ENTROPY = True
except ImportError:
    _HAS_MS_ENTROPY = False

# ----------------- Spectra reading & similarity functions -----------------

def read_spectrum(filepath, bin_width=1, basepeak=None):
    """Read mz,intensity CSV. Return dict {bin: intensity}.

    bin_width : int, m/z bin size (1 = standard 1 Da rounding, 10 = 10 Da bins)
    basepeak  : if set, rescale so max intensity == basepeak before returning
    """
    if not os.path.isfile(filepath):
        return None
    bins = {}
    with open(filepath, 'r') as f:
        for line in f:
            try:
                x, y = map(float, line.strip().split(','))
                b = int(x // bin_width) * bin_width if bin_width > 1 else int(round(x))
                bins[b] = bins.get(b, 0.0) + y
            except:
                continue
    if not bins:
        return None
    if basepeak is not None:
        max_i = max(bins.values())
        if max_i > 0:
            bins = {k: v * basepeak / max_i for k, v in bins.items()}
    return bins


def cosine_similarity(spec1, spec2):
    if spec1 is None or spec2 is None:
        return np.nan
    all_mz = sorted(set(spec1.keys()) | set(spec2.keys()))
    vec1 = np.array([spec1.get(mz, 0.0) for mz in all_mz])
    vec2 = np.array([spec2.get(mz, 0.0) for mz in all_mz])
    dot = np.dot(vec1, vec2)
    norm = np.linalg.norm(vec1) * np.linalg.norm(vec2)
    return 1000 * dot / norm if norm != 0 else np.nan


def weighted_dot(spec1, spec2, m_exp=3.0, i_exp=0.6):
    """
    Weighted dot-product similarity (Stein & Scott style).

    Dot product is computed over matched peaks only.
    Norms are computed over each spectrum's full peak list independently.

    Parameters:
        spec1, spec2 : dict {mz: intensity}
        m_exp        : float, exponent for m/z (default 3.0)
        i_exp        : float, exponent for intensity (default 0.6)

    Returns:
        float : weighted dot similarity (0-1000 scale)
    """
    if spec1 is None or spec2 is None:
        return np.nan

    matched_mz = set(spec1.keys()) & set(spec2.keys())
    if not matched_mz:
        return np.nan

    dot = sum(
        (spec1[mz] * spec2[mz]) ** i_exp * (mz ** 2) ** m_exp
        for mz in matched_mz
    )
    norm1 = sum((spec1[mz] ** i_exp * mz ** m_exp) ** 2 for mz in spec1)
    norm2 = sum((spec2[mz] ** i_exp * mz ** m_exp) ** 2 for mz in spec2)

    norm = (norm1 * norm2) ** 0.5
    return 1000 * dot / norm if norm != 0 else np.nan


def tanimoto_index(spec1, spec2):
    if spec1 is None or spec2 is None:
        return np.nan
    all_mz = sorted(set(spec1.keys()) | set(spec2.keys()))
    bin1 = np.array([1 if spec1.get(mz, 0.0) > 0 else 0 for mz in all_mz])
    bin2 = np.array([1 if spec2.get(mz, 0.0) > 0 else 0 for mz in all_mz])
    inter = np.sum(bin1 * bin2)
    union = np.sum(bin1) + np.sum(bin2) - inter
    return inter / union if union != 0 else np.nan


def fraction_sim_in_ref(spec_ref, spec_sim):
    if spec_ref is None or spec_sim is None:
        return np.nan
    ref_peaks = set(spec_ref.keys())
    sim_peaks = set(spec_sim.keys())
    if not sim_peaks:
        return np.nan
    return 100 * len(sim_peaks & ref_peaks) / len(sim_peaks)


def fraction_ref_in_sim(spec_ref, spec_sim):
    if spec_ref is None or spec_sim is None:
        return np.nan
    ref_peaks = set(spec_ref.keys())
    sim_peaks = set(spec_sim.keys())
    if not ref_peaks:
        return np.nan
    return 100 * len(ref_peaks & sim_peaks) / len(ref_peaks)


def entropy_similarity(spec_ref, spec_sim):
    """MS entropy similarity (0–1 scale). Requires ms_entropy package."""
    if not _HAS_MS_ENTROPY or spec_ref is None or spec_sim is None:
        return np.nan
    arr_ref = np.array([[mz, i] for mz, i in spec_ref.items()], dtype=np.float32)
    arr_sim = np.array([[mz, i] for mz, i in spec_sim.items()], dtype=np.float32)
    if len(arr_ref) == 0 or len(arr_sim) == 0:
        return np.nan
    return float(_ms_entropy.calculate_entropy_similarity(
        arr_ref, arr_sim, ms2_tolerance_in_da=0.05, clean_spectra=True
    ))


# ----------------- Main script -----------------

def main():
    parser = argparse.ArgumentParser(description="Compare spectra for multiple methods and molecules.")
    parser.add_argument("--base_dir", required=True, help="Base dataset directory")
    parser.add_argument("--include_all_peaks", action="store_true", help="Include all peaks spectra")
    parser.add_argument("--include_top20", action="store_true", help="Include top 20 peaks spectra")
    parser.add_argument("--include_top10", action="store_true", help="Include top 10 percent peaks spectra")
    parser.add_argument("--methods", nargs="+", default=None, help="List of methods to include (optional)")
    parser.add_argument("--bin_width", type=int, default=1, help="m/z bin width in Da (default 1; use 10 for sensitivity check)")
    parser.add_argument("--basepeak", type=float, default=None, help="Rescale spectra to this basepeak intensity before comparing (e.g. 999)")
    args = parser.parse_args()

    base_dir = os.path.abspath(args.base_dir)
    _suffix = f"_bin{args.bin_width}bp{int(args.basepeak)}" if (args.bin_width != 1 or args.basepeak is not None) else ""
    results_base = os.path.join(base_dir, f"results{_suffix}")
    os.makedirs(results_base, exist_ok=True)

    # Peak types to process
    peak_flags = []
    if args.include_all_peaks:
        peak_flags.append("all")
    if args.include_top20:
        peak_flags.append("top20")
    if args.include_top10:
        peak_flags.append("10pct")
    if not peak_flags:
        peak_flags = ["all", "top20", "10pct"]

    # Detect methods
    methods = args.methods or [
        d for d in os.listdir(base_dir)
        if os.path.isdir(os.path.join(base_dir, d))
        and d not in ["EXP"]
        and not d.startswith("results")
    ]

    # Detect molecules from EXP
    exp_dir = os.path.join(base_dir, "EXP")
    molecules = sorted(
        d for d in os.listdir(exp_dir)
        if os.path.isdir(os.path.join(exp_dir, d)) and d.isdigit()
    )

    for mol in molecules:
        mol_results_dir = os.path.join(results_base, mol)
        os.makedirs(mol_results_dir, exist_ok=True)

        ref_dir = os.path.join(exp_dir, mol, "spectra")
        ref_spectra = {
            "all":   read_spectrum(os.path.join(ref_dir, "spectra_all.csv"),   args.bin_width, args.basepeak),
            "top20": read_spectrum(os.path.join(ref_dir, "spectra_top20.csv"), args.bin_width, args.basepeak),
            "10pct": read_spectrum(os.path.join(ref_dir, "spectra_10pct.csv"), args.bin_width, args.basepeak),
        }

        for pt in peak_flags:
            rows = []
            for method in methods:
                sim_file = os.path.join(base_dir, method, mol, "spectra", f"spectra_{pt}.csv")
                print(f"[{mol}/{pt}] REF={os.path.join(ref_dir, f'spectra_{pt}.csv')}  SIM={sim_file}")
                spec_sim = read_spectrum(sim_file, args.bin_width, args.basepeak)

                row = {
                    "Method":            method,
                    "Cosine":            cosine_similarity(ref_spectra[pt], spec_sim),
                    "Weighted_Dot":      weighted_dot(ref_spectra[pt], spec_sim),
                    "Tanimoto":          tanimoto_index(ref_spectra[pt], spec_sim),
                    "%S_sim_in_ref":     fraction_sim_in_ref(ref_spectra[pt], spec_sim),
                    "%R_ref_in_sim":     fraction_ref_in_sim(ref_spectra[pt], spec_sim),
                    "Entropy_Similarity": entropy_similarity(ref_spectra[pt], spec_sim),
                }
                rows.append(row)

            output_file = os.path.join(mol_results_dir, f"spectra_{pt}_comparison{_suffix}.csv")
            df_new = pd.DataFrame(rows)

            # Carry over existing Entropy_Similarity values rather than overwriting with NaN.
            # ms_entropy requires Python 3.8+; run patch_entropy.py to compute fresh values.
            if os.path.isfile(output_file):
                try:
                    df_old = pd.read_csv(output_file)
                    if "Entropy_Similarity" in df_old.columns:
                        old_entropy = df_old.set_index("Method")["Entropy_Similarity"]
                        df_new["Entropy_Similarity"] = df_new["Method"].map(old_entropy)
                        print(f"[{mol}] {pt} — kept existing Entropy_Similarity values "
                              f"(run patch_entropy.py to recompute)")
                except Exception:
                    pass

            df_new.to_csv(output_file, index=False)
            print(f"[{mol}] {pt} results written -> {output_file}")


if __name__ == "__main__":
    main()
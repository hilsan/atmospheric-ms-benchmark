import argparse
import numpy as np
import os
import pandas as pd
from collections import defaultdict

def calculate_cosine_similarity(simulated_intensity, reference_intensity):
    sim = np.array(simulated_intensity)
    ref = np.array(reference_intensity)
    dot = np.dot(sim, ref)
    norm_sim = np.linalg.norm(sim)
    norm_ref = np.linalg.norm(ref)
    if norm_sim == 0 or norm_ref == 0:
        return 0.0
    return dot / (norm_sim * norm_ref) * 1000

def calculate_weighted_cosine_similarity(sim_mz, sim_int, ref_mz, ref_int):
    weighted_sim = np.array(sim_int) * np.array(sim_mz)
    weighted_ref = np.array(ref_int) * np.array(ref_mz)
    dot = np.dot(weighted_sim, weighted_ref)
    norm_sim = np.linalg.norm(weighted_sim)
    norm_ref = np.linalg.norm(weighted_ref)
    if norm_sim == 0 or norm_ref == 0:
        return 0.0
    return dot / (norm_sim * norm_ref) * 1000

def calculate_tanimoto_similarity(simulated_intensity, reference_intensity):
    sim = np.array(simulated_intensity)
    ref = np.array(reference_intensity)
    intersection = np.sum(sim * ref)
    union = np.sum(sim**2) + np.sum(ref**2) - intersection
    if union == 0:
        return 0.0
    return intersection / union * 1000

def bin_and_normalize(mz_values, intensities, threshold=1.0, top_n=None):
    binned = defaultdict(float)
    for mz, intensity in zip(mz_values, intensities):
        binned[round(mz)] += intensity  # Bin to nearest integer

    if not binned:
        return [], []

    max_int = max(binned.values())
    norm_binned = {mz: intensity / max_int * 999 for mz, intensity in binned.items()}
    norm_binned = {mz: intensity for mz, intensity in norm_binned.items() if intensity >= threshold}

    if top_n and top_n > 0:
        norm_binned = dict(sorted(norm_binned.items(), key=lambda x: -x[1])[:top_n])

    if not norm_binned:
        return [], []

    sorted_mz = sorted(norm_binned)
    sorted_int = [norm_binned[mz] for mz in sorted_mz]

    return sorted_mz, sorted_int


def parse_spectrum(filepath, threshold, top_n):
    ext = os.path.splitext(filepath)[1].lower()
    mz, intensity = [], []
    with open(filepath, 'r') as f:
        for line in f:
            try:
                if ext == '.csv':
                    x, y = map(float, line.strip().split(','))
                else:
                    x, y = map(float, line.strip().split())
                mz.append(x)
                intensity.append(y)
            except:
                continue
    return bin_and_normalize(mz, intensity, threshold, top_n)

def similarity_scores(spec1, spec2):
    if spec1 is None or spec2 is None:
        return {
            "Cosine": np.nan,
            "Weighted_Cosine": np.nan,
            "Tanimoto": np.nan
        }

    mz1, int1 = spec1
    mz2, int2 = spec2
    common_mz = sorted(set(mz1) | set(mz2))
    int1_aligned = [int1[mz1.index(mz)] if mz in mz1 else 0 for mz in common_mz]
    int2_aligned = [int2[mz2.index(mz)] if mz in mz2 else 0 for mz in common_mz]
    return {
        "Cosine": calculate_cosine_similarity(int1_aligned, int2_aligned),
        "Weighted_Cosine": calculate_weighted_cosine_similarity(common_mz, int1_aligned, common_mz, int2_aligned),
        "Tanimoto": calculate_tanimoto_similarity(int1_aligned, int2_aligned)
    }

def count_peaks(spectrum):
    return len(spectrum[0]) if spectrum else 0

def simulated_peaks_in_reference(reference, simulated):
    """
    Compute percentage of simulated peaks present in reference peaks.
    Both inputs are (mz, intensities) tuples after binning/thresholding.
    """
    if reference is None or simulated is None:
        return np.nan

    ref_mz, _ = reference
    sim_mz, _ = simulated

    if not sim_mz:
        return np.nan

    overlap = len(set(sim_mz) & set(ref_mz))
    return (overlap / len(sim_mz)) * 100.0


def peak_overlap_percentage(reference, simulated):
    """
    Compute percentage of reference peaks present in simulated peaks.
    Both inputs are (mz, intensities) tuples after binning/thresholding.
    """
    if reference is None or simulated is None:
        return np.nan

    ref_mz, _ = reference
    sim_mz, _ = simulated

    if not ref_mz:
        return np.nan

    ref_set = set(ref_mz)
    sim_set = set(sim_mz)

    overlap = len(ref_set & sim_set)  # intersection
    return (overlap / len(ref_set)) * 100.0


def main():
    parser = argparse.ArgumentParser(description="Compare mass spectra.")
    parser.add_argument("--simulated", required=True, nargs='+', help="Paths to simulated spectrum files")
    parser.add_argument("--reference", required=True, help="Path to the reference spectrum file")
    parser.add_argument("--output", required=True, help="Prefix for output files (CSV)")
    parser.add_argument("--threshold", type=float, default=1.0, help="Minimum normalized intensity threshold (default=1.0)")
    parser.add_argument("--top_n_peaks", type=int, default=0, help="Keep only the top N most intense peaks (default=20). Set to 0 to keep all peaks.")
    args = parser.parse_args()

    # ---------------- Load spectra ----------------
    spectra = {}

    # Reference spectrum
    if os.path.isfile(args.reference):
        reference_spectrum = parse_spectrum(args.reference, args.threshold, args.top_n_peaks)
    else:
        print(f"Warning: Reference file '{args.reference}' not found. Using empty spectrum.")
        reference_spectrum = None
    spectra["reference"] = reference_spectrum

    # Simulated spectra
    for f in args.simulated:
        base_name = os.path.basename(f)
        if os.path.isfile(f):
            spectra[base_name] = parse_spectrum(f, args.threshold, args.top_n_peaks)
        else:
            print(f"Warning: Simulated file '{f}' not found. Will use NaN for similarities.")
            spectra[base_name] = None

    names = list(spectra.keys())

    # ---------------- Similarity matrices ----------------
    cosine_matrix = pd.DataFrame(index=names, columns=names, dtype=float)
    weighted_cosine_matrix = pd.DataFrame(index=names, columns=names, dtype=float)
    tanimoto_matrix = pd.DataFrame(index=names, columns=names, dtype=float)

    for i, name1 in enumerate(names):
        for j, name2 in enumerate(names):
            scores = similarity_scores(spectra[name1], spectra[name2])
            cosine_matrix.loc[name1, name2] = scores["Cosine"]
            weighted_cosine_matrix.loc[name1, name2] = scores["Weighted_Cosine"]
            tanimoto_matrix.loc[name1, name2] = scores["Tanimoto"]

    # Save similarity matrices
    cosine_matrix.to_csv(f"{args.output}_cosine_thresh{int(args.threshold)}.csv")
    weighted_cosine_matrix.to_csv(f"{args.output}_weighted_cosine_thresh{int(args.threshold)}.csv")
    tanimoto_matrix.to_csv(f"{args.output}_tanimoto_thresh{int(args.threshold)}.csv")

    # ---------------- Peak counts ----------------
    peak_counts = {name: count_peaks(spectra[name]) for name in spectra}
    peak_df = pd.DataFrame.from_dict(peak_counts, orient='index', columns=['Peak_Count'])
    peak_df.to_csv(f"{args.output}_peak_counts_thresh{int(args.threshold)}.csv")

    # ---------------- Peak differences vs reference ----------------
    if "reference" in peak_counts:
        ref_peaks = peak_counts["reference"]
        diffs = {
            name: count - ref_peaks
            for name, count in peak_counts.items()
            if name != "reference"
        }
        diff_df = pd.DataFrame.from_dict(diffs, orient='index', columns=['Peak_Diff_vs_Reference'])
        diff_df['Reference_Peaks'] = ref_peaks
        diff_df['Simulated_Peaks'] = peak_df.loc[diff_df.index, 'Peak_Count']
        diff_df['Threshold'] = args.threshold
        diff_df.to_csv(f"{args.output}_peak_differences_thresh{int(args.threshold)}.csv")

        # ---------------- Peak overlap percentages ----------------
        # Fraction of reference peaks present in simulated spectrum
        pct_ref_overlap = {
            name: peak_overlap_percentage(reference_spectrum, spectra[name])
            for name in spectra if name != "reference"
        }

        # Fraction of simulated peaks present in reference spectrum (reverse)
        pct_sim_in_ref = {
            name: simulated_peaks_in_reference(reference_spectrum, spectra[name])
            for name in spectra if name != "reference"
        }

        overlap_df = pd.DataFrame({
            'Pct_Ref_Peaks_Recovered': pct_ref_overlap,
            'Pct_Sim_Peaks_in_Exp': pct_sim_in_ref
        })
        overlap_df['Reference_Peaks'] = ref_peaks
        overlap_df['Threshold'] = args.threshold

        overlap_df[['Pct_Ref_Peaks_Recovered']].to_csv(f"{args.output}_pct_ref_peaks.csv")
        overlap_df[['Pct_Sim_Peaks_in_Exp']].to_csv(f"{args.output}_pct_sim_in_ref.csv")


if __name__ == '__main__':
    main()


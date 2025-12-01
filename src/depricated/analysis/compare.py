import argparse
import numpy as np
import os
import pandas as pd

def parse_msp(filepath):
    mz_values, intensities = [], []
    with open(filepath, 'r') as file:
        in_ions_section = False
        for line in file:
            line = line.strip()
            if line.startswith('Begin Ions'):
                in_ions_section = True
                continue
            if line.startswith('End Ions'):
                in_ions_section = False
                continue
            if in_ions_section and line:
                try:
                    mz, intensity = map(float, line.split())
                    mz_values.append(mz)
                    intensities.append(intensity)
                except ValueError:
                    continue
    return mz_values, intensities

def parse_csv(filepath):
    mz_values, intensities = [], []
    with open(filepath, 'r') as file:
        for line in file:
            try:
                mz, intensity = map(float, line.strip().split(','))
                mz_values.append(mz)
                intensities.append(intensity)
            except ValueError:
                pass
    return mz_values, intensities

def parse_sdf(filepath):
    mz_values, intensities = [], []
    with open(filepath, 'r') as file:
        for line in file:
            if "V2000" in line or line.startswith("M  END"):
                continue
            try:
                mz, intensity = map(float, line.strip().split())
                mz_values.append(mz)
                intensities.append(intensity)
            except ValueError:
                pass
    return mz_values, intensities

def load_spectrum(filepath):
    _, ext = os.path.splitext(filepath)
    ext = ext.lower()
    if ext == ".msp":
        return parse_msp(filepath)
    elif ext == ".csv":
        return parse_csv(filepath)
    elif ext == ".sdf":
        return parse_sdf(filepath)
    else:
        raise ValueError(f"Unsupported file format: {ext}")

def bin_and_normalize(mz_values, intensities, max_intensity=999, intensity_threshold=1.0):
    # Bin to nearest integer m/z
    binned_dict = {}
    for mz, intensity in zip(mz_values, intensities):
        bin_mz = int(round(mz))
        if bin_mz in binned_dict:
            binned_dict[bin_mz] += intensity
        else:
            binned_dict[bin_mz] = intensity

    # Convert dict to sorted arrays
    binned_mz = np.array(sorted(binned_dict.keys()))
    binned_intensity = np.array([binned_dict[mz] for mz in binned_mz])

    # Normalize intensities
    max_val = binned_intensity.max() if binned_intensity.size > 0 else 0
    if max_val > 0:
        binned_intensity = (binned_intensity / max_val) * max_intensity
    else:
        binned_intensity = binned_intensity

    # Remove peaks below threshold
    mask = binned_intensity >= intensity_threshold
    return binned_mz[mask], binned_intensity[mask]

def build_intensity_vector(mz_bins, intensities, full_range):
    """
    Build intensity vector aligned on full_range m/z bins.
    """
    intensity_vector = np.zeros(len(full_range))
    mz_to_idx = {mz: i for i, mz in enumerate(full_range)}
    for mz, intensity in zip(mz_bins, intensities):
        if mz in mz_to_idx:
            intensity_vector[mz_to_idx[mz]] = intensity
    return intensity_vector

def cosine_similarity(vec1, vec2):
    dot = np.dot(vec1, vec2)
    norm1 = np.linalg.norm(vec1)
    norm2 = np.linalg.norm(vec2)
    if norm1 == 0 or norm2 == 0:
        return 0.0
    return (dot / (norm1 * norm2)) * 1000

def tanimoto_similarity(vec1, vec2):
    intersection = np.sum(vec1 * vec2)
    sum1 = np.sum(vec1 ** 2)
    sum2 = np.sum(vec2 ** 2)
    denom = sum1 + sum2 - intersection
    if denom == 0:
        return 0.0
    return (intersection / denom) * 1000

def compute_similarity_matrices(spectra_vectors):
    n = len(spectra_vectors)
    cosine_mat = np.zeros((n, n))
    tanimoto_mat = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            cosine_mat[i, j] = cosine_similarity(spectra_vectors[i], spectra_vectors[j])
            tanimoto_mat[i, j] = tanimoto_similarity(spectra_vectors[i], spectra_vectors[j])
    return cosine_mat, tanimoto_mat

def main():
    parser = argparse.ArgumentParser(description="Compute cosine and tanimoto similarity matrices for spectra.")
    parser.add_argument("--spectra", nargs="+", required=True, help="Paths to spectrum files (msp, csv, sdf).")
    parser.add_argument("--output_prefix", required=True, help="Prefix for output CSV files.")
    args = parser.parse_args()

    # Load and preprocess spectra
    processed_spectra = []
    all_bins = set()
    for filepath in args.spectra:
        mzs, intensities = load_spectrum(filepath)
        binned_mz, binned_intensity = bin_and_normalize(mzs, intensities)
        processed_spectra.append((binned_mz, binned_intensity))
        all_bins.update(binned_mz.tolist())

    full_range = np.array(sorted(all_bins))

    # Build aligned intensity vectors
    spectra_vectors = [build_intensity_vector(mz_bins, intensities, full_range)
                       for mz_bins, intensities in processed_spectra]

    # Compute similarity matrices
    cosine_mat, tanimoto_mat = compute_similarity_matrices(spectra_vectors)

    # Save to CSV
    spectrum_names = [os.path.basename(f) for f in args.spectra]
    cosine_df = pd.DataFrame(cosine_mat, index=spectrum_names, columns=spectrum_names)
    tanimoto_df = pd.DataFrame(tanimoto_mat, index=spectrum_names, columns=spectrum_names)

    cosine_csv = f"{args.output_prefix}_cosine_matrix.csv"
    tanimoto_csv = f"{args.output_prefix}_tanimoto_matrix.csv"

    cosine_df.to_csv(cosine_csv)
    tanimoto_df.to_csv(tanimoto_csv)

    print(f"Saved cosine similarity matrix to {cosine_csv}")
    print(f"Saved tanimoto similarity matrix to {tanimoto_csv}")

if __name__ == "__main__":
    main()


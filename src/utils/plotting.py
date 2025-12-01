"""
Shared utility functions for mass spectra visualization:
- Loading spectra from CSV, MSP, SDF
- Filtering, binning, and normalizing intensities
"""

import numpy as np
from collections import defaultdict
import os

def load_spectrum_csv(filepath):
    mz, intensity = [], []
    with open(filepath, 'r') as f:
        for line in f:
            try:
                m, i = map(float, line.strip().split(','))
                mz.append(m)
                intensity.append(i)
            except ValueError:
                continue
    return np.array(mz), np.array(intensity)

def load_spectrum_msp(filepath):
    mz, intensity = [], []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith("Name") or line.startswith("Comment") or "Ions" in line:
                continue
            try:
                m, i = map(float, line.strip().split())
                mz.append(m)
                intensity.append(i)
            except ValueError:
                continue
    return np.array(mz), np.array(intensity)

def load_spectrum_sdf(filepath):
    mz, intensity = [], []
    reading = False
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if "PREDICTED SPECTRUM" in line:
                reading = True
                continue
            if reading and line == "":
                break
            if reading:
                try:
                    m, i = map(float, line.split())
                    mz.append(m)
                    intensity.append(i)
                except ValueError:
                    continue
    return np.array(mz), np.array(intensity)

def filter_peaks(mz_values, intensities, threshold=1.0):
    mask = intensities >= threshold
    return mz_values[mask], intensities[mask]

def bin_and_normalize(mz_values, intensities, threshold=1.0, top_n=20):
    # Bin to nearest integer
    binned = defaultdict(float)
    for m, i in zip(mz_values, intensities):
        binned[round(m)] += i

    if not binned:
        return np.array([]), np.array([])

    max_int = max(binned.values())
    norm_binned = {mz_: inten / max_int * 999 for mz_, inten in binned.items() if inten / max_int * 999 >= threshold}

    if top_n > 0:
        norm_binned = dict(sorted(norm_binned.items(), key=lambda x: -x[1])[:top_n])

    sorted_mz = np.array(sorted(norm_binned.keys()))
    sorted_int = np.array([norm_binned[m] for m in sorted_mz])

    return sorted_mz, sorted_int

def get_file_format(filepath):
    _, ext = os.path.splitext(filepath)
    ext = ext.lower()
    if ext == ".csv":
        return "csv"
    elif ext == ".msp":
        return "msp"
    elif ext == ".sdf":
        return "sdf"
    else:
        raise ValueError(f"Unsupported file format: {ext}")


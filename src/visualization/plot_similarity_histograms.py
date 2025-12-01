#!/usr/bin/env python3
"""
Plot similarity score histograms for multiple methods.
Usage:
python plot_similarity_histograms.py --dir data_dir --output histograms.png --qcxms_flag --qcxms2_flag
"""

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def collect_data_for_methods(directory, methods_files):
    """Collect data from CSV files across subdirectories for each method."""
    data_dict = {}
    for method, file_name in methods_files.items():
        method_data = []
        for subdir in sorted(os.listdir(directory)):
            subdir_path = os.path.join(directory, subdir)
            file_path = os.path.join(subdir_path, file_name)
            if os.path.isdir(subdir_path) and os.path.exists(file_path):
                try:
                    df = pd.read_csv(file_path)
                    method_data.append(df)
                except Exception as e:
                    print(f"Could not read {file_path}: {e}")
        if method_data:
            data_dict[method] = pd.concat(method_data, ignore_index=True)
    return data_dict

def plot_similarity_histograms(data_dict, output_dir):
    metrics = ["Cosine", "Weighted_Cosine", "Tanimoto", "Dot_Product"]
    bins = np.linspace(0, 1000, 21)
    colors = {"NEIMS": "#E059B0", "QCxMS": "#8CD3FA", "QCxMS2": "#A9CD6F", "QCxMS2-dft": "#FF7F0E"}

    for metric in metrics:
        fig, ax = plt.subplots(figsize=(10,7.5), dpi=100)
        for method, df in data_dict.items():
            if metric not in df.columns:
                continue
            counts, _ = np.histogram(df[metric], bins=bins)
            normalized_counts = counts / counts.sum()
            ax.bar(bins[:-1], normalized_counts, width=np.diff(bins), label=method,
                   color=colors.get(method, "gray"), alpha=0.6)

        ax.set_title(f"{metric} Similarity", fontsize=18)
        ax.set_xlabel("Similarity Score", fontsize=16)
        ax.set_ylabel("Frequency (Percentage)", fontsize=16)
        ax.set_xlim(0,1000)
        ax.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{metric}_similarity_histogram.png"), dpi=300, transparent=True)
        plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description="Plot similarity histograms.")
    parser.add_argument("--dir", required=True, help="Parent directory with subdirectories")
    parser.add_argument("--output", required=True, help="Output directory for histograms")
    parser.add_argument("--qcxms_flag", action="store_true")
    parser.add_argument("--qcxms2_flag", action="store_true")
    parser.add_argument("--qcxms2_dft_flag", action="store_true")
    parser.add_argument("--neims", default="exp2neims_similarity_scores.csv")
    parser.add_argument("--qcxms", default="exp2qcxms_similarity_scores.csv")
    parser.add_argument("--qcxms2", default="exp2qcxms2_similarity_scores.csv")
    parser.add_argument("--qcxms2_dft", default="exp2qcxms-dft_similarity_scores.csv")
    args = parser.parse_args()

    methods_files = {"NEIMS": args.neims}
    if args.qcxms_flag: methods_files["QCxMS"] = args.qcxms
    if args.qcxms2_flag: methods_files["QCxMS2"] = args.qcxms2
    if args.qcxms2_dft_flag: methods_files["QCxMS2-dft"] = args.qcxms2_dft

    data = collect_data_for_methods(args.dir, methods_files)
    plot_similarity_histograms(data, args.output)

if __name__ == "__main__":
    main()


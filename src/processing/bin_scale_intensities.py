#!/usr/bin/env python3
import argparse
import pandas as pd
import os
import sys
import re


def read_input_file(input_file):
    """Detect file type and return DataFrame with columns mz, intensity."""
    ext = os.path.splitext(input_file)[1].lower()

    try:
        if ext == ".csv":
            df = pd.read_csv(input_file, header=None, names=["mz", "intensity"])

        elif ext == ".sdf":
            with open(input_file) as f:
                lines = f.readlines()
            spectrum_lines = []
            in_spectrum = False
            for line in lines:
                if line.startswith(">  <PREDICTED SPECTRUM>"):
                    in_spectrum = True
                    continue
                if in_spectrum:
                    if line.strip() == "$$$$":
                        break
                    if line.strip():
                        spectrum_lines.append(line.strip())
            mz_list, intensity_list = zip(*[map(float, l.split()) for l in spectrum_lines])
            df = pd.DataFrame({"mz": mz_list, "intensity": intensity_list})

        elif ext == ".msp":
            with open(input_file) as f:
                lines = f.readlines()
            mz_list, intensity_list = [], []
            in_peaks = False
            for line in lines:
                line = line.strip()
                # Handle both "Num Peaks:" and "Begin Ions" as block starters
                if line.lower().startswith("num peaks") or line.lower() == "begin ions":
                    in_peaks = True
                    continue
                # Handle "End Ions" as block ender
                if line.lower() == "end ions":
                    in_peaks = False
                    continue
                if not in_peaks:
                    continue
                # Handle both "mz intensity" and "mz intensity; mz intensity; ..." formats
                for pair in line.split(";"):
                    pair = pair.strip()
                    if pair:
                        parts = pair.split()
                        if len(parts) == 2:
                            try:
                                mz_list.append(float(parts[0]))
                                intensity_list.append(float(parts[1]))
                            except ValueError:
                                continue
            df = pd.DataFrame({"mz": mz_list, "intensity": intensity_list})

        elif ext == ".log":
            mz_list, intensity_list = [], []
            with open(input_file) as f:
                for line in f:
                    line = line.strip()
                    if not line or not re.match(r"^\d", line):
                        continue
                    parts = line.split()
                    if len(parts) == 2:
                        mz_list.append(float(parts[0]))
                        intensity_list.append(float(parts[1]))
            df = pd.DataFrame({"mz": mz_list, "intensity": intensity_list})

        else:
            print(f"Unsupported file type: {input_file}", file=sys.stderr)
            sys.exit(1)

    except Exception as e:
        print(f"Error reading {input_file}: {e}", file=sys.stderr)
        sys.exit(1)

    return df


def bin_and_scale(df):
    """
    Bin to integer m/z (summing colliding peaks), scale base peak to 999,
    keep intensities as floats, and drop zero-intensity peaks.
    """
    df = df.copy()

    # Integer bin m/z
    df["mz"] = df["mz"].round().astype(int)

    # Sum intensities for any colliding bins
    df = df.groupby("mz", as_index=False)["intensity"].sum()

    # Scale to base peak = 999 (keep as float)
    base_peak = df["intensity"].max()
    if base_peak <= 0:
        print("Warning: all intensities are zero or negative.", file=sys.stderr)
        sys.exit(1)
    df["intensity"] = df["intensity"] * 999.0 / base_peak

    # Drop zero or negligible peaks (rounding artifacts from source data)
    df = df[df["intensity"] > 0].reset_index(drop=True)

    return df.sort_values("mz").reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(
        description="Bin spectrum to integer m/z, scale intensities to 999, and export filtered spectra."
    )
    parser.add_argument("-i", "--input", required=True, help="Input spectrum file (.csv, .sdf, .msp, .log)")
    parser.add_argument("-o", "--output", required=True, help="Output prefix (e.g. path/to/spectra)")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    df = read_input_file(args.input)
    df = bin_and_scale(df)

    if df.empty:
        print("Error: no peaks remaining after filtering.", file=sys.stderr)
        sys.exit(1)

    # All peaks
    all_df = df.copy()

    # Top 10% of base peak intensity
    p10_df = df[df["intensity"] >= 0.1 * df["intensity"].max()].copy()

    # Top 20 most intense peaks, sorted by m/z
    top20_df = (
        df.sort_values("intensity", ascending=False)
        .head(20)
        .sort_values("mz")
        .copy()
    )

    fmt = "%.4f"
    all_df.to_csv(f"{args.output}_all.csv",     index=False, header=False, float_format=fmt)
    p10_df.to_csv(f"{args.output}_10pct.csv",   index=False, header=False, float_format=fmt)
    top20_df.to_csv(f"{args.output}_top20.csv", index=False, header=False, float_format=fmt)

    print(f"Written: {args.output}_all.csv     ({len(all_df)} peaks)")
    print(f"Written: {args.output}_10pct.csv   ({len(p10_df)} peaks)")
    print(f"Written: {args.output}_top20.csv   ({len(top20_df)} peaks)")


if __name__ == "__main__":
    main()

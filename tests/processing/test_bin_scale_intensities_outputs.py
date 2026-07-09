"""Tests for derive_output_spectra() and the full bin_scale_intensities.py
script, against a 25-peak worked example supplied and confirmed by the
domain owner.

The 25 raw peaks are a representative QCxMS-style result.csv (no binning
collisions, so this isolates the all/10pct/top20 selection logic from
the binning/collision-summing behaviour already covered in
test_bin_scale_intensities.py). The domain owner independently derived
which peaks belong in the 10%-of-base-peak filter and the top-20 filter;
both were verified here to match the actual code's scaled-intensity
output exactly before being locked in as ground truth (their first
attempt had a m/z rounding mistake, corrected and reconfirmed — the
intensity values and 10pct/top20 selection pattern were already correct
throughout).
"""

import subprocess
import sys

import pandas as pd
import pytest

from src.processing.bin_scale_intensities import bin_and_scale, derive_output_spectra

RAW_PEAKS_25 = [
    (26.8, 12.4), (29.3, 44.7), (30.9, 7.6), (39.2, 64.8), (40.6, 209.5),
    (43.4, 499.6), (44.7, 14.9), (50.3, 89.5), (52.6, 22.3), (55.4, 130.2),
    (57.2, 339.8), (59.8, 5.4), (63.4, 47.6), (64.6, 174.6), (67.3, 3.2),
    (69.1, 399.7), (71.4, 59.6), (72.7, 24.8), (75.3, 109.6), (76.6, 2.3),
    (78.6, 379.9), (81.4, 50.4), (83.1, 48.9), (84.8, 7.1), (86.5, 1.3),
]

# (mz, scaled_intensity, in_10pct, in_top20) — confirmed by domain owner
EXPECTED_ROWS = [
    (27, 24.79503603, False, True),
    (29, 89.38210568, False, True),
    (31, 15.19695757, False, True),
    (39, 129.5740592, True, True),
    (41, 418.9161329, True, True),
    (43, 999.0, True, True),
    (45, 29.79403523, False, True),
    (50, 178.9641713, True, True),
    (53, 44.59107286, False, True),
    (55, 260.3478783, True, True),
    (57, 679.4639712, True, True),
    (60, 10.79783827, False, False),
    (63, 95.18094476, False, True),
    (65, 349.1301041, True, True),
    (67, 6.398718975, False, False),
    (69, 799.239992, True, True),
    (71, 119.1761409, True, True),
    (73, 49.59007206, False, True),
    (75, 219.1561249, True, True),
    (77, 4.599079263, False, False),
    (79, 759.6479183, True, True),
    (81, 100.7798239, True, True),
    (83, 97.78042434, False, True),
    (85, 14.19715773, False, False),
    (87, 2.599479584, False, False),
]

EXPECTED_ALL_MZ = [mz for mz, _, _, _ in EXPECTED_ROWS]
EXPECTED_10PCT_MZ = [mz for mz, _, in10, _ in EXPECTED_ROWS if in10]
EXPECTED_TOP20_MZ = [mz for mz, _, _, in20 in EXPECTED_ROWS if in20]


def _binned_scaled_df():
    df = pd.DataFrame(RAW_PEAKS_25, columns=["mz", "intensity"])
    return bin_and_scale(df)


class TestDeriveOutputSpectra:
    def test_all_contains_every_peak(self):
        outputs = derive_output_spectra(_binned_scaled_df())

        assert outputs["all"]["mz"].tolist() == EXPECTED_ALL_MZ

    def test_10pct_matches_confirmed_selection(self):
        outputs = derive_output_spectra(_binned_scaled_df())

        assert outputs["10pct"]["mz"].tolist() == EXPECTED_10PCT_MZ
        assert len(outputs["10pct"]) == 12

    def test_top20_matches_confirmed_selection_and_is_sorted_by_mz(self):
        outputs = derive_output_spectra(_binned_scaled_df())

        assert outputs["top20"]["mz"].tolist() == EXPECTED_TOP20_MZ
        assert len(outputs["top20"]) == 20
        assert outputs["top20"]["mz"].is_monotonic_increasing

    def test_10pct_boundary_case(self):
        # mz=81 (raw 81.4, scaled ~100.78) sits just above the 10% threshold
        # (99.9); mz=83 (raw 83.1, scaled ~97.78) sits just below it. Both
        # were part of the domain owner's deliberately close boundary test.
        outputs = derive_output_spectra(_binned_scaled_df())

        assert 81 in outputs["10pct"]["mz"].tolist()
        assert 83 not in outputs["10pct"]["mz"].tolist()


class TestFullScriptEndToEnd:
    """Runs bin_scale_intensities.py exactly as process_spectra_batch.py
    invokes it in production (subprocess with -i/-o), then reads back the
    three real output CSV files it writes."""

    def test_writes_correct_all_10pct_top20_csv_files(self, tmp_path):
        result_csv = tmp_path / "result.csv"
        result_csv.write_text(
            "\n".join(f"{mz},{intensity}" for mz, intensity in RAW_PEAKS_25) + "\n"
        )
        output_prefix = tmp_path / "spectra"

        script_path = (
            __file__.rsplit("/tests/", 1)[0] + "/src/processing/bin_scale_intensities.py"
        )
        result = subprocess.run(
            [sys.executable, script_path, "-i", str(result_csv), "-o", str(output_prefix)],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr

        all_df = pd.read_csv(f"{output_prefix}_all.csv", header=None, names=["mz", "intensity"])
        p10_df = pd.read_csv(f"{output_prefix}_10pct.csv", header=None, names=["mz", "intensity"])
        top20_df = pd.read_csv(f"{output_prefix}_top20.csv", header=None, names=["mz", "intensity"])

        assert all_df["mz"].tolist() == EXPECTED_ALL_MZ
        assert p10_df["mz"].tolist() == EXPECTED_10PCT_MZ
        assert top20_df["mz"].tolist() == EXPECTED_TOP20_MZ

        # Written with %.4f formatting, so compare with matching tolerance.
        for mz, expected_intensity, _, _ in EXPECTED_ROWS:
            written = all_df.loc[all_df["mz"] == mz, "intensity"].iloc[0]
            assert written == pytest.approx(expected_intensity, abs=1e-4)

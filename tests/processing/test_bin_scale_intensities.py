"""Tests for bin_and_scale() against a hand-verified worked example.

Ground truth (raw peaks, binning rule, and expected output) supplied
and confirmed by the domain owner from a representative QCxMS-style
result.csv: (14, 0.005), (16.1, 0.1), (16.5, 0.2), (17.2, 3), (17.8, 21).

The 16.5 peak is a deliberate exact-half m/z tie: confirmed round-half-up
(16.5 -> 17), not pandas' default banker's rounding (which would give
16.5 -> 16). Expected scaled intensities are computed from the raw
grouped sums via the documented "scale max peak to 999" rule, and were
cross-checked to match the domain owner's independently hand-calculated
values to 9 decimal places:
    mz=14 -> 0.237857143
    mz=16 -> 4.757142857
    mz=17 -> 152.2285714
    mz=18 -> 999

NOTE: as of writing, this test fails against bin_and_scale()'s current
implementation, which rounds 16.5 -> 16 (banker's rounding via pandas
.round()). It will pass once bin_and_scale() is switched to round-half-up
(np.floor(mz + 0.5)) per the agreed fix.
"""

import pandas as pd
import pytest

from src.processing.bin_scale_intensities import bin_and_scale, read_input_file

RAW_PEAKS = [
    (14.0, 0.005),
    (16.1, 0.1),
    (16.5, 0.2),   # exact-half tie: must round up to 17, not down to 16
    (17.2, 3.0),
    (17.8, 21.0),
]

# Grouped sums after round-half-up binning: {14: 0.005, 16: 0.1, 17: 3.2, 18: 21}
# Scale factor = 999 / max(grouped sums) = 999 / 21
_SCALE = 999.0 / 21.0
EXPECTED_MZ = [14, 16, 17, 18]
EXPECTED_INTENSITY = [
    pytest.approx(0.005 * _SCALE),
    pytest.approx(0.1 * _SCALE),
    pytest.approx(3.2 * _SCALE),
    pytest.approx(21.0 * _SCALE),
]


class TestBinAndScale:
    def test_matches_hand_verified_example(self):
        df = pd.DataFrame(RAW_PEAKS, columns=["mz", "intensity"])

        result = bin_and_scale(df)

        assert result["mz"].tolist() == EXPECTED_MZ
        assert result["intensity"].tolist() == EXPECTED_INTENSITY

    def test_base_peak_scales_to_exactly_999(self):
        df = pd.DataFrame(RAW_PEAKS, columns=["mz", "intensity"])

        result = bin_and_scale(df)

        assert result["intensity"].max() == pytest.approx(999.0)


class TestBinAndScaleFromQcxmsResultCsv:
    """End-to-end: a real headerless result.csv, as produced under
    GS-opt/MS-run/result.csv by QCxMS/QCxMS2, read via read_input_file()
    and fed through bin_and_scale() exactly as process_spectra_batch.py
    does for the QCxMS family of datasets."""

    def test_qcxms_result_csv_end_to_end(self, tmp_path):
        result_csv = tmp_path / "result.csv"
        result_csv.write_text(
            "\n".join(f"{mz},{intensity}" for mz, intensity in RAW_PEAKS) + "\n"
        )

        df = read_input_file(str(result_csv))
        result = bin_and_scale(df)

        assert result["mz"].tolist() == EXPECTED_MZ
        assert result["intensity"].tolist() == EXPECTED_INTENSITY

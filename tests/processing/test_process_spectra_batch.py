"""Tests for get_input_file() / spectra_complete() routing logic.

Both functions are pure path-dispatch/existence checks with no
subprocess or SLURM dependency, so expected values are derived
directly from the function's own documented routing table in
process_spectra_batch.py rather than external domain ground truth.
"""

import os

import pytest

from src.processing.process_spectra_batch import (
    EXPECTED_SPECTRA,
    get_input_file,
    spectra_complete,
)


class TestGetInputFile:
    def test_neims_routes_to_annotated_sdf(self, tmp_path):
        subdir_path = str(tmp_path)
        assert get_input_file("NEIMS", "0000", subdir_path) == os.path.join(
            subdir_path, "annotated.sdf"
        )

    def test_neims_is_case_insensitive(self, tmp_path):
        subdir_path = str(tmp_path)
        assert get_input_file("neims", "0000", subdir_path) == os.path.join(
            subdir_path, "annotated.sdf"
        )

    def test_exp_routes_to_exp_msp(self, tmp_path):
        subdir_path = str(tmp_path)
        assert get_input_file("EXP", "0000", subdir_path) == os.path.join(
            subdir_path, "exp.msp"
        )

    @pytest.mark.parametrize("dataset", ["CFM-ID", "cfmid", "CFMID"])
    def test_cfmid_routes_to_subdir_named_log(self, tmp_path, dataset):
        subdir_path = str(tmp_path)
        assert get_input_file(dataset, "0042", subdir_path) == os.path.join(
            subdir_path, "0042.log"
        )

    def test_qcxms2_returns_peaks_csv_when_auxiliary_archive_exists(self, tmp_path):
        subdir_path = str(tmp_path)
        gfn2_dir = tmp_path / "qcxms2_gfn2"
        gfn2_dir.mkdir()
        (gfn2_dir / "qcxms2_auxiliary.tar.gz").write_text("placeholder")

        result = get_input_file("qcxms2", "0000", subdir_path)

        assert result == str(gfn2_dir / "peaks.csv")

    def test_qcxms2_returns_none_when_auxiliary_archive_missing(self, tmp_path):
        # A missing auxiliary archive is the project's proxy for an
        # incomplete/still-running QCxMS2 job, so this must return
        # None rather than a path that doesn't exist yet.
        subdir_path = str(tmp_path)

        result = get_input_file("qcxms2", "0000", subdir_path)

        assert result is None

    @pytest.mark.parametrize(
        "dataset", ["QCxMS_10_ps", "QCxMS_25_ps", "QCxMS_10_ps_iee03", "anything_else"]
    )
    def test_qcxms_variants_route_to_result_csv(self, tmp_path, dataset):
        subdir_path = str(tmp_path)
        expected = os.path.join(subdir_path, "GS-opt", "MS-run", "result.csv")
        assert get_input_file(dataset, "0000", subdir_path) == expected


class TestSpectraComplete:
    def test_false_when_directory_missing(self, tmp_path):
        assert spectra_complete(str(tmp_path / "does_not_exist")) is False

    def test_false_when_some_expected_files_missing(self, tmp_path):
        (tmp_path / EXPECTED_SPECTRA[0]).write_text("mz,intensity\n")
        assert spectra_complete(str(tmp_path)) is False

    def test_true_when_all_expected_files_present(self, tmp_path):
        for fname in EXPECTED_SPECTRA:
            (tmp_path / fname).write_text("mz,intensity\n")
        assert spectra_complete(str(tmp_path)) is True

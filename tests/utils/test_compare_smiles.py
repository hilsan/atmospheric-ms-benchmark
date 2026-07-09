"""Tests for canonical_smiles() and read_smiles_from_file() in
src/utils/compare_smiles.py.

canonical_smiles() here uses RDKit's plain default MolToSmiles(),
producing aromatic lowercase notation - a genuinely different
convention from remove_duplicate_SMILES_entries.py's Kekule-forcing
version (see tests/processing/test_remove_duplicate_smiles_entries.py),
confirmed intentional by the project owner for each function's own
purpose. Not the same function, not meant to agree.

Ground truth: three notations for benzene, confirmed by the project
owner to canonicalize to c1ccccc1.
"""

import pandas as pd
import pytest

from src.utils.compare_smiles import canonical_smiles, read_smiles_from_file

BENZENE_VARIANTS = ["c1ccccc1", "C1=CC=CC=C1", "c1ccc(cc1)"]
BENZENE_AROMATIC = "c1ccccc1"


class TestCanonicalSmilesAromatic:
    @pytest.mark.parametrize("smi", BENZENE_VARIANTS)
    def test_all_benzene_notations_canonicalize_to_aromatic_form(self, smi):
        assert canonical_smiles(smi) == BENZENE_AROMATIC

    def test_invalid_smiles_returns_none(self):
        assert canonical_smiles("not_a_smiles!!!") is None


class TestReadSmilesFromFile:
    def test_reads_and_canonicalizes_detected_smiles_column(self, tmp_path):
        f = tmp_path / "smiles.csv"
        pd.DataFrame({"SMILES": BENZENE_VARIANTS}).to_csv(f, index=False)

        result = read_smiles_from_file(str(f))

        assert result == [BENZENE_AROMATIC] * 3

    def test_detects_smiles_column_case_insensitively(self, tmp_path):
        f = tmp_path / "data.csv"
        pd.DataFrame({"Compound_Smiles": ["c1ccccc1"]}).to_csv(f, index=False)

        assert read_smiles_from_file(str(f)) == [BENZENE_AROMATIC]

    def test_missing_smiles_column_returns_empty_list(self, tmp_path):
        f = tmp_path / "data.csv"
        pd.DataFrame({"name": ["benzene"]}).to_csv(f, index=False)

        assert read_smiles_from_file(str(f)) == []

    def test_missing_file_returns_empty_list(self, tmp_path):
        assert read_smiles_from_file(str(tmp_path / "does_not_exist.csv")) == []

    def test_invalid_smiles_silently_skipped(self, tmp_path):
        f = tmp_path / "smiles.csv"
        pd.DataFrame({"SMILES": ["c1ccccc1", "invalid!!!"]}).to_csv(f, index=False)

        assert read_smiles_from_file(str(f)) == [BENZENE_AROMATIC]

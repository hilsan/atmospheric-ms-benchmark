"""Tests for canonical_smiles() and find_duplicates() in
remove_duplicate_SMILES_entries.py.

canonical_smiles() here deliberately produces Kekule (explicit
alternating bond) SMILES, not aromatic lowercase notation - confirmed
intentional by the project owner (needed for downstream structural
work, e.g. aprl_ssp/derivatization), NOT the same convention as
compare_smiles.py's canonical_smiles() (see
tests/utils/test_compare_smiles.py), which is also intentional and
tested separately against its own correct behavior.

Ground truth: three different SMILES notations for benzene, all
representing the same molecule, confirmed by the project owner to
canonicalize (aromatic form) to c1ccccc1. Verified directly against
RDKit that this function's Kekule-forcing gives C1=CC=CC=C1 for all
three - and confirmed empirically that "canonicalize then kekulize"
and "kekulize directly" are provably equivalent (Kekulization acts on
the molecular graph, not the input string), so this one-step
implementation is a valid shortcut, not a bug.
"""

import pandas as pd
import pytest

from src.processing.remove_duplicate_SMILES_entries import (
    canonical_smiles,
    find_duplicates,
)

# Three notations for benzene -> all confirmed to canonicalize to the
# same Kekule form.
BENZENE_VARIANTS = ["c1ccccc1", "C1=CC=CC=C1", "c1ccc(cc1)"]
BENZENE_KEKULE = "C1=CC=CC=C1"


class TestCanonicalSmilesKekule:
    @pytest.mark.parametrize("smi", BENZENE_VARIANTS)
    def test_all_benzene_notations_canonicalize_to_same_kekule_form(self, smi):
        assert canonical_smiles(smi) == BENZENE_KEKULE

    def test_invalid_smiles_returns_none(self):
        assert canonical_smiles("not_a_smiles!!!") is None


class TestFindDuplicates:
    def test_detects_duplicates_across_equivalent_smiles_notations(self):
        df = pd.DataFrame({"SMILES": BENZENE_VARIANTS})

        to_drop, duplicate_info, canonical_series = find_duplicates(df)

        # First occurrence (index 0) is kept; the other two equivalent
        # notations are flagged as duplicates of it.
        assert to_drop == [1, 2]
        assert len(duplicate_info) == 2
        assert all(d["original_kept_index"] == 0 for d in duplicate_info)
        assert all(canonical_series[i] == BENZENE_KEKULE for i in range(3))

    def test_no_duplicates_for_distinct_molecules(self):
        df = pd.DataFrame({"SMILES": ["c1ccccc1", "CCO", "CC(=O)O"]})

        to_drop, duplicate_info, _ = find_duplicates(df)

        assert to_drop == []
        assert duplicate_info == []

    def test_invalid_smiles_excluded_not_treated_as_duplicates_of_each_other(self):
        df = pd.DataFrame({"SMILES": ["c1ccccc1", "invalid!!!", "also_invalid!!!"]})

        to_drop, duplicate_info, _ = find_duplicates(df)

        # None-canonical entries are skipped entirely, not grouped together.
        assert to_drop == []

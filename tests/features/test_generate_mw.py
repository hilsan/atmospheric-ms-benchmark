"""Tests for count_atoms() in src/features/generate_mw.py.

Ground truth confirmed by the project owner for a real molecule from
this project's data (a xanthophyll-type epoxide-ring polyene):
    CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/C=C/C=C(\\C)/C=C/C=C(\\C)/C=C/C23C(CCCC2(O3)C)(C)C)/C)/C
    -> 41 non-hydrogen (heavy) atoms, 97 atoms total with hydrogens.
Verified directly against RDKit before locking in as the expected value.

NOTE: generate_mw.py originally had no `if __name__ == "__main__":`
guard - argparse/file-processing ran at module import time, which
would crash any import (args.data_path has no default, so
os.path.join(None, ...) raised TypeError). Fixed by wrapping the
script body in main(), which doesn't change CLI behavior at all.
"""

import pytest

from src.features.generate_mw import count_atoms

COMPLEX_POLYENE = (
    r"CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/C=C/C=C(\C)/C=C/C=C(\C)"
    r"/C=C/C23C(CCCC2(O3)C)(C)C)/C)/C"
)


class TestCountAtoms:
    def test_heavy_atom_count(self):
        assert count_atoms(COMPLEX_POLYENE, include_hydrogens=False) == 41

    def test_total_atom_count_with_hydrogens(self):
        assert count_atoms(COMPLEX_POLYENE, include_hydrogens=True) == 97

    def test_default_matches_include_hydrogens_true(self):
        # count_atoms()'s own default is include_hydrogens=True, matching
        # every real usage found in the codebase (always passes
        # --include_hydrogens explicitly).
        assert count_atoms(COMPLEX_POLYENE) == 97

    def test_invalid_smiles_returns_none(self):
        assert count_atoms("not_a_smiles!!!") is None

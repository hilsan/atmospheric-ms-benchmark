"""Tests for copy_unique_folders() — pure file-mapping logic driven by
an old_index -> new_index CSV, no domain ground truth needed beyond
the function's own documented contract.
"""

import pandas as pd
import pytest

from src.processing.copy_unique_folders import copy_unique_folders


def _make_source_folder(base, index, content="placeholder"):
    folder = base / f"{index:04d}"
    folder.mkdir()
    (folder / "structure.sdf").write_text(content)
    return folder


class TestCopyUniqueFolders:
    def test_copies_and_renumbers_according_to_mapping(self, tmp_path):
        source_dir = tmp_path / "source"
        target_dir = tmp_path / "target"
        source_dir.mkdir()

        _make_source_folder(source_dir, 3, content="molecule-3")
        _make_source_folder(source_dir, 7, content="molecule-7")

        mapping_file = tmp_path / "mapping.csv"
        pd.DataFrame({"old_index": [3, 7], "new_index": [0, 1]}).to_csv(
            mapping_file, index=False
        )

        copy_unique_folders(mapping_file, source_dir, target_dir)

        assert (target_dir / "0000" / "structure.sdf").read_text() == "molecule-3"
        assert (target_dir / "0001" / "structure.sdf").read_text() == "molecule-7"

    def test_creates_target_dir_if_missing(self, tmp_path):
        source_dir = tmp_path / "source"
        target_dir = tmp_path / "does" / "not" / "exist" / "yet"
        source_dir.mkdir()
        _make_source_folder(source_dir, 0)

        mapping_file = tmp_path / "mapping.csv"
        pd.DataFrame({"old_index": [0], "new_index": [0]}).to_csv(
            mapping_file, index=False
        )

        copy_unique_folders(mapping_file, source_dir, target_dir)

        assert target_dir.exists()
        assert (target_dir / "0000").exists()

    def test_skips_missing_source_folder_without_raising(self, tmp_path, capsys):
        source_dir = tmp_path / "source"
        target_dir = tmp_path / "target"
        source_dir.mkdir()
        # Deliberately no folder "0005" created under source_dir.

        mapping_file = tmp_path / "mapping.csv"
        pd.DataFrame({"old_index": [5], "new_index": [0]}).to_csv(
            mapping_file, index=False
        )

        copy_unique_folders(mapping_file, source_dir, target_dir)

        assert not (target_dir / "0000").exists()
        assert "Warning" in capsys.readouterr().out

    def test_raises_on_malformed_mapping_csv(self, tmp_path):
        source_dir = tmp_path / "source"
        target_dir = tmp_path / "target"
        source_dir.mkdir()

        mapping_file = tmp_path / "mapping.csv"
        pd.DataFrame({"wrong_column": [1, 2]}).to_csv(mapping_file, index=False)

        with pytest.raises(ValueError, match="old_index.*new_index"):
            copy_unique_folders(mapping_file, source_dir, target_dir)

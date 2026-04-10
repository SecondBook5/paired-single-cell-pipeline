from __future__ import annotations

import pandas as pd
import pytest

from paired_sc import ManifestTable, WorkflowConfig


def test_manifest_loads(smoke_inputs) -> None:
    manifest = ManifestTable.from_csv(smoke_inputs["manifest"])
    assert len(manifest.dataframe) == 4
    assert set(manifest.dataframe["condition"]) == {"Normal", "LE"}


def test_manifest_requires_columns(tmp_path) -> None:
    path = tmp_path / "manifest.csv"
    pd.DataFrame({"sample_id": ["A"]}).to_csv(path, index=False)
    with pytest.raises(ValueError):
        ManifestTable.from_csv(path)


def test_manifest_validates_expected_conditions(smoke_inputs) -> None:
    manifest = ManifestTable.from_csv(smoke_inputs["manifest"])
    config = WorkflowConfig.from_yaml(smoke_inputs["project"])
    manifest.validate_conditions(config)


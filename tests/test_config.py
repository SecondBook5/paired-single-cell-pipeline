from __future__ import annotations

import pytest

from paired_sc import WorkflowConfig


def test_workflow_config_loads(smoke_inputs) -> None:
    config = WorkflowConfig.from_yaml(smoke_inputs["project"])
    assert config.project_name == "smoke_project"
    assert config.condition_order == ["Normal", "LE"]
    assert "trajectory" in config.domains.enabled


def test_workflow_config_rejects_identical_conditions(tmp_path) -> None:
    bad = tmp_path / "bad.yaml"
    bad.write_text(
        "project_name: bad\ncase_condition: LE\ncontrol_condition: LE\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError):
        WorkflowConfig.from_yaml(bad)


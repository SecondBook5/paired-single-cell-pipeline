from __future__ import annotations

from typer.testing import CliRunner

from paired_sc.cli.app import app


runner = CliRunner()


def test_cli_validate(smoke_inputs) -> None:
    result = runner.invoke(
        app,
        [
            "validate",
            "--project",
            str(smoke_inputs["project"]),
            "--manifest",
            str(smoke_inputs["manifest"]),
        ],
    )
    assert result.exit_code == 0
    assert "smoke_project" in result.stdout


def test_cli_domain_status() -> None:
    result = runner.invoke(app, ["domain", "status"])
    assert result.exit_code == 0
    assert "trajectory" in result.stdout


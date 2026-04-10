"""Typer CLI for paired_sc."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import typer

from ..config.models import WorkflowConfig
from ..domains import get_domain_status
from ..io.manifest import ManifestTable
from ..paths import RunPaths
from ..report.core import build_report
from ..workflow import run_advanced_domains, run_core_pipeline

app = typer.Typer(add_completion=False, help="paired_sc CLI")
run_app = typer.Typer(help="Run core workflows and callable domains.")
domain_app = typer.Typer(help="Invoke or inspect callable analysis domains.")
export_app = typer.Typer(help="Export starter config material.")
app.add_typer(run_app, name="run")
app.add_typer(domain_app, name="domain")
app.add_typer(export_app, name="export")


def _load_inputs(project: str | Path, manifest: str | Path) -> tuple[WorkflowConfig, ManifestTable]:
    config = WorkflowConfig.from_yaml(project)
    manifest_table = ManifestTable.from_csv(manifest)
    manifest_table.validate_conditions(config)
    return config, manifest_table


@app.command()
def validate(
    project: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    manifest: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
) -> None:
    """Validate project and manifest inputs."""
    config, manifest_table = _load_inputs(project, manifest)
    typer.echo(
        json.dumps(
            {
                "project_name": config.project_name,
                "manifest": str(manifest_table.path),
                "n_samples": len(manifest_table.dataframe),
                "conditions": config.condition_order,
                "domains_enabled": config.domains.enabled,
            },
            indent=2,
        )
    )


@run_app.command("core")
def run_core(
    project: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    manifest: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    workdir: Path = typer.Option(...),
) -> None:
    """Run the standardized core workflow."""
    outputs = run_core_pipeline(project, manifest, workdir)
    typer.echo(json.dumps(outputs, indent=2))


@run_app.command("advanced")
def run_advanced(
    project: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    manifest: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    workdir: Path = typer.Option(...),
    domains: list[str] = typer.Option(None, help="Specific domain names to run."),
) -> None:
    """Run one or more advanced callable domains."""
    results = run_advanced_domains(project, manifest, workdir, domains=domains or None)
    typer.echo(json.dumps(results, indent=2))


@run_app.command("all")
def run_all(
    project: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    manifest: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    workdir: Path = typer.Option(...),
) -> None:
    """Run the core workflow followed by configured domains."""
    core_outputs = run_core_pipeline(project, manifest, workdir)
    domain_outputs = run_advanced_domains(project, manifest, workdir)
    typer.echo(json.dumps({"core": core_outputs, "domains": domain_outputs}, indent=2))


@app.command()
def report(
    project: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    manifest: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    workdir: Path = typer.Option(...),
) -> None:
    """Rebuild the run report from an existing work directory."""
    config, manifest_table = _load_inputs(project, manifest)
    paths = RunPaths.from_config(workdir, config.outputs)
    core_outputs = {}
    for candidate in sorted(paths.results.glob("*")):
        if candidate.is_file():
            core_outputs[candidate.stem] = str(candidate)
    report_outputs = build_report(config, manifest_table, paths, core_outputs)
    typer.echo(json.dumps(report_outputs, indent=2))


@export_app.command("example-config")
def export_example_config(output_dir: Path = typer.Option(...)) -> None:
    """Copy the flagship example profile into a target directory."""
    package_root = Path(__file__).resolve().parents[3]
    example_root = package_root / "examples" / "skin_paired_human"
    output_dir.mkdir(parents=True, exist_ok=True)
    for item in example_root.iterdir():
        target = output_dir / item.name
        if item.is_dir():
            shutil.copytree(item, target, dirs_exist_ok=True)
        else:
            shutil.copy2(item, target)
    typer.echo(f"Exported example config to {output_dir}")


@domain_app.command("status")
def domain_status() -> None:
    """Show callable domain availability."""
    typer.echo(json.dumps(get_domain_status(), indent=2))


def _run_single_domain(name: str, project: Path, manifest: Path, workdir: Path) -> None:
    results = run_advanced_domains(project, manifest, workdir, domains=[name])
    typer.echo(json.dumps(results, indent=2))


@domain_app.command("liana")
def domain_liana(
    project: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    manifest: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    workdir: Path = typer.Option(...),
) -> None:
    _run_single_domain("liana", project, manifest, workdir)


@domain_app.command("magic")
def domain_magic(
    project: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    manifest: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    workdir: Path = typer.Option(...),
) -> None:
    _run_single_domain("magic", project, manifest, workdir)


@domain_app.command("trajectory")
def domain_trajectory(
    project: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    manifest: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    workdir: Path = typer.Option(...),
) -> None:
    _run_single_domain("trajectory", project, manifest, workdir)


@domain_app.command("latent")
def domain_latent(
    project: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    manifest: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    workdir: Path = typer.Option(...),
) -> None:
    _run_single_domain("latent", project, manifest, workdir)


@domain_app.command("regulatory")
def domain_regulatory(
    project: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    manifest: Path = typer.Option(..., exists=True, file_okay=True, dir_okay=False),
    workdir: Path = typer.Option(...),
) -> None:
    _run_single_domain("regulatory", project, manifest, workdir)


def main() -> None:
    app()


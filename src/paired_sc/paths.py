"""Run path management for paired_sc."""

from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel, ConfigDict


class RunPaths(BaseModel):
    """Resolved run-time directories for a workflow execution."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    root: Path
    results: Path
    figures: Path
    reports: Path
    logs: Path
    domains: Path

    @classmethod
    def from_config(cls, workdir: Path, outputs) -> "RunPaths":
        root = workdir.resolve()
        results = root / outputs.results_dir
        figures = root / outputs.figures_dir
        reports = root / outputs.reports_dir
        logs = root / outputs.logs_dir
        domains = results / "domains"
        for directory in [root, results, figures, reports, logs, domains]:
            directory.mkdir(parents=True, exist_ok=True)
        return cls(
            root=root,
            results=results,
            figures=figures,
            reports=reports,
            logs=logs,
            domains=domains,
        )



"""Report generation."""

from __future__ import annotations

import json


def build_report(config, manifest, paths, core_outputs: dict, domain_results: list[dict] | None = None) -> dict:
    report_md = paths.reports / "run_report.md"
    report_json = paths.reports / "run_report.json"
    domain_results = domain_results or []

    lines = [
        f"# {config.project_name}",
        "",
        config.project_subtitle,
        "",
        "## Run summary",
        f"- Manifest: `{manifest.path}`",
        f"- Samples: `{len(manifest.dataframe)}`",
        f"- Conditions: `{', '.join(config.condition_order)}`",
        f"- Output root: `{paths.root}`",
        "",
        "## Core outputs",
    ]
    for key, value in core_outputs.items():
        lines.append(f"- {key}: `{value}`")

    if domain_results:
        lines.extend(["", "## Domain results"])
        for result in domain_results:
            lines.append(f"- {result['name']}: `{result['status']}` - {result['message']}")

    report_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    payload = {
        "project": config.to_serializable(),
        "manifest_path": str(manifest.path),
        "core_outputs": core_outputs,
        "domain_results": domain_results,
    }
    report_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return {
        "report_markdown": str(report_md),
        "report_json": str(report_json),
    }


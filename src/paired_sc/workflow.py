"""Top-level workflow runners."""

from __future__ import annotations

import json
from pathlib import Path

import scanpy as sc

from .annotate.core import cluster_and_annotate
from .config.models import WorkflowConfig
from .domains.base import DomainContext
from .domains.registry import get_domain_registry
from .io.manifest import ManifestTable, load_manifest_adata
from .paths import RunPaths
from .plotting.core import (
    plot_composition,
    plot_differential_abundance,
    plot_qc_summary,
    plot_umap_annotation,
    plot_umap_condition,
)
from .preprocess.core import preprocess_and_integrate
from .qc.core import compute_qc_metrics, filter_qc
from .report.core import build_report
from .stats.paired import (
    compute_cell_type_composition,
    compute_differential_abundance,
    compute_pseudobulk_profiles,
)


def _save_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def run_core_pipeline(project: str | Path, manifest_csv: str | Path, workdir: str | Path) -> dict:
    """Run the standardized core workflow."""
    config = WorkflowConfig.from_yaml(project)
    manifest = ManifestTable.from_csv(manifest_csv)
    manifest.validate_conditions(config)
    paths = RunPaths.from_config(Path(workdir), config.outputs)

    adata = load_manifest_adata(manifest, config)
    adata = compute_qc_metrics(adata)
    n_cells_raw = int(adata.n_obs)
    n_genes_raw = int(adata.n_vars)
    raw_path = paths.results / "adata_raw_with_qc.h5ad"
    adata.write_h5ad(raw_path)

    qc_figures = plot_qc_summary(adata, config, paths.figures)

    adata = filter_qc(adata, config.qc)
    filtered_path = paths.results / "adata_filtered.h5ad"
    adata.write_h5ad(filtered_path)

    adata, preprocess_notes = preprocess_and_integrate(adata, config)
    preprocessed_path = paths.results / "adata_preprocessed.h5ad"
    adata.write_h5ad(preprocessed_path)

    adata, annotation_summary, annotation_notes = cluster_and_annotate(adata, config)
    annotated_path = paths.results / "adata_annotated.h5ad"
    adata.write_h5ad(annotated_path)

    sample_level, patient_level = compute_cell_type_composition(adata, config)
    sample_level.to_csv(paths.results / "cell_type_composition_sample_level.csv", index=False)
    patient_level.to_csv(paths.results / "cell_type_composition_pair_level.csv", index=False)

    da_df = compute_differential_abundance(patient_level, config)
    da_df.to_csv(paths.results / "differential_abundance.csv", index=False)

    pseudobulk_df = compute_pseudobulk_profiles(adata, config)
    pseudobulk_df.to_csv(paths.results / "pseudobulk_profiles.csv", index=False)
    annotation_summary.to_csv(paths.results / "annotation_summary.csv", index=False)

    figure_outputs = {
        "qc_summary_raw": qc_figures,
        "umap_condition": plot_umap_condition(adata, config, paths.figures),
        "umap_annotation": plot_umap_annotation(adata, config, paths.figures),
        "cell_type_composition": plot_composition(patient_level, config, paths.figures),
        "differential_abundance": plot_differential_abundance(da_df, paths.figures),
    }

    run_summary = {
        "preprocess_notes": preprocess_notes,
        "annotation_notes": annotation_notes,
        "n_cells_raw": n_cells_raw,
        "n_genes_raw": n_genes_raw,
        "n_cells_filtered": int(adata.n_obs),
        "n_genes_filtered": int(adata.n_vars),
        "figure_outputs": figure_outputs,
    }
    _save_json(paths.logs / "core_run_summary.json", run_summary)

    core_outputs = {
        "adata_raw_with_qc": str(raw_path),
        "adata_filtered": str(filtered_path),
        "adata_preprocessed": str(preprocessed_path),
        "adata_annotated": str(annotated_path),
        "cell_type_composition_sample_level": str(paths.results / "cell_type_composition_sample_level.csv"),
        "cell_type_composition_pair_level": str(paths.results / "cell_type_composition_pair_level.csv"),
        "differential_abundance": str(paths.results / "differential_abundance.csv"),
        "pseudobulk_profiles": str(paths.results / "pseudobulk_profiles.csv"),
        "annotation_summary": str(paths.results / "annotation_summary.csv"),
        "core_run_summary": str(paths.logs / "core_run_summary.json"),
    }
    core_outputs.update(build_report(config, manifest, paths, core_outputs))
    return core_outputs


def run_advanced_domains(
    project: str | Path,
    manifest_csv: str | Path,
    workdir: str | Path,
    domains: list[str] | None = None,
) -> list[dict]:
    """Run callable analysis domains against an annotated object."""
    config = WorkflowConfig.from_yaml(project)
    manifest = ManifestTable.from_csv(manifest_csv)
    manifest.validate_conditions(config)
    paths = RunPaths.from_config(Path(workdir), config.outputs)
    annotated_path = paths.results / "adata_annotated.h5ad"
    if not annotated_path.exists():
        raise FileNotFoundError(f"Expected annotated checkpoint at {annotated_path}")

    adata = sc.read_h5ad(annotated_path)
    context = DomainContext(config=config, manifest=manifest, paths=paths, adata=adata)
    registry = get_domain_registry()
    selected = domains or config.domains.enabled

    results = []
    for name in selected:
        if name not in registry:
            results.append(
                {
                    "name": name,
                    "status": "unknown",
                    "message": "Domain is not registered.",
                    "outputs": [],
                }
            )
            continue
        result = registry[name](context)
        results.append(result.model_dump())

    advanced_path = paths.results / "adata_advanced.h5ad"
    context.adata.write_h5ad(advanced_path)
    _save_json(paths.logs / "advanced_domain_summary.json", {"results": results})

    core_outputs = {
        "adata_annotated": str(annotated_path),
        "adata_advanced": str(advanced_path),
        "advanced_domain_summary": str(paths.logs / "advanced_domain_summary.json"),
    }
    build_report(config, manifest, paths, core_outputs, domain_results=results)
    return results

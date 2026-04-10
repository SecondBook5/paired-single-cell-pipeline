"""Public API for paired_sc."""

from .config.models import WorkflowConfig
from .domains import (
    run_cell_cycle,
    run_cellot,
    run_differential_expression,
    run_integration_quality,
    run_latent,
    run_liana,
    run_magic,
    run_pathway_activity,
    run_pathway_enrichment,
    run_pseudobulk_de,
    run_regulatory,
    run_robustness,
    run_target_population,
    run_target_subclustering,
    run_trajectory,
)
from .io.manifest import ManifestTable
from .report.core import build_report
from .workflow import run_advanced_domains, run_core_pipeline

__all__ = [
    "WorkflowConfig",
    "ManifestTable",
    "build_report",
    "run_advanced_domains",
    "run_core_pipeline",
    "run_cell_cycle",
    "run_differential_expression",
    "run_pathway_enrichment",
    "run_liana",
    "run_magic",
    "run_trajectory",
    "run_latent",
    "run_regulatory",
    "run_cellot",
    "run_pseudobulk_de",
    "run_pathway_activity",
    "run_robustness",
    "run_integration_quality",
    "run_target_population",
    "run_target_subclustering",
]

__version__ = "0.1.0"


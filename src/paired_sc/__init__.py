"""Public API for paired_sc."""

from .config.models import WorkflowConfig
from .domains import run_latent, run_liana, run_magic, run_regulatory, run_trajectory
from .io.manifest import ManifestTable
from .report.core import build_report
from .workflow import run_advanced_domains, run_core_pipeline

__all__ = [
    "WorkflowConfig",
    "ManifestTable",
    "build_report",
    "run_advanced_domains",
    "run_core_pipeline",
    "run_liana",
    "run_magic",
    "run_trajectory",
    "run_latent",
    "run_regulatory",
]

__version__ = "0.1.0"


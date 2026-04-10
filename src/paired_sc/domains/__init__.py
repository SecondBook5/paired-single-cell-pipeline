"""Callable advanced domains."""

from .cell_cycle import run as run_cell_cycle
from .cellot import run as run_cellot
from .differential_expression import run as run_differential_expression
from .integration_quality import run as run_integration_quality
from .latent import run as run_latent
from .liana import run as run_liana
from .magic import run as run_magic
from .pathway_activity import run as run_pathway_activity
from .pathway_enrichment import run as run_pathway_enrichment
from .pseudobulk_de import run as run_pseudobulk_de
from .regulatory import run as run_regulatory
from .robustness import run as run_robustness
from .target_population import run as run_target_population
from .target_subclustering import run as run_target_subclustering
from .registry import get_domain_registry, get_domain_status
from .trajectory import run as run_trajectory

__all__ = [
    "get_domain_registry",
    "get_domain_status",
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

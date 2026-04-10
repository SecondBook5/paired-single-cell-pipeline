"""Domain registry and availability status."""

from __future__ import annotations

from .base import dependency_available
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
from .trajectory import run as run_trajectory


DOMAIN_DEPENDENCIES = {
    "cell_cycle": [],
    "differential_expression": [],
    "pathway_enrichment": ["gseapy"],
    "liana": ["liana"],
    "magic": ["magic"],
    "trajectory": [],
    "latent": ["scvi"],
    "regulatory": ["pyscenic"],
    "cellot": ["ot"],
    "pseudobulk_de": [],
    "pathway_activity": ["gseapy"],
    "robustness": [],
    "integration_quality": ["sklearn"],
    "target_population": [],
    "target_subclustering": [],
}


def get_domain_registry():
    return {
        "cell_cycle": run_cell_cycle,
        "differential_expression": run_differential_expression,
        "pathway_enrichment": run_pathway_enrichment,
        "liana": run_liana,
        "magic": run_magic,
        "trajectory": run_trajectory,
        "latent": run_latent,
        "regulatory": run_regulatory,
        "cellot": run_cellot,
        "pseudobulk_de": run_pseudobulk_de,
        "pathway_activity": run_pathway_activity,
        "robustness": run_robustness,
        "integration_quality": run_integration_quality,
        "target_population": run_target_population,
        "target_subclustering": run_target_subclustering,
    }


def get_domain_status() -> list[dict]:
    rows = []
    for name, deps in DOMAIN_DEPENDENCIES.items():
        missing = [dep for dep in deps if not dependency_available(dep)]
        rows.append(
            {
                "name": name,
                "available": len(missing) == 0,
                "missing_dependencies": missing,
            }
        )
    return rows

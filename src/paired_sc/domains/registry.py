"""Domain registry and availability status."""

from __future__ import annotations

from .base import dependency_available
from .latent import run as run_latent
from .liana import run as run_liana
from .magic import run as run_magic
from .regulatory import run as run_regulatory
from .trajectory import run as run_trajectory


DOMAIN_DEPENDENCIES = {
    "liana": ["liana"],
    "magic": ["magic"],
    "trajectory": [],
    "latent": ["scvi"],
    "regulatory": ["pyscenic"],
}


def get_domain_registry():
    return {
        "liana": run_liana,
        "magic": run_magic,
        "trajectory": run_trajectory,
        "latent": run_latent,
        "regulatory": run_regulatory,
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


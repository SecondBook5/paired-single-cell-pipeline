"""Callable advanced domains."""

from .latent import run as run_latent
from .liana import run as run_liana
from .magic import run as run_magic
from .regulatory import run as run_regulatory
from .registry import get_domain_registry, get_domain_status
from .trajectory import run as run_trajectory

__all__ = [
    "get_domain_registry",
    "get_domain_status",
    "run_liana",
    "run_magic",
    "run_trajectory",
    "run_latent",
    "run_regulatory",
]

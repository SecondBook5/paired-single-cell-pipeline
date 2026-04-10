"""Statistical summaries."""

from .paired import (
    compute_cell_type_composition,
    compute_differential_abundance,
    compute_pseudobulk_profiles,
)

__all__ = [
    "compute_cell_type_composition",
    "compute_differential_abundance",
    "compute_pseudobulk_profiles",
]


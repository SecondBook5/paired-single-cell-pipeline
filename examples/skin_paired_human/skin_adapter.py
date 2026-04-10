"""Example CellTypist label adapter for paired human skin studies."""

SKIN_CELLTYPIST_MAP = {
    "Undifferentiated_KC": "Keratinocytes",
    "Differentiated_KC": "Keratinocytes",
    "F1": "Fibroblasts",
    "F2": "Fibroblasts",
    "F3": "Fibroblasts",
    "VE1": "Endothelial",
    "VE2": "Endothelial",
    "VE3": "Endothelial",
    "LE1": "Endothelial",
    "LE2": "Endothelial",
    "Pericyte_1": "Pericytes",
    "Pericyte_2": "Pericytes",
    "Melanocyte": "Melanocytes",
    "Schwann_1": "Schwann cells",
    "Schwann_2": "Schwann cells",
    "Th": "CD4+ T cells",
    "Tc": "CD8+ T cells",
    "Treg": "Regulatory T cells",
    "ILC1_3": "NK cells",
    "ILC1_NK": "NK cells",
    "ILC2": "NK cells",
    "NK": "NK cells",
    "Macro_1": "Macrophages",
    "Macro_2": "Macrophages",
    "Inf_mac": "Macrophages",
    "Mono_mac": "Monocytes",
    "DC1": "Dendritic cells",
    "DC2": "Dendritic cells",
    "MigDC": "Dendritic cells",
    "moDC": "Dendritic cells",
    "LC": "Dendritic cells",
    "migLC": "Dendritic cells",
    "Mast_cell": "Mast cells",
    "Plasma": "Plasma cells",
}


def remap_labels(labels: list[str]) -> list[str]:
    """Map CellTypist fine labels to a higher-level skin schema."""
    return [SKIN_CELLTYPIST_MAP.get(label, label) for label in labels]


"""Quality control helpers."""

from __future__ import annotations

import scanpy as sc


def compute_qc_metrics(adata: sc.AnnData) -> sc.AnnData:
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]", regex=True)
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    return adata


def filter_qc(adata: sc.AnnData, qc_config) -> sc.AnnData:
    keep = (
        (adata.obs["n_genes_by_counts"] >= qc_config.min_genes)
        & (adata.obs["total_counts"] >= qc_config.min_counts)
        & (adata.obs["pct_counts_mt"] <= qc_config.max_pct_mt)
    )
    if qc_config.max_counts is not None:
        keep &= adata.obs["total_counts"] <= qc_config.max_counts
    out = adata[keep, :].copy()
    sc.pp.filter_genes(out, min_cells=qc_config.min_cells_per_gene)
    return out


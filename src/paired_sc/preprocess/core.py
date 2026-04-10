"""Normalization, PCA, Harmony integration, and manifold construction."""

from __future__ import annotations

import numpy as np
import scanpy as sc


def preprocess_and_integrate(adata: sc.AnnData, config) -> tuple[sc.AnnData, list[str]]:
    notes: list[str] = []

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=config.preprocess.target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata

    try:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=min(config.preprocess.n_top_genes, adata.n_vars),
            batch_key=config.batch_key if adata.obs[config.batch_key].nunique() > 1 else None,
            flavor=config.preprocess.hvg_flavor,
            layer="counts" if config.preprocess.hvg_flavor == "seurat_v3" and "counts" in adata.layers else None,
            subset=False,
        )
        if "highly_variable" not in adata.var or not bool(adata.var["highly_variable"].fillna(False).any()):
            adata.var["highly_variable"] = True
            notes.append("Highly variable gene selection returned no genes; all genes were retained.")
    except Exception as exc:
        adata.var["highly_variable"] = True
        notes.append(f"Highly variable gene selection failed; all genes were retained instead ({exc}).")
    sc.pp.scale(adata, max_value=config.preprocess.scale_max_value, zero_center=True)
    n_comps = max(1, min(config.preprocess.pca_n_comps, adata.n_vars - 1, adata.n_obs - 1))
    sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack")

    try:
        import harmonypy

        harmony_out = harmonypy.run_harmony(
            adata.obsm["X_pca"],
            adata.obs,
            config.batch_key,
            max_iter_harmony=config.preprocess.harmony_max_iter,
            verbose=False,
        )
        adata.obsm["X_pca_harmony"] = np.asarray(harmony_out.Z_corr)
    except Exception as exc:  # pragma: no cover - fallback path
        adata.obsm["X_pca_harmony"] = np.asarray(adata.obsm["X_pca"])
        notes.append(f"Harmony unavailable; using raw PCA representation instead ({exc}).")

    sc.pp.neighbors(
        adata,
        n_neighbors=config.preprocess.neighbors_n_neighbors,
        n_pcs=max(1, min(config.preprocess.neighbors_n_pcs, n_comps)),
        use_rep="X_pca_harmony",
        method="umap",
        metric="euclidean",
    )
    sc.tl.umap(
        adata,
        min_dist=config.preprocess.umap_min_dist,
        spread=config.preprocess.umap_spread,
    )
    return adata, notes

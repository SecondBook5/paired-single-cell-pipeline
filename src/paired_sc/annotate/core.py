"""Leiden clustering and configurable annotation."""

from __future__ import annotations

import pandas as pd
import scanpy as sc

from .adapters import load_label_adapter


def cluster_and_annotate(adata: sc.AnnData, config) -> tuple[sc.AnnData, pd.DataFrame, list[str]]:
    notes: list[str] = []

    for resolution in config.annotation.leiden_resolutions:
        sc.tl.leiden(adata, resolution=resolution, key_added=f"leiden_res{resolution}")

    if config.annotation.backend == "none":
        primary_key = config.annotation.primary_leiden_key
        adata.obs[config.annotation.cell_type_key] = adata.obs[primary_key].astype(str)
        summary = (
            adata.obs[config.annotation.cell_type_key]
            .value_counts()
            .rename_axis(config.annotation.cell_type_key)
            .reset_index(name="n_cells")
        )
        notes.append("Annotation backend set to 'none'; primary Leiden clusters were used as cell labels.")
        return adata, summary, notes

    from celltypist import annotate as celltypist_annotate
    from celltypist import models

    models.Model.load(config.annotation.model)
    annotation_input = adata.raw.to_adata() if adata.raw is not None else adata.copy()
    annotation_input.var_names_make_unique()
    annotation_input.obs = adata.obs.copy()
    predictions = celltypist_annotate(
        annotation_input,
        model=config.annotation.model,
        majority_voting=True,
        over_clustering=config.annotation.primary_leiden_key,
    )
    pred_df = predictions.predicted_labels.copy()
    if not pred_df.index.equals(adata.obs_names):
        pred_df = pred_df.reindex(adata.obs_names)
    raw_labels = pred_df["majority_voting"].astype(str).tolist()

    adapter = load_label_adapter(config.resolve_optional_path(config.annotation.label_map_path))
    mapped = adapter(raw_labels) if adapter is not None else raw_labels
    adata.obs["predicted_labels"] = pred_df["predicted_labels"].astype(str).values
    adata.obs["majority_voting"] = raw_labels
    adata.obs[config.annotation.cell_type_key] = mapped

    summary = (
        adata.obs[config.annotation.cell_type_key]
        .value_counts()
        .rename_axis(config.annotation.cell_type_key)
        .reset_index(name="n_cells")
    )
    return adata, summary, notes

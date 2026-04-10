"""Target-population overview domain."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
from scipy import sparse

from ..plotting.core import save_figure_bundle
from .base import DomainContext, DomainResult


def _expr_vector(adata, gene: str) -> np.ndarray | None:
    if gene not in adata.var_names:
        return None
    x = adata[:, gene].X
    if sparse.issparse(x):
        return x.toarray().ravel()
    return np.asarray(x).ravel()


def run(context: DomainContext) -> DomainResult:
    name = "target_population"
    outdir = context.domain_dir(name)
    groupby = context.config.domains.target_groupby or context.config.annotation.cell_type_key
    target = context.config.domains.target_group
    if groupby not in context.adata.obs.columns or not target:
        return DomainResult(name=name, status="skipped", message="Configure domains.target_groupby and domains.target_group to run target-population analysis.")

    mask = context.adata.obs[groupby].astype(str).eq(target).to_numpy()
    if not mask.any():
        return DomainResult(name=name, status="skipped", message=f"No cells matched {groupby}={target!r}.")
    target_adata = context.adata[mask].copy()
    outputs: list[str] = []

    total_by_condition = (
        context.adata.obs.groupby(context.config.condition_key, observed=False)
        .size()
        .rename("total_cells")
        .reset_index()
    )
    target_by_condition = (
        target_adata.obs.groupby(context.config.condition_key, observed=False)
        .size()
        .rename("target_cells")
        .reset_index()
    )
    freq_df = total_by_condition.merge(target_by_condition, on=context.config.condition_key, how="left")
    freq_df["target_cells"] = freq_df["target_cells"].fillna(0).astype(int)
    freq_df["pct_cells"] = 100 * freq_df["target_cells"] / freq_df["total_cells"]
    freq_df.to_csv(outdir / "target_population_frequency.csv", index=False)
    outputs.append(str(outdir / "target_population_frequency.csv"))

    marker_rows = []
    for gene in context.config.domains.target_marker_genes:
        expr = _expr_vector(target_adata, gene)
        if expr is None:
            continue
        tmp = target_adata.obs[[context.config.condition_key]].copy()
        tmp["expr"] = expr
        summary = (
            tmp.groupby(context.config.condition_key, observed=False)["expr"]
            .agg(mean="mean", pct_positive=lambda x: 100 * np.mean(np.asarray(x) > 0))
            .reset_index()
        )
        summary["gene"] = gene
        marker_rows.append(summary)
    if marker_rows:
        marker_df = pd.concat(marker_rows, ignore_index=True)
        marker_df.to_csv(outdir / "target_population_marker_summary.csv", index=False)
        outputs.append(str(outdir / "target_population_marker_summary.csv"))

    if context.config.domains.target_score_genes:
        genes = [gene for gene in context.config.domains.target_score_genes if gene in target_adata.var_names]
        if genes:
            sc.tl.score_genes(target_adata, gene_list=genes, score_name="target_program_score")
            score_df = (
                target_adata.obs.groupby(context.config.condition_key, observed=False)["target_program_score"]
                .agg(["mean", "median"])
                .reset_index()
            )
            score_df.to_csv(outdir / "target_population_program_score.csv", index=False)
            outputs.append(str(outdir / "target_population_program_score.csv"))

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0))
    sns.barplot(data=freq_df, x=context.config.condition_key, y="pct_cells", order=context.config.condition_order, ax=axes[0], color="#4C78A8")
    axes[0].set_title(f"{target} frequency")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Percent of all cells")

    if "X_umap" in context.adata.obsm:
        umap = np.asarray(context.adata.obsm["X_umap"])
        axes[1].scatter(umap[:, 0], umap[:, 1], c="#DDDDDD", s=2, alpha=0.25, linewidths=0)
        axes[1].scatter(umap[mask, 0], umap[mask, 1], c="#C41E3A", s=4, alpha=0.8, linewidths=0)
        axes[1].set_title(f"{target} subset on UMAP")
        axes[1].set_xticks([])
        axes[1].set_yticks([])
    else:
        axes[1].text(0.5, 0.5, "UMAP not available", ha="center", va="center", transform=axes[1].transAxes)
        axes[1].axis("off")
    outputs.extend(save_figure_bundle(fig, outdir / "target_population_overview"))

    return DomainResult(
        name=name,
        status="completed",
        message="Target-population overview completed.",
        outputs=outputs,
        metadata={"groupby": groupby, "target_group": target, "n_target_cells": int(target_adata.n_obs)},
    )

"""Target-population subclustering domain."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns

from ..plotting.core import save_figure_bundle
from .base import DomainContext, DomainResult, dependency_available


def run(context: DomainContext) -> DomainResult:
    name = "target_subclustering"
    outdir = context.domain_dir(name)
    groupby = context.config.domains.target_groupby or context.config.annotation.cell_type_key
    target = context.config.domains.target_group
    if groupby not in context.adata.obs.columns or not target:
        return DomainResult(name=name, status="skipped", message="Configure domains.target_groupby and domains.target_group to run target-population subclustering.")

    mask = context.adata.obs[groupby].astype(str).eq(target).to_numpy()
    if not mask.any():
        return DomainResult(name=name, status="skipped", message=f"No cells matched {groupby}={target!r}.")

    sub = context.adata[mask].copy()
    if sub.n_obs < 50:
        return DomainResult(name=name, status="skipped", message="Too few target cells were available for subclustering.")

    counts_layer = sub.layers["counts"].copy() if "counts" in sub.layers else None
    if counts_layer is not None:
        sc.pp.highly_variable_genes(sub, n_top_genes=context.config.domains.target_subcluster_hvgs, flavor="seurat_v3", layer="counts")
    else:
        sc.pp.highly_variable_genes(sub, n_top_genes=context.config.domains.target_subcluster_hvgs, flavor="seurat_v3")
    hvg = sub[:, sub.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=context.config.preprocess.scale_max_value)
    sc.tl.pca(hvg, n_comps=min(30, hvg.n_vars - 1), svd_solver="arpack")

    rep_key = "X_pca"
    if dependency_available("harmonypy") and context.config.batch_key in hvg.obs.columns:
        try:
            import harmonypy

            harmony = harmonypy.run_harmony(np.asarray(hvg.obsm["X_pca"]), hvg.obs, [context.config.batch_key], verbose=False, random_state=42)
            hvg.obsm["X_pca_harmony"] = np.asarray(harmony.Z_corr)
            rep_key = "X_pca_harmony"
        except Exception:
            rep_key = "X_pca"

    sc.pp.neighbors(hvg, use_rep=rep_key, n_neighbors=context.config.domains.target_subcluster_neighbors)
    sc.tl.umap(hvg, min_dist=0.35)
    sc.tl.leiden(hvg, resolution=context.config.domains.target_subcluster_resolution, key_added="target_subcluster")

    sub.obs["target_subcluster"] = hvg.obs["target_subcluster"].astype(str).values
    sc.tl.rank_genes_groups(sub, "target_subcluster", method="wilcoxon", n_genes=25)
    markers = sc.get.rank_genes_groups_df(sub, None)
    markers.to_csv(outdir / "target_subcluster_markers.csv", index=False)

    context.adata.obs["target_subcluster"] = pd.Series(index=context.adata.obs_names, dtype="object")
    context.adata.obs.loc[sub.obs_names, "target_subcluster"] = sub.obs["target_subcluster"].astype(str).values
    outputs = [str(outdir / "target_subcluster_markers.csv")]

    abundance = (
        sub.obs.groupby([context.config.donor_key, context.config.condition_key, "target_subcluster"], observed=False)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    totals = (
        sub.obs.groupby([context.config.donor_key, context.config.condition_key], observed=False)
        .size()
        .rename("total_target_cells")
        .reset_index()
    )
    abundance = abundance.merge(totals, on=[context.config.donor_key, context.config.condition_key], how="left")
    abundance["pct_within_target"] = 100 * abundance["n_cells"] / abundance["total_target_cells"]
    abundance.to_csv(outdir / "target_subcluster_abundance.csv", index=False)
    outputs.append(str(outdir / "target_subcluster_abundance.csv"))

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2))
    sc.pl.umap(hvg, color="target_subcluster", ax=axes[0], show=False, frameon=False, legend_loc="on data", legend_fontsize=6)
    axes[0].set_title(f"{target} subclusters")

    plot_df = abundance.pivot_table(
        index="target_subcluster",
        columns=context.config.condition_key,
        values="pct_within_target",
        aggfunc="mean",
        observed=False,
    ).fillna(0)
    sns.heatmap(plot_df, cmap="RdBu_r", center=plot_df.mean().mean(), ax=axes[1], linewidths=0.2, linecolor="#EEEEEE")
    axes[1].set_title("Mean abundance by condition")
    axes[1].set_xlabel("Condition")
    axes[1].set_ylabel("Subcluster")
    outputs.extend(save_figure_bundle(fig, outdir / "target_subcluster_summary"))

    return DomainResult(
        name=name,
        status="completed",
        message="Target-population subclustering completed.",
        outputs=outputs,
        metadata={"target_group": target, "n_subclusters": int(sub.obs["target_subcluster"].nunique())},
    )

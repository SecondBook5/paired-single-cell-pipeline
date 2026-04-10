"""Trajectory and topology domain."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy import sparse

from ..plotting.core import save_figure_bundle
from .base import DomainContext, DomainResult


def _connectivity_table(adata, groupby: str) -> pd.DataFrame:
    groups = adata.obs[groupby].astype("category").cat.categories.tolist()
    conn = adata.uns["paga"]["connectivities"]
    matrix = conn.toarray() if sparse.issparse(conn) else np.asarray(conn)
    return pd.DataFrame(matrix, index=groups, columns=groups)


def run(context: DomainContext) -> DomainResult:
    name = "trajectory"
    outdir = context.domain_dir(name)
    groupby = context.config.domains.trajectory_groupby or context.config.annotation.cell_type_key

    if groupby not in context.adata.obs.columns:
        return DomainResult(name=name, status="skipped", message=f"Grouping key {groupby!r} is not present in adata.obs.")

    outputs: list[str] = []
    counts = (
        context.adata.obs[groupby]
        .astype(str)
        .value_counts()
        .rename_axis(groupby)
        .reset_index(name="n_cells")
    )

    def _save_group_sizes(note: str) -> DomainResult:
        counts_path = outdir / "trajectory_group_sizes.csv"
        counts.to_csv(counts_path, index=False)
        outputs.append(str(counts_path))

        fig, ax = plt.subplots(figsize=(5.0, 3.6))
        sns.barplot(data=counts, x=groupby, y="n_cells", color="#4C78A8", ax=ax)
        ax.set_title(f"{groupby} group sizes")
        ax.set_xlabel(groupby)
        ax.set_ylabel("Cells")
        outputs.extend(save_figure_bundle(fig, outdir / "trajectory_group_sizes"))
        return DomainResult(
            name=name,
            status="completed",
            message=note,
            outputs=outputs,
            metadata={"groupby": groupby},
        )

    n_groups = int(counts[groupby].nunique())
    if n_groups < 2:
        return _save_group_sizes(
            f"Only {n_groups} {groupby} group was available; saved group-size summaries instead of PAGA topology."
        )

    try:
        try:
            sc.tl.paga(context.adata, groups=groupby)
            conn_df = _connectivity_table(context.adata, groupby)
        except Exception as exc:
            return _save_group_sizes(
                f"PAGA topology could not be constructed for {groupby}; saved group-size summaries instead ({exc})."
            )
        conn_path = outdir / "paga_connectivities.csv"
        conn_df.to_csv(conn_path)
        outputs.append(str(conn_path))

        fig_w = max(5.2, len(conn_df.columns) * 0.42 + 1.8)
        fig_h = max(4.4, len(conn_df.index) * 0.36 + 1.6)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        sns.heatmap(
            conn_df,
            cmap="YlOrRd",
            linewidths=0.3,
            linecolor="#EEEEEE",
            cbar_kws={"label": "PAGA connectivity"},
            ax=ax,
        )
        ax.set_title(f"PAGA connectivity by {groupby}")
        ax.set_xlabel(groupby)
        ax.set_ylabel(groupby)
        outputs.extend(save_figure_bundle(fig, outdir / "paga_connectivity"))

        root_group = context.config.domains.trajectory_root_group
        if root_group is not None:
            mask = context.adata.obs[groupby].astype(str).eq(root_group).to_numpy()
            if mask.any():
                sc.tl.diffmap(context.adata)
                context.adata.uns["iroot"] = int(np.flatnonzero(mask)[0])
                sc.tl.dpt(context.adata)
                pseudotime = context.adata.obs[["dpt_pseudotime"]].copy()
                pseudo_path = outdir / "dpt_pseudotime.csv"
                pseudotime.to_csv(pseudo_path)
                outputs.append(str(pseudo_path))

                if "X_umap" in context.adata.obsm:
                    fig, ax = plt.subplots(figsize=(5.2, 4.4))
                    sc.pl.umap(
                        context.adata,
                        color="dpt_pseudotime",
                        cmap="viridis",
                        ax=ax,
                        show=False,
                        frameon=False,
                    )
                    ax.set_title("Diffusion pseudotime")
                    outputs.extend(save_figure_bundle(fig, outdir / "dpt_umap"))

        return DomainResult(
            name=name,
            status="completed",
            message="Trajectory domain completed.",
            outputs=outputs,
            metadata={"groupby": groupby, "root_group": root_group},
        )
    except Exception as exc:
        return DomainResult(name=name, status="failed", message=f"Trajectory domain failed ({exc}).")

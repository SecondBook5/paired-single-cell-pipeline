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
from .base import DomainContext, DomainResult, dependency_available


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

        if context.config.domains.trajectory_enable_palantir and dependency_available("palantir") and "X_umap" in context.adata.obsm:
            try:
                import palantir

                rep_key = "X_scvi" if "X_scvi" in context.adata.obsm else "X_pca_harmony" if "X_pca_harmony" in context.adata.obsm else "X_pca"
                rep = np.asarray(context.adata.obsm[rep_key])
                max_cells = min(context.config.domains.trajectory_max_cells, context.adata.n_obs)
                if context.adata.n_obs > max_cells:
                    rng = np.random.default_rng(42)
                    if root_group is not None:
                        root_mask = context.adata.obs[groupby].astype(str).eq(root_group).to_numpy()
                    else:
                        root_mask = np.zeros(context.adata.n_obs, dtype=bool)
                    root_idx = np.flatnonzero(root_mask)
                    keep_root = root_idx[: min(200, len(root_idx))]
                    other_idx = np.setdiff1d(np.arange(context.adata.n_obs), keep_root, assume_unique=False)
                    n_remaining = max_cells - len(keep_root)
                    sampled = rng.choice(other_idx, size=min(n_remaining, len(other_idx)), replace=False) if n_remaining > 0 else np.array([], dtype=int)
                    idx = np.unique(np.concatenate([keep_root, sampled]))
                    adata_pal = context.adata[idx].copy()
                    rep = rep[idx]
                else:
                    adata_pal = context.adata.copy()
                rep_df = pd.DataFrame(rep, index=adata_pal.obs_names)
                dm_res = palantir.utils.run_diffusion_maps(rep_df, n_components=min(10, rep_df.shape[1]))
                ms_data = palantir.utils.determine_multiscale_space(dm_res)
                if root_group is not None and groupby in adata_pal.obs.columns and adata_pal.obs[groupby].astype(str).eq(root_group).any():
                    start_cell = adata_pal.obs_names[adata_pal.obs[groupby].astype(str).eq(root_group).to_numpy()][0]
                else:
                    start_cell = adata_pal.obs_names[0]
                pr_res = palantir.core.run_palantir(ms_data, early_cell=start_cell)
                pt = pd.Series(pr_res.pseudotime.values, index=adata_pal.obs_names)
                entropy = pd.Series(pr_res.entropy.values, index=adata_pal.obs_names)
                context.adata.obs["palantir_pseudotime"] = pt.reindex(context.adata.obs_names)
                context.adata.obs["palantir_entropy"] = entropy.reindex(context.adata.obs_names)
                pal_df = pd.DataFrame(
                    {
                        "palantir_pseudotime": context.adata.obs["palantir_pseudotime"],
                        "palantir_entropy": context.adata.obs["palantir_entropy"],
                    },
                    index=context.adata.obs_names,
                )
                pal_df.to_csv(outdir / "palantir_metrics.csv")
                outputs.append(str(outdir / "palantir_metrics.csv"))
                fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2))
                sc.pl.umap(context.adata, color="palantir_pseudotime", cmap="viridis", ax=axes[0], show=False, frameon=False)
                axes[0].set_title("Palantir pseudotime")
                sc.pl.umap(context.adata, color="palantir_entropy", cmap="magma", ax=axes[1], show=False, frameon=False)
                axes[1].set_title("Palantir entropy")
                outputs.extend(save_figure_bundle(fig, outdir / "palantir_umap"))
            except Exception as exc:
                fail_path = outdir / "palantir_failed.txt"
                fail_path.write_text(str(exc), encoding="utf-8")
                outputs.append(str(fail_path))

        if context.config.domains.trajectory_enable_cellrank and dependency_available("cellrank"):
            try:
                from cellrank.kernels import ConnectivityKernel, PseudotimeKernel
                from cellrank.estimators import GPCCA

                adata_cr = context.adata.copy()
                rep_key = "X_scvi" if "X_scvi" in adata_cr.obsm else "X_pca_harmony" if "X_pca_harmony" in adata_cr.obsm else "X_pca"
                sc.pp.neighbors(
                    adata_cr,
                    use_rep=rep_key,
                    n_neighbors=context.config.preprocess.neighbors_n_neighbors,
                    key_added="neighbors_cellrank",
                )
                ck = ConnectivityKernel(adata_cr, conn_key="neighbors_cellrank_connectivities").compute_transition_matrix()
                if "palantir_pseudotime" in adata_cr.obs.columns and adata_cr.obs["palantir_pseudotime"].notna().any():
                    pk = PseudotimeKernel(adata_cr, time_key="palantir_pseudotime").compute_transition_matrix()
                    kernel = 0.8 * pk + 0.2 * ck
                else:
                    kernel = ck
                estimator = GPCCA(kernel)
                estimator.compute_schur(n_components=min(20, max(3, n_groups)))
                estimator.compute_macrostates(n_states=min(max(3, n_groups), 8), cluster_key=groupby)
                estimator.predict_terminal_states()
                estimator.compute_fate_probabilities()
                fate = estimator.fate_probabilities
                fate_names = list(getattr(fate, "names", [f"state_{idx}" for idx in range(fate.shape[1])]))
                fate_df = pd.DataFrame(fate.X, index=adata_cr.obs_names, columns=fate_names)
                fate_df.to_csv(outdir / "cellrank_fate_probabilities.csv")
                outputs.append(str(outdir / "cellrank_fate_probabilities.csv"))
                max_prob = fate_df.max(axis=1).reindex(context.adata.obs_names)
                context.adata.obs["cellrank_max_fate_probability"] = max_prob
                if "X_umap" in context.adata.obsm:
                    fig, ax = plt.subplots(figsize=(5.2, 4.4))
                    sc.pl.umap(context.adata, color="cellrank_max_fate_probability", cmap="RdYlBu_r", ax=ax, show=False, frameon=False)
                    ax.set_title("CellRank max fate probability")
                    outputs.extend(save_figure_bundle(fig, outdir / "cellrank_umap"))
            except Exception as exc:
                fail_path = outdir / "cellrank_failed.txt"
                fail_path.write_text(str(exc), encoding="utf-8")
                outputs.append(str(fail_path))

        return DomainResult(
            name=name,
            status="completed",
            message="Trajectory domain completed.",
            outputs=outputs,
            metadata={"groupby": groupby, "root_group": root_group},
        )
    except Exception as exc:
        return DomainResult(name=name, status="failed", message=f"Trajectory domain failed ({exc}).")

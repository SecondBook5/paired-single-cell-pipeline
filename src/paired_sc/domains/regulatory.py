"""Regulatory screening domain."""

from __future__ import annotations

import json
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import sparse

from ..plotting.core import save_figure_bundle
from .base import DomainContext, DomainResult, dependency_available


def _resolve_gene_matrix(adata):
    source = adata.raw.to_adata() if adata.raw is not None else adata
    matrix = source.X.toarray() if sparse.issparse(source.X) else source.X
    return pd.DataFrame(matrix, index=adata.obs_names, columns=source.var_names)


def run(context: DomainContext) -> DomainResult:
    name = "regulatory"
    outdir = context.domain_dir(name)

    tf_list_path = context.config.resolve_optional_path(context.config.domains.regulatory_tf_list)
    rankings_path = context.config.resolve_optional_path(context.config.domains.regulatory_rankings)
    motifs_path = context.config.resolve_optional_path(context.config.domains.regulatory_motifs)

    if tf_list_path is None:
        return DomainResult(
            name=name,
            status="skipped",
            message="No regulatory_tf_list was configured for the regulatory domain.",
        )
    if not tf_list_path.exists():
        return DomainResult(name=name, status="skipped", message=f"TF list was not found at {tf_list_path}.")

    with tf_list_path.open("r", encoding="utf-8") as handle:
        tf_candidates = [line.strip() for line in handle if line.strip()]

    expr = _resolve_gene_matrix(context.adata)
    tf_genes = [gene for gene in tf_candidates if gene in expr.columns]
    if not tf_genes:
        return DomainResult(name=name, status="skipped", message="No TFs from the configured list were present in the expression matrix.")

    tf_expr = expr[tf_genes]
    per_gene_var = tf_expr.var(axis=0).sort_values(ascending=False)
    top_tfs = per_gene_var.head(context.config.domains.regulatory_top_variable_tfs).index.tolist()
    tf_summary = (
        tf_expr[top_tfs]
        .join(context.adata.obs[[context.config.annotation.cell_type_key, context.config.condition_key]])
        .groupby([context.config.annotation.cell_type_key, context.config.condition_key], observed=False)
        .mean(numeric_only=True)
        .reset_index()
    )

    outputs = []
    summary_path = outdir / "regulatory_tf_summary.csv"
    tf_summary.to_csv(summary_path, index=False)
    outputs.append(str(summary_path))

    notes_path = outdir / "regulatory_domain_metadata.json"
    notes_path.write_text(
        json.dumps(
            {
                "tf_list_path": str(tf_list_path),
                "rankings_path": str(rankings_path) if rankings_path else None,
                "motifs_path": str(motifs_path) if motifs_path else None,
                "pyscenic_available": dependency_available("pyscenic"),
                "n_tf_candidates": len(tf_candidates),
                "n_tf_in_data": len(tf_genes),
                "top_tfs": top_tfs,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    outputs.append(str(notes_path))

    heatmap_df = tf_summary.pivot_table(
        index=context.config.annotation.cell_type_key,
        values=top_tfs,
        aggfunc="mean",
        observed=False,
    )
    fig_w = max(6.0, len(top_tfs) * 0.32 + 1.8)
    fig_h = max(4.0, len(heatmap_df) * 0.32 + 1.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        heatmap_df,
        cmap="RdBu_r",
        center=0,
        linewidths=0.3,
        linecolor="#EEEEEE",
        cbar_kws={"label": "Mean TF expression"},
        ax=ax,
    )
    ax.set_title("Regulatory screening summary")
    ax.set_xlabel("Transcription factor")
    ax.set_ylabel("Annotation")
    outputs.extend(save_figure_bundle(fig, outdir / "regulatory_tf_heatmap"))

    message = "Regulatory screening summary completed."
    if dependency_available("pyscenic") and rankings_path and motifs_path:
        try:
            from arboreto.algo import grnboost2
            from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
            from pyscenic.aucell import aucell
            from pyscenic.prune import prune2df
            from pyscenic.transform import df2regulons
            from pyscenic.utils import modules_from_adjacencies

            expr_full = _resolve_gene_matrix(context.adata)
            if len(expr_full) > context.config.domains.regulatory_max_cells:
                expr_full = expr_full.sample(context.config.domains.regulatory_max_cells, random_state=42)
            tf_for_grn = [gene for gene in tf_candidates if gene in expr_full.columns]
            if tf_for_grn:
                adj = grnboost2(expr_full, tf_names=tf_for_grn, verbose=False, seed=42)
                adj_path = outdir / "scenic_adjacencies.csv"
                adj.to_csv(adj_path, index=False)
                outputs.append(str(adj_path))

                dbs = [RankingDatabase(fname=str(rankings_path), name=rankings_path.stem)]
                modules = list(modules_from_adjacencies(adj, expr_full))
                pruned = prune2df(dbs, modules, str(motifs_path))
                regulons = df2regulons(pruned)
                if regulons:
                    auc = aucell(expr_full, regulons, num_workers=1)
                    auc_path = outdir / "scenic_auc_matrix.csv"
                    auc.to_csv(auc_path)
                    outputs.append(str(auc_path))

                    common = auc.index.intersection(context.adata.obs_names)
                    if len(common):
                        auc_ct = auc.loc[common].groupby(context.adata.obs.loc[common, context.config.annotation.cell_type_key]).mean()
                        top_regs = auc.var().sort_values(ascending=False).head(min(25, auc.shape[1])).index.tolist()
                        auc_plot = auc_ct[top_regs].copy()
                        auc_z = auc_plot.apply(lambda col: (col - col.mean()) / col.std(ddof=0) if col.std(ddof=0) > 0 else col * 0, axis=0).fillna(0)
                        fig_w = max(6.0, len(top_regs) * 0.32 + 1.8)
                        fig_h = max(4.0, len(auc_z) * 0.32 + 1.8)
                        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
                        sns.heatmap(
                            auc_z,
                            cmap="RdBu_r",
                            center=0,
                            linewidths=0.3,
                            linecolor="#EEEEEE",
                            cbar_kws={"label": "Regulon activity (z-score)"},
                            ax=ax,
                        )
                        ax.set_title("pySCENIC regulon activity")
                        ax.set_xlabel("Regulon")
                        ax.set_ylabel("Annotation")
                        outputs.extend(save_figure_bundle(fig, outdir / "scenic_regulon_heatmap"))
                message += " Full pySCENIC execution completed."
            else:
                message += " No configured TFs were present for pySCENIC expansion."
        except Exception as exc:
            message += f" Full pySCENIC execution was attempted but fell back to the screening summary ({exc})."
    else:
        message += " Full pySCENIC execution is available when rankings and motif resources are configured."

    return DomainResult(
        name=name,
        status="completed",
        message=message,
        outputs=outputs,
        metadata={"n_top_tfs": len(top_tfs)},
    )

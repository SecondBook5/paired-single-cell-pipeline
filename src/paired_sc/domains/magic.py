"""MAGIC targeted imputation domain."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse

from ..plotting.core import save_figure_bundle
from .base import DomainContext, DomainResult, dependency_available


def _expression_frame(adata, genes: list[str]) -> pd.DataFrame:
    source = adata.raw.to_adata() if adata.raw is not None else adata
    gene_index = pd.Index(source.var_names)
    keep = [gene for gene in genes if gene in gene_index]
    if not keep:
        return pd.DataFrame(index=adata.obs_names)
    idx = [int(gene_index.get_loc(gene)) for gene in keep]
    matrix = source.X[:, idx]
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return pd.DataFrame(np.asarray(matrix), index=adata.obs_names, columns=keep)


def run(context: DomainContext) -> DomainResult:
    name = "magic"
    outdir = context.domain_dir(name)

    if not dependency_available("magic"):
        return DomainResult(
            name=name,
            status="skipped",
            message="MAGIC is not installed in the active environment.",
        )

    genes = [gene for gene in context.config.domains.magic_genes if gene]
    if not genes:
        return DomainResult(name=name, status="skipped", message="No magic_genes were configured.")

    expr = _expression_frame(context.adata, genes)
    if expr.empty:
        return DomainResult(name=name, status="skipped", message="Configured MAGIC genes were not found in the AnnData object.")

    try:
        import magic as magic_impute

        operator = magic_impute.MAGIC(
            knn=context.config.domains.magic_knn,
            solver=context.config.domains.magic_solver,
            random_state=42,
            verbose=0,
        )
        imputed = operator.fit_transform(expr)
        if not isinstance(imputed, pd.DataFrame):
            imputed = pd.DataFrame(np.asarray(imputed), index=expr.index, columns=expr.columns)

        outputs = []
        imputed_path = outdir / "magic_imputed_selected_genes.csv"
        imputed.to_csv(imputed_path)
        outputs.append(str(imputed_path))

        summary = (
            imputed.join(context.adata.obs[[context.config.condition_key, context.config.annotation.cell_type_key]])
            .groupby([context.config.annotation.cell_type_key, context.config.condition_key])
            .mean(numeric_only=True)
            .reset_index()
        )
        summary_path = outdir / "magic_imputed_summary.csv"
        summary.to_csv(summary_path, index=False)
        outputs.append(str(summary_path))

        heatmap_df = summary.pivot_table(
            index=context.config.annotation.cell_type_key,
            columns=context.config.condition_key,
            values=imputed.columns[0],
            aggfunc="mean",
        )
        fig, ax = plt.subplots(figsize=(5.4, max(3.4, len(heatmap_df) * 0.35 + 1.6)))
        sns.heatmap(
            heatmap_df,
            cmap="viridis",
            linewidths=0.3,
            linecolor="#EEEEEE",
            cbar_kws={"label": f"{imputed.columns[0]} MAGIC-imputed mean"},
            ax=ax,
        )
        ax.set_title(f"MAGIC summary for {imputed.columns[0]}")
        ax.set_xlabel("Condition")
        ax.set_ylabel("Annotation")
        outputs.extend(save_figure_bundle(fig, outdir / "magic_summary"))

        for gene in imputed.columns:
            context.adata.obs[f"{gene}_magic"] = imputed[gene].reindex(context.adata.obs_names).values

        return DomainResult(
            name=name,
            status="completed",
            message="MAGIC targeted imputation completed.",
            outputs=outputs,
            metadata={"genes": list(imputed.columns)},
        )
    except Exception as exc:
        return DomainResult(name=name, status="failed", message=f"MAGIC domain failed ({exc}).")

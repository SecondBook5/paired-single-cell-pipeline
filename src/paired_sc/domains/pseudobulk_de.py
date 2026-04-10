"""Pseudobulk differential-expression domain."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse, stats as scipy_stats

from ..plotting.core import save_figure_bundle
from .base import DomainContext, DomainResult


def _pseudobulk_table(context: DomainContext) -> tuple[pd.DataFrame, list[str]]:
    matrix = context.adata.raw.X if context.adata.raw is not None else (
        context.adata.layers["counts"] if "counts" in context.adata.layers else context.adata.X
    )
    gene_names = context.adata.raw.var_names if context.adata.raw is not None else context.adata.var_names
    cell_type_key = context.config.annotation.cell_type_key
    obs = context.adata.obs[[context.config.donor_key, context.config.condition_key, cell_type_key]].astype(str).copy()
    labels = (
        obs[context.config.donor_key]
        + "__"
        + obs[context.config.condition_key]
        + "__"
        + obs[cell_type_key]
    )
    codes, unique_labels = pd.factorize(labels, sort=False)
    indicator = sparse.coo_matrix(
        (np.ones(len(codes), dtype=np.float32), (codes, np.arange(len(codes)))),
        shape=(len(unique_labels), context.adata.n_obs),
    ).tocsr()
    summed = indicator @ matrix
    summed = summed.toarray() if sparse.issparse(summed) else np.asarray(summed)
    meta = (
        obs.assign(_group=labels)
        .drop_duplicates("_group")
        .set_index("_group")
        .reindex(unique_labels)
        .reset_index(drop=True)
    )
    wide = pd.DataFrame(summed, columns=gene_names)
    wide.insert(0, "cell_type", meta[cell_type_key].values)
    wide.insert(0, "condition", meta[context.config.condition_key].values)
    wide.insert(0, "donor_id", meta[context.config.donor_key].values)
    return wide, list(gene_names)


def run(context: DomainContext) -> DomainResult:
    name = "pseudobulk_de"
    outdir = context.domain_dir(name)
    pb_df, genes = _pseudobulk_table(context)
    outputs: list[str] = []
    results: list[pd.DataFrame] = []
    cell_type_key = context.config.annotation.cell_type_key

    for cell_type, sub in pb_df.groupby("cell_type", observed=False):
        case = sub[sub["condition"].eq(context.config.case_condition)].copy()
        control = sub[sub["condition"].eq(context.config.control_condition)].copy()
        if len(case) < context.config.domains.pseudobulk_min_replicates_per_condition:
            continue
        if len(control) < context.config.domains.pseudobulk_min_replicates_per_condition:
            continue

        gene_mask = ((sub[genes] >= context.config.domains.pseudobulk_min_count).sum(axis=0) >= 2).values
        filtered_genes = [gene for gene, keep in zip(genes, gene_mask) if keep]
        if not filtered_genes:
            continue

        case_log = np.log1p(case[filtered_genes].astype(float))
        control_log = np.log1p(control[filtered_genes].astype(float))
        rows = []
        for gene in filtered_genes:
            statistic, pvalue = scipy_stats.ttest_ind(case_log[gene], control_log[gene], equal_var=False, nan_policy="omit")
            rows.append(
                {
                    "names": gene,
                    "logfoldchanges": float(case_log[gene].mean() - control_log[gene].mean()),
                    "pvals": float(pvalue) if np.isfinite(pvalue) else 1.0,
                    cell_type_key: cell_type,
                    "n_case": int(len(case)),
                    "n_control": int(len(control)),
                }
            )
        result = pd.DataFrame(rows).sort_values("pvals", kind="stable")
        ranks = np.arange(1, len(result) + 1)
        result["pvals_adj"] = np.minimum.accumulate((result["pvals"].to_numpy()[::-1] * len(result) / ranks[::-1]))[::-1]
        result.to_csv(outdir / f"pseudobulk_de_{str(cell_type).replace(' ', '_')}.csv", index=False)
        outputs.append(str(outdir / f"pseudobulk_de_{str(cell_type).replace(' ', '_')}.csv"))
        results.append(result)

    if not results:
        return DomainResult(name=name, status="skipped", message="No annotation groups had enough donor-level pseudobulk replicates.")

    combined = pd.concat(results, ignore_index=True)
    combined.to_csv(outdir / "pseudobulk_de_combined.csv", index=False)
    outputs.append(str(outdir / "pseudobulk_de_combined.csv"))

    top_hits = (
        combined[combined["pvals_adj"].fillna(1.0) < 0.05]
        .groupby(cell_type_key, observed=False)
        .size()
        .rename("n_significant")
        .reset_index()
        .sort_values("n_significant", ascending=False)
    )
    top_hits.to_csv(outdir / "pseudobulk_de_summary.csv", index=False)
    outputs.append(str(outdir / "pseudobulk_de_summary.csv"))

    fig, ax = plt.subplots(figsize=(7.4, max(3.2, 0.42 * len(top_hits))))
    sns.barplot(data=top_hits, x="n_significant", y=cell_type_key, color="#2E8B57", ax=ax)
    ax.set_title("Significant pseudobulk differential-expression hits")
    ax.set_xlabel("Genes with FDR < 0.05")
    ax.set_ylabel("")
    outputs.extend(save_figure_bundle(fig, outdir / "pseudobulk_de_summary"))

    return DomainResult(
        name=name,
        status="completed",
        message="Pseudobulk differential-expression analysis completed.",
        outputs=outputs,
        metadata={"n_annotations_tested": int(top_hits.shape[0])},
    )

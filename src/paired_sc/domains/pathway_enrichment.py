"""Pathway enrichment domain based on differential-expression outputs."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ..plotting.core import save_figure_bundle
from .base import DomainContext, DomainResult, dependency_available


def run(context: DomainContext) -> DomainResult:
    name = "pathway_enrichment"
    outdir = context.domain_dir(name)
    de_path = context.paths.domains / "differential_expression" / "differential_expression_combined.csv"
    if not dependency_available("gseapy"):
        return DomainResult(name=name, status="skipped", message="gseapy is not installed in the active environment.")
    if not de_path.exists():
        return DomainResult(name=name, status="skipped", message="Differential-expression outputs were not found; run the differential_expression domain first.")

    import gseapy as gp

    de_df = pd.read_csv(de_path)
    cell_type_key = context.config.annotation.cell_type_key
    outputs: list[str] = []
    results: list[pd.DataFrame] = []

    for cell_type, sub in de_df.groupby(cell_type_key, observed=False):
        sig = sub[(sub["pvals_adj"].fillna(1.0) < 0.05)].copy()
        if sig.empty:
            continue
        up_genes = sig[sig["logfoldchanges"] > 0]["names"].astype(str).dropna().tolist()[: context.config.domains.pathway_top_genes]
        down_genes = sig[sig["logfoldchanges"] < 0]["names"].astype(str).dropna().tolist()[: context.config.domains.pathway_top_genes]
        for direction, genes in [("up_case", up_genes), ("up_control", down_genes)]:
            if len(genes) < 5:
                continue
            for gene_set in context.config.domains.pathway_gene_sets:
                try:
                    enr = gp.enrichr(gene_list=genes, gene_sets=[gene_set], organism="human", outdir=None, cutoff=1.0)
                    df = enr.results
                except Exception:
                    continue
                if df is None or df.empty:
                    continue
                df = df[df["Adjusted P-value"] < 0.05].copy()
                if df.empty:
                    continue
                df["cell_type"] = cell_type
                df["direction"] = direction
                df["gene_set"] = gene_set
                df["log10p"] = -np.log10(df["Adjusted P-value"].clip(lower=1e-300))
                results.append(df)

    if not results:
        return DomainResult(name=name, status="skipped", message="No significant enrichment results were recovered from the differential-expression signatures.")

    combined = pd.concat(results, ignore_index=True)
    combined.to_csv(outdir / "pathway_enrichment_combined.csv", index=False)
    outputs.append(str(outdir / "pathway_enrichment_combined.csv"))

    summary = (
        combined.groupby(["cell_type", "direction", "gene_set"], observed=False)
        .head(context.config.domains.pathway_top_terms)
        .copy()
    )
    summary.to_csv(outdir / "pathway_enrichment_top_terms.csv", index=False)
    outputs.append(str(outdir / "pathway_enrichment_top_terms.csv"))

    plot_df = (
        summary.groupby(["cell_type", "direction"], observed=False)
        .head(1)
        .loc[:, ["cell_type", "direction", "Term", "log10p"]]
        .copy()
    )
    plot_df["label"] = plot_df["cell_type"].astype(str) + " | " + plot_df["direction"].astype(str)
    fig_h = max(4.0, 0.35 * len(plot_df))
    fig, ax = plt.subplots(figsize=(7.8, fig_h))
    sns.scatterplot(data=plot_df, x="log10p", y="label", size="log10p", hue="log10p", palette="magma", legend=False, ax=ax)
    for _, row in plot_df.iterrows():
        ax.text(float(row["log10p"]) + 0.08, row["label"], str(row["Term"]), va="center", fontsize=7)
    ax.set_title("Top enriched pathways by annotation and direction")
    ax.set_xlabel("-log10 adjusted p-value")
    ax.set_ylabel("")
    outputs.extend(save_figure_bundle(fig, outdir / "pathway_enrichment_summary"))

    return DomainResult(
        name=name,
        status="completed",
        message="Pathway enrichment completed.",
        outputs=outputs,
        metadata={"n_enrichment_rows": int(len(combined))},
    )

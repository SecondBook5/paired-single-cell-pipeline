"""Pathway activity scoring domain using ssGSEA on pseudobulk profiles."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ..plotting.core import save_figure_bundle
from .base import DomainContext, DomainResult, dependency_available


def run(context: DomainContext) -> DomainResult:
    name = "pathway_activity"
    outdir = context.domain_dir(name)
    pb_path = context.paths.results / "pseudobulk_profiles.csv"
    if not dependency_available("gseapy"):
        return DomainResult(name=name, status="skipped", message="gseapy is not installed in the active environment.")
    if not pb_path.exists():
        return DomainResult(name=name, status="skipped", message="Core pseudobulk profiles were not found.")

    import gseapy as gp

    pb = pd.read_csv(pb_path)
    meta_cols = {"donor_id", "condition", "cell_type", "pseudobulk_total_counts"}
    gene_cols = [col for col in pb.columns if col not in meta_cols]
    if not gene_cols:
        return DomainResult(name=name, status="skipped", message="No gene columns were present in the pseudobulk profiles table.")

    expression = pb[gene_cols].T
    expression.columns = (
        pb["donor_id"].astype(str)
        + "__"
        + pb["condition"].astype(str)
        + "__"
        + pb["cell_type"].astype(str)
    ).tolist()
    metadata = pd.DataFrame(
        {
            "sample": expression.columns,
            "donor_id": pb["donor_id"].astype(str).values,
            "condition": pb["condition"].astype(str).values,
            "cell_type": pb["cell_type"].astype(str).values,
        }
    ).set_index("sample")

    try:
        library = gp.get_library(name=context.config.domains.pathway_activity_library, organism="Human")
    except Exception as exc:
        return DomainResult(name=name, status="failed", message=f"Failed to retrieve pathway library {context.config.domains.pathway_activity_library!r} ({exc}).")

    genes_in_data = set(expression.index.astype(str))
    library = {
        key: [gene for gene in genes if gene in genes_in_data]
        for key, genes in library.items()
    }
    library = {
        key: genes
        for key, genes in library.items()
        if len(genes) >= context.config.domains.pathway_activity_min_genes
    }
    if not library:
        return DomainResult(name=name, status="skipped", message="No pathway gene sets overlapped the pseudobulk expression matrix.")

    scores = gp.ssgsea(
        data=expression,
        gene_sets=library,
        outdir=None,
        sample_norm_method="rank",
        no_plot=True,
        processes=1,
    ).res2d.pivot(index="Term", columns="Name", values="NES").fillna(0)

    scores.to_csv(outdir / "pathway_activity_scores.csv")
    metadata.to_csv(outdir / "pathway_activity_metadata.csv")
    outputs = [str(outdir / "pathway_activity_scores.csv"), str(outdir / "pathway_activity_metadata.csv")]

    grouped = []
    for sample in scores.columns:
        donor_id, condition, cell_type = sample.split("__", 2)
        grouped.append({"sample": sample, "condition": condition, "cell_type": cell_type})
    grouped_df = pd.DataFrame(grouped).set_index("sample")
    mean_scores = (
        scores.T.join(grouped_df)
        .groupby(["cell_type", "condition"], observed=False)
        .mean(numeric_only=True)
        .T
    )
    top_terms = mean_scores.var(axis=1).sort_values(ascending=False).head(25).index.tolist()
    plot_df = mean_scores.loc[top_terms]

    fig_w = max(6.5, 1.2 + 0.35 * plot_df.shape[1])
    fig_h = max(5.0, 1.5 + 0.22 * plot_df.shape[0])
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(plot_df, cmap="RdBu_r", center=0, linewidths=0.2, linecolor="#EEEEEE", ax=ax)
    ax.set_title("Pathway activity across donor-aware pseudobulk profiles")
    ax.set_xlabel("Cell type / condition")
    ax.set_ylabel("Pathway")
    outputs.extend(save_figure_bundle(fig, outdir / "pathway_activity_heatmap"))

    return DomainResult(
        name=name,
        status="completed",
        message="Pathway activity scoring completed.",
        outputs=outputs,
        metadata={"n_pathways": int(scores.shape[0]), "n_profiles": int(scores.shape[1])},
    )

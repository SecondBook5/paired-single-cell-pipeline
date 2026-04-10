"""Cell-type differential-expression domain."""

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


def _safe_name(value: str) -> str:
    return value.replace(" ", "_").replace("/", "_").replace("+", "pos")


def _vectorized_stats(adata_ct, result: pd.DataFrame, case_mask, control_mask) -> pd.DataFrame:
    genes = [gene for gene in result["names"].astype(str).tolist() if gene in adata_ct.var_names]
    if not genes:
        return result

    x_case = adata_ct[case_mask, genes].X
    x_ctrl = adata_ct[control_mask, genes].X
    if sparse.issparse(x_case):
        mean_case = np.asarray(x_case.mean(axis=0)).ravel()
        mean_ctrl = np.asarray(x_ctrl.mean(axis=0)).ravel()
        pct_case = np.asarray((x_case > 0).mean(axis=0)).ravel() * 100
        pct_ctrl = np.asarray((x_ctrl > 0).mean(axis=0)).ravel() * 100
    else:
        mean_case = np.asarray(x_case.mean(axis=0)).ravel()
        mean_ctrl = np.asarray(x_ctrl.mean(axis=0)).ravel()
        pct_case = (np.asarray(x_case) > 0).mean(axis=0) * 100
        pct_ctrl = (np.asarray(x_ctrl) > 0).mean(axis=0) * 100

    lookup = {gene: idx for idx, gene in enumerate(genes)}
    out = result.copy()
    out["mean_case"] = out["names"].map(lambda gene: float(mean_case[lookup[gene]]) if gene in lookup else np.nan)
    out["mean_control"] = out["names"].map(lambda gene: float(mean_ctrl[lookup[gene]]) if gene in lookup else np.nan)
    out["pct_case"] = out["names"].map(lambda gene: float(pct_case[lookup[gene]]) if gene in lookup else np.nan)
    out["pct_control"] = out["names"].map(lambda gene: float(pct_ctrl[lookup[gene]]) if gene in lookup else np.nan)
    return out


def run(context: DomainContext) -> DomainResult:
    name = "differential_expression"
    outdir = context.domain_dir(name)
    cell_type_key = context.config.annotation.cell_type_key
    condition_key = context.config.condition_key
    if cell_type_key not in context.adata.obs.columns:
        return DomainResult(name=name, status="skipped", message=f"{cell_type_key!r} was not present in adata.obs.")

    outputs: list[str] = []
    combined_rows: list[pd.DataFrame] = []
    summary_rows: list[dict] = []

    for cell_type in sorted(context.adata.obs[cell_type_key].astype(str).unique().tolist()):
        adata_ct = context.adata[context.adata.obs[cell_type_key].astype(str).eq(cell_type)].copy()
        case_mask = adata_ct.obs[condition_key].astype(str).eq(context.config.case_condition).to_numpy()
        control_mask = adata_ct.obs[condition_key].astype(str).eq(context.config.control_condition).to_numpy()
        n_case = int(case_mask.sum())
        n_control = int(control_mask.sum())
        if min(n_case, n_control) < context.config.domains.de_min_cells_per_condition:
            continue

        sc.tl.rank_genes_groups(
            adata_ct,
            groupby=condition_key,
            groups=[context.config.case_condition],
            reference=context.config.control_condition,
            method="wilcoxon",
            corr_method="benjamini-hochberg",
            use_raw=False,
            key_added="de_domain",
        )
        result = sc.get.rank_genes_groups_df(adata_ct, group=context.config.case_condition, key="de_domain")
        result = _vectorized_stats(adata_ct, result, case_mask, control_mask)
        result[cell_type_key] = cell_type
        result["n_case"] = n_case
        result["n_control"] = n_control
        result.to_csv(outdir / f"de_{_safe_name(cell_type)}.csv", index=False)
        outputs.append(str(outdir / f"de_{_safe_name(cell_type)}.csv"))
        combined_rows.append(result)

        sig = result[result["pvals_adj"].fillna(1.0) < 0.05]
        summary_rows.append(
            {
                cell_type_key: cell_type,
                "n_case": n_case,
                "n_control": n_control,
                "n_significant": int(len(sig)),
                "n_up_case": int((sig["logfoldchanges"] > 0).sum()),
                "n_up_control": int((sig["logfoldchanges"] < 0).sum()),
            }
        )

    if not combined_rows:
        return DomainResult(name=name, status="skipped", message="No annotation groups met the minimum cells-per-condition threshold.")

    combined = pd.concat(combined_rows, ignore_index=True)
    summary = pd.DataFrame(summary_rows).sort_values("n_significant", ascending=False)
    combined.to_csv(outdir / "differential_expression_combined.csv", index=False)
    summary.to_csv(outdir / "differential_expression_summary.csv", index=False)
    outputs.extend(
        [
            str(outdir / "differential_expression_combined.csv"),
            str(outdir / "differential_expression_summary.csv"),
        ]
    )

    fig, ax = plt.subplots(figsize=(7.4, max(3.6, 0.42 * len(summary))))
    sns.barplot(data=summary, x="n_significant", y=cell_type_key, color="#4C78A8", ax=ax)
    ax.set_title("Significant differential-expression hits by annotation")
    ax.set_xlabel("Genes with FDR < 0.05")
    ax.set_ylabel("")
    outputs.extend(save_figure_bundle(fig, outdir / "differential_expression_summary"))

    shared = combined[
        (combined["pvals_adj"].fillna(1.0) < 0.05)
        & (combined["logfoldchanges"].abs() >= context.config.domains.de_logfc_threshold)
    ].copy()
    shared_counts = shared.groupby("names", observed=False)[cell_type_key].nunique().sort_values(ascending=False)
    top_genes = shared_counts.head(30).index.tolist()
    if top_genes and summary.shape[0] > 1:
        fc_df = (
            combined[combined["names"].isin(top_genes)]
            .pivot_table(index="names", columns=cell_type_key, values="logfoldchanges", aggfunc="first", observed=False)
            .fillna(0.0)
        )
        fig_w = max(6.0, 1.0 + 0.45 * len(fc_df.columns))
        fig_h = max(4.0, 1.0 + 0.22 * len(fc_df.index))
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        vmax = float(np.nanpercentile(np.abs(fc_df.to_numpy()), 95)) if fc_df.size else 1.0
        vmax = max(vmax, 1.0)
        sns.heatmap(fc_df, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax, linewidths=0.2, linecolor="#EEEEEE", ax=ax)
        ax.set_title("Shared differential-expression patterns")
        ax.set_xlabel("Annotation")
        ax.set_ylabel("Gene")
        outputs.extend(save_figure_bundle(fig, outdir / "differential_expression_heatmap"))

    return DomainResult(
        name=name,
        status="completed",
        message="Differential-expression analysis completed.",
        outputs=outputs,
        metadata={"n_annotations_tested": int(summary.shape[0])},
    )

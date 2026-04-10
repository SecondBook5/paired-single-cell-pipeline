"""LIANA ligand-receptor domain."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from ..plotting.core import save_figure_bundle
from .base import DomainContext, DomainResult, dependency_available


def _prep_adata(adata):
    if adata.raw is not None:
        out = adata.raw.to_adata()
        out.obs = adata.obs.copy()
    else:
        out = adata.copy()
    out.var_names_make_unique()
    return out


def _interaction_heatmap(df: pd.DataFrame, outdir: Path) -> list[str]:
    sig = df.loc[df["magnitude_rank"] < 0.05].copy() if "magnitude_rank" in df.columns else df.copy()
    if sig.empty and "magnitude_rank" in df.columns:
        sig = df.nsmallest(min(len(df), 50), "magnitude_rank").copy()
    elif sig.empty:
        sig = df.head(50).copy()

    if sig.empty:
        fig, ax = plt.subplots(figsize=(4.8, 3.6))
        ax.text(0.5, 0.5, "No LIANA interactions available", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        return save_figure_bundle(fig, outdir / "interaction_heatmap")

    mat = sig.groupby(["source", "target"]).size().unstack(fill_value=0)
    fig_w = max(4.8, len(mat.columns) * 0.45 + 1.8)
    fig_h = max(4.0, len(mat.index) * 0.38 + 1.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        mat,
        cmap="YlOrRd",
        linewidths=0.3,
        linecolor="#EEEEEE",
        cbar_kws={"label": "Significant interactions"},
        ax=ax,
    )
    ax.set_title("LIANA interaction counts")
    ax.set_xlabel("Target")
    ax.set_ylabel("Source")
    return save_figure_bundle(fig, outdir / "interaction_heatmap")


def run(context: DomainContext) -> DomainResult:
    name = "liana"
    outdir = context.domain_dir(name)

    if not dependency_available("liana"):
        return DomainResult(
            name=name,
            status="skipped",
            message="LIANA is not installed in the active environment.",
        )

    try:
        import liana as li
    except Exception as exc:
        return DomainResult(name=name, status="skipped", message=f"Unable to import LIANA ({exc}).")

    try:
        adata = _prep_adata(context.adata)
        groupby = context.config.annotation.cell_type_key
        li.mt.rank_aggregate(
            adata,
            groupby=groupby,
            expr_prop=context.config.domains.liana_expr_prop,
            min_cells=context.config.domains.liana_min_cells,
            verbose=False,
            use_raw=False,
            n_perms=100,
        )
        combined = adata.uns["liana_res"].copy()
        outputs = []
        combined_path = outdir / "liana_rank_aggregate.csv"
        combined.to_csv(combined_path, index=False)
        outputs.append(str(combined_path))
        outputs.extend(_interaction_heatmap(combined, outdir))

        condition_results = []
        for condition in context.config.condition_order:
            mask = context.adata.obs[context.config.condition_key].astype(str).eq(condition).to_numpy()
            if mask.sum() < context.config.domains.liana_min_cells:
                continue
            sub = _prep_adata(context.adata[mask].copy())
            li.mt.rank_aggregate(
                sub,
                groupby=groupby,
                expr_prop=context.config.domains.liana_expr_prop,
                min_cells=context.config.domains.liana_min_cells,
                verbose=False,
                use_raw=False,
                n_perms=100,
            )
            res = sub.uns["liana_res"].copy()
            res["condition"] = condition
            res_path = outdir / f"liana_rank_aggregate_{condition}.csv"
            res.to_csv(res_path, index=False)
            outputs.append(str(res_path))
            condition_results.append(res)

        if len(condition_results) >= 2:
            combo = pd.concat(condition_results, ignore_index=True)
            combo["interaction_id"] = (
                combo["source"].astype(str)
                + "__"
                + combo["target"].astype(str)
                + "__"
                + combo["ligand_complex"].astype(str)
                + "__"
                + combo["receptor_complex"].astype(str)
            )
            comparison = combo.pivot_table(
                index="interaction_id",
                columns="condition",
                values="magnitude_rank",
                aggfunc="mean",
            ).reset_index()
            if set(context.config.condition_order).issubset(comparison.columns):
                comparison["delta_rank"] = (
                    comparison[context.config.case_condition] - comparison[context.config.control_condition]
                )
                comparison_path = outdir / "liana_condition_comparison.csv"
                comparison.to_csv(comparison_path, index=False)
                outputs.append(str(comparison_path))

        return DomainResult(
            name=name,
            status="completed",
            message="LIANA interaction analysis completed.",
            outputs=outputs,
            metadata={"n_interactions": int(len(combined))},
        )
    except Exception as exc:
        return DomainResult(name=name, status="failed", message=f"LIANA domain failed ({exc}).")

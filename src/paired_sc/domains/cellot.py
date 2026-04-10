"""Optimal-transport domain for paired target populations."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ..plotting.core import save_figure_bundle
from .base import DomainContext, DomainResult, dependency_available


def run(context: DomainContext) -> DomainResult:
    name = "cellot"
    outdir = context.domain_dir(name)
    groupby = context.config.domains.target_groupby or context.config.annotation.cell_type_key
    target = context.config.domains.target_group
    if groupby not in context.adata.obs.columns or not target:
        return DomainResult(name=name, status="skipped", message="Configure domains.target_groupby and domains.target_group to run the CellOT/OT domain.")
    if not dependency_available("ot"):
        return DomainResult(name=name, status="skipped", message="POT (`ot`) is not installed in the active environment.")

    import ot

    target_adata = context.adata[context.adata.obs[groupby].astype(str).eq(target)].copy()
    case_mask = target_adata.obs[context.config.condition_key].astype(str).eq(context.config.case_condition).to_numpy()
    control_mask = target_adata.obs[context.config.condition_key].astype(str).eq(context.config.control_condition).to_numpy()
    if case_mask.sum() < 10 or control_mask.sum() < 10:
        return DomainResult(name=name, status="skipped", message="Too few target-population cells were available in one or both conditions.")

    rep_key = next((key for key in ["X_scvi", "X_pca_harmony", "X_pca"] if key in target_adata.obsm), None)
    if rep_key is None:
        return DomainResult(name=name, status="skipped", message="No compatible latent representation was available for optimal-transport analysis.")

    x_case = np.asarray(target_adata[case_mask].obsm[rep_key])
    x_control = np.asarray(target_adata[control_mask].obsm[rep_key])
    n_sub = min(context.config.domains.cellot_subsample_n, len(x_case), len(x_control))
    if n_sub < 10:
        return DomainResult(name=name, status="skipped", message="Too few target-population cells remained after subsampling for optimal transport.")

    rng = np.random.default_rng(42)
    idx_case = rng.choice(len(x_case), size=n_sub, replace=False)
    idx_control = rng.choice(len(x_control), size=n_sub, replace=False)
    x_case_sub = x_case[idx_case]
    x_control_sub = x_control[idx_control]
    a = np.ones(n_sub) / n_sub
    b = np.ones(n_sub) / n_sub
    cost = ot.dist(x_case_sub, x_control_sub, metric="sqeuclidean")
    if np.max(cost) > 0:
        cost = cost / np.max(cost)
    plan = ot.sinkhorn(a, b, cost, reg=0.05, numItermax=1000, stopThr=1e-7)
    wasserstein = float(np.sum(plan * cost))

    pd.DataFrame({"wasserstein_distance": [wasserstein], "representation": [rep_key], "target_group": [target]}).to_csv(
        outdir / "cellot_summary.csv", index=False
    )
    outputs = [str(outdir / "cellot_summary.csv")]

    if "X_umap" in target_adata.obsm:
        umap = np.asarray(target_adata.obsm["X_umap"])
        case_umap = umap[case_mask][idx_case]
        control_umap = umap[control_mask][idx_control]
        transported = plan @ control_umap / (plan.sum(axis=1, keepdims=True) + 1e-8)
        displacement = transported - case_umap
        disp_df = pd.DataFrame(
            {
                "case_umap_1": case_umap[:, 0],
                "case_umap_2": case_umap[:, 1],
                "transport_umap_1": transported[:, 0],
                "transport_umap_2": transported[:, 1],
                "dx": displacement[:, 0],
                "dy": displacement[:, 1],
            }
        )
        disp_df.to_csv(outdir / "cellot_displacement.csv", index=False)
        outputs.append(str(outdir / "cellot_displacement.csv"))

        fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2))
        axes[0].scatter(umap[:, 0], umap[:, 1], c="#DDDDDD", s=4, alpha=0.25, linewidths=0)
        axes[0].scatter(case_umap[:, 0], case_umap[:, 1], c="#C41E3A", s=9, alpha=0.75, linewidths=0, label=context.config.case_condition)
        axes[0].scatter(control_umap[:, 0], control_umap[:, 1], c="#1B4F8A", s=9, alpha=0.75, linewidths=0, label=context.config.control_condition)
        axes[0].legend(frameon=False, fontsize=8)
        axes[0].set_title(f"{target}: case vs control")
        axes[0].set_xticks([])
        axes[0].set_yticks([])

        axes[1].scatter(umap[:, 0], umap[:, 1], c="#EEEEEE", s=3, alpha=0.2, linewidths=0)
        axes[1].scatter(case_umap[:, 0], case_umap[:, 1], c="#C41E3A", s=8, alpha=0.5, linewidths=0)
        step = max(1, n_sub // 60)
        for i in range(0, n_sub, step):
            axes[1].annotate("", xy=(transported[i, 0], transported[i, 1]), xytext=(case_umap[i, 0], case_umap[i, 1]), arrowprops=dict(arrowstyle="->", color="#1B4F8A", lw=0.6))
        axes[1].set_title(f"Optimal transport (W={wasserstein:.3f})")
        axes[1].set_xticks([])
        axes[1].set_yticks([])
        outputs.extend(save_figure_bundle(fig, outdir / "cellot_displacement"))

    return DomainResult(
        name=name,
        status="completed",
        message="Optimal-transport analysis completed.",
        outputs=outputs,
        metadata={"target_group": target, "representation": rep_key, "wasserstein_distance": wasserstein},
    )

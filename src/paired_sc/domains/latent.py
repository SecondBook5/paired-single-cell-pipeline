"""Latent representation domain built around scVI."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

from ..plotting.core import save_figure_bundle
from .base import DomainContext, DomainResult, dependency_available


def run(context: DomainContext) -> DomainResult:
    name = "latent"
    outdir = context.domain_dir(name)

    if not dependency_available("scvi"):
        return DomainResult(
            name=name,
            status="skipped",
            message="scvi-tools is not installed in the active environment.",
        )

    try:
        import scvi

        adata = context.adata.copy()
        if "counts" in context.adata.layers:
            adata.X = context.adata.layers["counts"].copy()
        elif context.adata.raw is not None:
            adata = context.adata.raw.to_adata()
            adata.obs = context.adata.obs.copy()
        if "highly_variable" in adata.var.columns and adata.var["highly_variable"].any():
            adata = adata[:, adata.var["highly_variable"]].copy()

        scvi.model.SCVI.setup_anndata(adata, batch_key=context.config.batch_key)
        model = scvi.model.SCVI(adata, n_layers=2, n_latent=context.config.domains.latent_n_latent)
        model.train(
            max_epochs=context.config.domains.latent_max_epochs,
            train_size=0.9,
            early_stopping=True,
            batch_size=256,
            check_val_every_n_epoch=1,
        )

        latent = model.get_latent_representation()
        context.adata.obsm["X_scvi"] = latent
        sc.pp.neighbors(
            context.adata,
            use_rep="X_scvi",
            n_neighbors=context.config.preprocess.neighbors_n_neighbors,
            key_added="neighbors_scvi",
        )
        sc.tl.umap(context.adata, neighbors_key="neighbors_scvi")
        context.adata.obsm["X_umap_scvi"] = context.adata.obsm["X_umap"].copy()

        outputs = []
        latent_path = outdir / "scvi_latent.csv"
        latent_df = pd.DataFrame(latent, index=context.adata.obs_names)
        latent_df.to_csv(latent_path)
        outputs.append(str(latent_path))

        fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.4))
        sc.pl.umap(
            context.adata,
            color=context.config.condition_key,
            palette=[context.config.condition_colors.get(cond, "#666666") for cond in context.config.condition_order],
            ax=axes[0],
            show=False,
            frameon=False,
        )
        axes[0].set_title("Condition in scVI latent UMAP")
        sc.pl.umap(
            context.adata,
            color=context.config.annotation.cell_type_key,
            ax=axes[1],
            show=False,
            frameon=False,
            legend_loc="on data",
            legend_fontsize=5,
        )
        axes[1].set_title("Annotation in scVI latent UMAP")
        outputs.extend(save_figure_bundle(fig, outdir / "scvi_umap"))

        model_dir = outdir / "scvi_model"
        model.save(model_dir, overwrite=True)
        outputs.append(str(model_dir))

        return DomainResult(
            name=name,
            status="completed",
            message="Latent scVI domain completed.",
            outputs=outputs,
            metadata={"n_latent": context.config.domains.latent_n_latent},
        )
    except Exception as exc:
        return DomainResult(name=name, status="failed", message=f"Latent domain failed ({exc}).")

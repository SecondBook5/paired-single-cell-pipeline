"""Generic plotting utilities and core figure builders."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns


def save_figure_bundle(fig, stem: Path) -> list[str]:
    outputs = []
    for suffix in [".png", ".pdf"]:
        path = stem.with_suffix(suffix)
        fig.savefig(path, dpi=300, bbox_inches="tight")
        outputs.append(str(path))
    plt.close(fig)
    return outputs


def plot_qc_summary(adata, config, figures_dir: Path) -> list[str]:
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    fig, axes = plt.subplots(1, len(metrics), figsize=(4.2 * len(metrics), 3.4))
    for ax, metric in zip(axes, metrics):
        sns.boxplot(
            data=adata.obs,
            x=config.condition_key,
            y=metric,
            order=config.condition_order,
            ax=ax,
        )
        ax.set_title(metric)
        ax.set_xlabel("")
    fig.suptitle("QC summaries by condition", y=1.02)
    return save_figure_bundle(fig, figures_dir / "qc_summary")


def plot_umap_condition(adata, config, figures_dir: Path) -> list[str]:
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    sc.pl.umap(
        adata,
        color=config.condition_key,
        palette=[config.condition_colors.get(cond, "#666666") for cond in config.condition_order],
        ax=ax,
        show=False,
        frameon=False,
    )
    ax.set_title("UMAP by condition")
    return save_figure_bundle(fig, figures_dir / "umap_condition")


def plot_umap_annotation(adata, config, figures_dir: Path) -> list[str]:
    fig, ax = plt.subplots(figsize=(5.8, 4.8))
    sc.pl.umap(
        adata,
        color=config.annotation.cell_type_key,
        ax=ax,
        show=False,
        frameon=False,
        legend_loc="on data",
        legend_fontsize=5,
    )
    ax.set_title("UMAP by annotation")
    return save_figure_bundle(fig, figures_dir / "umap_annotation")


def plot_composition(patient_level: pd.DataFrame, config, figures_dir: Path) -> list[str]:
    plot_df = patient_level.copy()
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    if plot_df.empty:
        ax.text(0.5, 0.5, "No paired composition rows available", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        return save_figure_bundle(fig, figures_dir / "cell_type_composition")
    sns.barplot(
        data=plot_df,
        x=config.annotation.cell_type_key,
        y="pct_cells",
        hue=config.condition_key,
        hue_order=config.condition_order,
        ax=ax,
    )
    ax.tick_params(axis="x", rotation=45)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("right")
    ax.set_ylabel("Percent of cells")
    ax.set_title("Patient-condition cell-type composition")
    fig.tight_layout()
    return save_figure_bundle(fig, figures_dir / "cell_type_composition")


def plot_differential_abundance(da_df: pd.DataFrame, figures_dir: Path) -> list[str]:
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    if da_df.empty:
        ax.text(0.5, 0.5, "No differential abundance rows available", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        return save_figure_bundle(fig, figures_dir / "differential_abundance")
    x_key = "cell_type" if "cell_type" in da_df.columns else da_df.columns[0]
    sns.barplot(
        data=da_df,
        x=x_key,
        y="mean_delta_pct_points",
        ax=ax,
        color="#4C78A8",
    )
    ax.axhline(0, color="#222222", linewidth=0.8)
    ax.tick_params(axis="x", rotation=45)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("right")
    ax.set_ylabel("Case - control (pct points)")
    ax.set_title("Differential abundance by cell type")
    fig.tight_layout()
    return save_figure_bundle(fig, figures_dir / "differential_abundance")

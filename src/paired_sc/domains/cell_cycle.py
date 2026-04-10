"""Cell-cycle scoring domain."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

from ..plotting.core import save_figure_bundle
from .base import DomainContext, DomainResult

S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
    "CDCA7", "DTL", "PRIM1", "UHRF1", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1",
    "GMNN", "WDR76", "SLBP", "CCNE2", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45",
    "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "USP1", "CLSPN", "POLA1", "CHAF1B",
    "BRIP1", "E2F8",
]
G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
    "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2",
    "CKAP2L", "AURKB", "BUB1", "KIF11", "GTSE1", "KIF20B", "HJURP", "CDCA3", "CDC20",
    "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2",
    "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "CKAP5", "CENPE", "NEK2", "GAS2L3", "CENPA",
]


def run(context: DomainContext) -> DomainResult:
    name = "cell_cycle"
    outdir = context.domain_dir(name)
    s_genes = [gene for gene in S_GENES if gene in context.adata.var_names]
    g2m_genes = [gene for gene in G2M_GENES if gene in context.adata.var_names]
    if len(s_genes) < 5 or len(g2m_genes) < 5:
        return DomainResult(
            name=name,
            status="skipped",
            message="Too few canonical S/G2M genes were present for reliable cell-cycle scoring.",
        )

    sc.tl.score_genes_cell_cycle(context.adata, s_genes=s_genes, g2m_genes=g2m_genes)
    outputs: list[str] = []

    phase_summary = (
        context.adata.obs.groupby([context.config.condition_key, "phase"], observed=False)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    phase_summary.to_csv(outdir / "cell_cycle_phase_summary.csv", index=False)
    outputs.append(str(outdir / "cell_cycle_phase_summary.csv"))

    if context.config.annotation.cell_type_key in context.adata.obs.columns:
        phase_celltype = (
            context.adata.obs.groupby([context.config.annotation.cell_type_key, "phase"], observed=False)
            .size()
            .rename("n_cells")
            .reset_index()
        )
        phase_celltype.to_csv(outdir / "cell_cycle_phase_by_annotation.csv", index=False)
        outputs.append(str(outdir / "cell_cycle_phase_by_annotation.csv"))

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8))
    sns.countplot(
        data=context.adata.obs,
        x="phase",
        order=["G1", "S", "G2M"],
        color="#4C78A8",
        ax=axes[0],
    )
    axes[0].set_title("Cell-cycle phase distribution")
    axes[0].set_xlabel("Phase")
    axes[0].set_ylabel("Cells")

    sns.countplot(
        data=context.adata.obs,
        x=context.config.condition_key,
        hue="phase",
        order=context.config.condition_order,
        hue_order=["G1", "S", "G2M"],
        ax=axes[1],
    )
    axes[1].set_title("Cell-cycle phase by condition")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("Cells")
    outputs.extend(save_figure_bundle(fig, outdir / "cell_cycle_summary"))

    if "X_umap" in context.adata.obsm:
        fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.8))
        sc.pl.umap(context.adata, color="phase", ax=axes[0], show=False, frameon=False)
        axes[0].set_title("Phase")
        sc.pl.umap(context.adata, color="S_score", cmap="viridis", ax=axes[1], show=False, frameon=False)
        axes[1].set_title("S score")
        sc.pl.umap(context.adata, color="G2M_score", cmap="viridis", ax=axes[2], show=False, frameon=False)
        axes[2].set_title("G2/M score")
        outputs.extend(save_figure_bundle(fig, outdir / "cell_cycle_umap"))

    return DomainResult(
        name=name,
        status="completed",
        message="Cell-cycle scoring completed.",
        outputs=outputs,
        metadata={"n_s_genes": len(s_genes), "n_g2m_genes": len(g2m_genes)},
    )

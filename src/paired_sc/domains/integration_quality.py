"""Integration quality metric domain."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ..plotting.core import save_figure_bundle
from .base import DomainContext, DomainResult, dependency_available


def _stratified_subsample(labels: pd.Series, n_total: int, seed: int = 42) -> np.ndarray:
    rng = np.random.default_rng(seed)
    labels = labels.astype(str)
    counts = labels.value_counts()
    frac = min(1.0, n_total / len(labels))
    chosen = []
    for label, count in counts.items():
        idx = np.flatnonzero(labels.to_numpy() == label)
        n_pick = min(len(idx), max(10, int(round(count * frac))))
        chosen.append(rng.choice(idx, size=n_pick, replace=False))
    merged = np.unique(np.concatenate(chosen))
    if len(merged) > n_total:
        merged = rng.choice(merged, size=n_total, replace=False)
    return np.sort(merged)


def _compute_lisi(embedding: np.ndarray, labels: np.ndarray, n_neighbors: int) -> np.ndarray:
    from sklearn.neighbors import NearestNeighbors

    k = max(1, min(int(n_neighbors), int(len(labels) - 1)))
    nbrs = NearestNeighbors(n_neighbors=k, metric="euclidean")
    nbrs.fit(embedding)
    idx = nbrs.kneighbors(return_distance=False)[:, 1:]
    out = np.empty(len(labels), dtype=float)
    for i, neighbors in enumerate(idx):
        _, counts = np.unique(labels[neighbors], return_counts=True)
        probs = counts / counts.sum()
        out[i] = 1.0 / np.sum(probs ** 2)
    return out


def run(context: DomainContext) -> DomainResult:
    name = "integration_quality"
    outdir = context.domain_dir(name)
    if not dependency_available("sklearn"):
        return DomainResult(name=name, status="skipped", message="scikit-learn is not installed in the active environment.")

    required = {
        "PCA": np.asarray(context.adata.obsm["X_pca"])[:, :30],
        "Harmony": np.asarray(context.adata.obsm["X_pca_harmony"])[:, :30],
    }
    if "X_scvi" in context.adata.obsm:
        required["scVI"] = np.asarray(context.adata.obsm["X_scvi"])[:, :30]

    if context.config.annotation.cell_type_key not in context.adata.obs.columns:
        return DomainResult(name=name, status="skipped", message="Cell-type annotations are required for integration-quality metrics.")

    sub_idx = _stratified_subsample(
        context.adata.obs[context.config.annotation.cell_type_key],
        n_total=min(context.config.domains.integration_subsample_n, context.adata.n_obs),
    )
    obs_sub = context.adata.obs.iloc[sub_idx].copy()
    outputs: list[str] = []
    rows: list[dict] = []
    for embedding_name, embedding in required.items():
        emb_sub = embedding[sub_idx]
        batch_lisi = _compute_lisi(emb_sub, obs_sub[context.config.batch_key].astype(str).to_numpy(), context.config.domains.integration_n_neighbors)
        sample_lisi = _compute_lisi(emb_sub, obs_sub[context.config.sample_key].astype(str).to_numpy(), context.config.domains.integration_n_neighbors)
        celltype_lisi = _compute_lisi(emb_sub, obs_sub[context.config.annotation.cell_type_key].astype(str).to_numpy(), context.config.domains.integration_n_neighbors)
        rows.extend({"embedding": embedding_name, "metric": "batch_iLISI", "value": float(v)} for v in batch_lisi)
        rows.extend({"embedding": embedding_name, "metric": "sample_iLISI", "value": float(v)} for v in sample_lisi)
        rows.extend({"embedding": embedding_name, "metric": "celltype_cLISI", "value": float(v)} for v in celltype_lisi)

    metrics = pd.DataFrame(rows)
    summary = (
        metrics.groupby(["embedding", "metric"], observed=False)["value"]
        .agg(median="median", mean="mean", q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75))
        .reset_index()
    )
    metrics.to_csv(outdir / "integration_quality_metrics_long.csv", index=False)
    summary.to_csv(outdir / "integration_quality_summary.csv", index=False)
    outputs.extend(
        [
            str(outdir / "integration_quality_metrics_long.csv"),
            str(outdir / "integration_quality_summary.csv"),
        ]
    )

    fig, axes = plt.subplots(1, 3, figsize=(9.0, 3.4))
    metric_order = ["batch_iLISI", "sample_iLISI", "celltype_cLISI"]
    for ax, metric in zip(axes, metric_order):
        sns.boxplot(data=metrics[metrics["metric"].eq(metric)], x="embedding", y="value", ax=ax, width=0.55, fliersize=0)
        ax.set_title(metric)
        ax.set_xlabel("")
        ax.set_ylabel("LISI")
    outputs.extend(save_figure_bundle(fig, outdir / "integration_quality"))

    return DomainResult(
        name=name,
        status="completed",
        message="Integration-quality metrics completed.",
        outputs=outputs,
        metadata={"n_embeddings": int(len(required))},
    )

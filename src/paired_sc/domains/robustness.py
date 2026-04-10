"""Donor leave-one-out robustness domain."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

from ..plotting.core import save_figure_bundle
from ..stats.paired import compute_cell_type_composition
from .base import DomainContext, DomainResult


def _paired_wilcoxon(values: np.ndarray) -> float:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) < 2 or np.allclose(vals, 0):
        return np.nan
    try:
        return float(stats.wilcoxon(vals, alternative="two-sided", zero_method="wilcox").pvalue)
    except ValueError:
        return np.nan


def run(context: DomainContext) -> DomainResult:
    name = "robustness"
    outdir = context.domain_dir(name)
    _, patient_level = compute_cell_type_composition(context.adata, context.config)
    cell_type_key = context.config.annotation.cell_type_key
    outputs: list[str] = []
    iteration_rows: list[dict] = []
    summary_rows: list[dict] = []

    for cell_type, sub in patient_level.groupby(cell_type_key, observed=False):
        wide = (
            sub.pivot_table(
                index=context.config.donor_key,
                columns=context.config.condition_key,
                values="pct_cells",
                aggfunc="first",
                observed=False,
            )
            .dropna()
            .reset_index()
        )
        if len(wide) < 3:
            continue
        delta = wide[context.config.case_condition] - wide[context.config.control_condition]
        full_delta = float(delta.mean())
        full_p = _paired_wilcoxon(delta.to_numpy())
        loo_values = []
        for _, row in wide.iterrows():
            donor_id = row[context.config.donor_key]
            sub_wide = wide[wide[context.config.donor_key] != donor_id].copy()
            if len(sub_wide) < 2:
                continue
            sub_delta = sub_wide[context.config.case_condition] - sub_wide[context.config.control_condition]
            loo_delta = float(sub_delta.mean())
            loo_p = _paired_wilcoxon(sub_delta.to_numpy())
            loo_values.append(loo_delta)
            iteration_rows.append(
                {
                    cell_type_key: cell_type,
                    "omitted_donor": donor_id,
                    "loo_delta_pct_points": loo_delta,
                    "loo_p": loo_p,
                    "full_delta_pct_points": full_delta,
                    "full_p": full_p,
                }
            )
        if loo_values:
            summary_rows.append(
                {
                    cell_type_key: cell_type,
                    "n_pairs": int(len(wide)),
                    "full_delta_pct_points": full_delta,
                    "full_p": full_p,
                    "loo_min": float(np.min(loo_values)),
                    "loo_max": float(np.max(loo_values)),
                    "same_direction_fraction": float(np.mean(np.sign(loo_values) == np.sign(full_delta))),
                }
            )

    if not summary_rows:
        return DomainResult(name=name, status="skipped", message="Too few paired donors were available for leave-one-out robustness checks.")

    iterations = pd.DataFrame(iteration_rows)
    summary = pd.DataFrame(summary_rows).sort_values("same_direction_fraction", ascending=False)
    iterations.to_csv(outdir / "robustness_leave_one_out_iterations.csv", index=False)
    summary.to_csv(outdir / "robustness_summary.csv", index=False)
    outputs.extend(
        [
            str(outdir / "robustness_leave_one_out_iterations.csv"),
            str(outdir / "robustness_summary.csv"),
        ]
    )

    top_n = min(context.config.domains.robustness_top_n, len(summary))
    plot_df = summary.head(top_n).copy()
    fig, ax = plt.subplots(figsize=(7.0, max(3.6, 0.45 * top_n)))
    sns.barplot(data=plot_df, x="full_delta_pct_points", y=cell_type_key, color="#7A5195", ax=ax)
    for idx, row in plot_df.reset_index(drop=True).iterrows():
        ax.plot([row["loo_min"], row["loo_max"]], [idx, idx], color="#333333", lw=1.1)
    ax.axvline(0, color="#222222", lw=0.8)
    ax.set_title("Donor leave-one-out robustness")
    ax.set_xlabel("Case - control (pct points)")
    ax.set_ylabel("")
    outputs.extend(save_figure_bundle(fig, outdir / "robustness_summary"))

    return DomainResult(
        name=name,
        status="completed",
        message="Donor leave-one-out robustness analysis completed.",
        outputs=outputs,
        metadata={"n_annotations": int(len(summary))},
    )

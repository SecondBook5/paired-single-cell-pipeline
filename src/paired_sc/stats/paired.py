"""Donor-aware summaries for paired tissue designs."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse, stats
import statsmodels.formula.api as smf


def _extract_counts_sum(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.sum(axis=1)).ravel()
    return np.asarray(matrix.sum(axis=1)).ravel()


def compute_cell_type_composition(adata, config) -> tuple[pd.DataFrame, pd.DataFrame]:
    obs = adata.obs.copy()
    cell_type_key = config.annotation.cell_type_key
    total_by_sample = obs.groupby(config.sample_key, observed=False).size().rename("total_cells").reset_index()
    sample_counts = (
        obs.groupby(
            [config.sample_key, config.donor_key, config.condition_key, config.batch_key, cell_type_key],
            observed=False,
        )
        .size()
        .rename("n_cells")
        .reset_index()
    )
    sample_counts = sample_counts.merge(total_by_sample, on=config.sample_key, how="left")
    sample_counts["proportion"] = sample_counts["n_cells"] / sample_counts["total_cells"]
    sample_counts["pct_cells"] = sample_counts["proportion"] * 100

    patient_counts = (
        sample_counts.groupby([config.donor_key, config.condition_key, cell_type_key], as_index=False, observed=False)
        .agg(n_cells=("n_cells", "sum"), total_cells=("total_cells", "sum"))
    )
    patient_counts["proportion"] = patient_counts["n_cells"] / patient_counts["total_cells"]
    patient_counts["pct_cells"] = patient_counts["proportion"] * 100
    return sample_counts, patient_counts


def compute_pseudobulk_profiles(adata, config) -> pd.DataFrame:
    obs = adata.obs.copy()
    gene_index = adata.raw.var_names if adata.raw is not None else adata.var_names
    matrix = adata.raw.X if adata.raw is not None else (adata.layers["counts"] if "counts" in adata.layers else adata.X)

    group_meta = obs[
        [config.donor_key, config.condition_key, config.annotation.cell_type_key]
    ].astype(str).copy()
    group_labels = (
        group_meta[config.donor_key]
        + "__"
        + group_meta[config.condition_key]
        + "__"
        + group_meta[config.annotation.cell_type_key]
    )
    codes, unique_labels = pd.factorize(group_labels, sort=False)
    n_groups = len(unique_labels)

    indicator = sparse.coo_matrix(
        (np.ones(len(codes), dtype=np.float32), (codes, np.arange(len(codes)))),
        shape=(n_groups, adata.n_obs),
    ).tocsr()
    summed = indicator @ matrix
    summed = summed.toarray() if sparse.issparse(summed) else np.asarray(summed)

    meta = (
        group_meta.assign(_group=group_labels)
        .drop_duplicates("_group")
        .set_index("_group")
        .reindex(unique_labels)
        .reset_index(drop=True)
    )
    total_counts = summed.sum(axis=1)

    wide = pd.DataFrame(summed, columns=gene_index)
    wide.insert(0, "pseudobulk_total_counts", total_counts)
    wide.insert(0, "cell_type", meta[config.annotation.cell_type_key].values)
    wide.insert(0, "condition", meta[config.condition_key].values)
    wide.insert(0, "donor_id", meta[config.donor_key].values)
    return wide


def _paired_wilcoxon(values: np.ndarray) -> float:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) < 2 or np.allclose(vals, 0):
        return np.nan
    try:
        return float(stats.wilcoxon(vals, alternative="two-sided", zero_method="wilcox").pvalue)
    except ValueError:
        return np.nan


def _fit_mixed_effect(df: pd.DataFrame, case_label: str) -> tuple[float, float, str]:
    sub = df.copy()
    sub["is_case"] = sub["condition"].eq(case_label).astype(int)
    sub["asin_sqrt_prop"] = np.arcsin(np.sqrt(np.clip(sub["proportion"], 1e-8, 1 - 1e-8)))
    try:
        fit = smf.mixedlm(
            "asin_sqrt_prop ~ is_case",
            data=sub,
            groups=sub["donor_id"].astype(str),
        ).fit(reml=False, method="lbfgs", disp=False)
        return float(fit.params["is_case"]), float(fit.pvalues["is_case"]), "mixedlm"
    except Exception:
        fit = smf.ols("asin_sqrt_prop ~ is_case + C(donor_id)", data=sub).fit()
        return float(fit.params["is_case"]), float(fit.pvalues["is_case"]), "blocked_ols"


def compute_differential_abundance(patient_level: pd.DataFrame, config) -> pd.DataFrame:
    cell_type_key = config.annotation.cell_type_key
    rows = []
    for cell_type, sub in patient_level.groupby(cell_type_key, observed=True):
        wide = (
            sub.pivot_table(
                index=config.donor_key,
                columns=config.condition_key,
                values="pct_cells",
                aggfunc="first",
                observed=False,
            )
            .dropna()
            .reset_index()
        )
        if {config.case_condition, config.control_condition}.issubset(wide.columns):
            deltas = wide[config.case_condition] - wide[config.control_condition]
            paired_p = _paired_wilcoxon(deltas.to_numpy())
            effect, model_p, model_name = _fit_mixed_effect(
                sub[[config.donor_key, config.condition_key, "proportion"]]
                .rename(columns={config.donor_key: "donor_id", config.condition_key: "condition"}),
                config.case_condition,
            )
            rows.append(
                {
                    cell_type_key: cell_type,
                    "n_donors": int(len(wide)),
                    "mean_control_pct": float(wide[config.control_condition].mean()),
                    "mean_case_pct": float(wide[config.case_condition].mean()),
                    "mean_delta_pct_points": float(deltas.mean()),
                    "paired_wilcoxon_p": paired_p,
                    "model_effect": effect,
                    "model_p": model_p,
                    "model_name": model_name,
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values("paired_wilcoxon_p", na_position="last")

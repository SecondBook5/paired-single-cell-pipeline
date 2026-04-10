"""Manifest validation and 10x H5 loading."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import scanpy as sc
from pydantic import BaseModel, ConfigDict, Field


REQUIRED_COLUMNS = ["sample_id", "donor_id", "condition", "batch", "matrix_h5"]
OPTIONAL_COLUMNS = ["sample_name", "replicate_group", "metrics_csv", "metadata_json"]


class ManifestTable(BaseModel):
    """Validated manifest wrapper."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    path: Path
    dataframe: pd.DataFrame = Field(repr=False)

    @classmethod
    def from_csv(cls, path: str | Path) -> "ManifestTable":
        csv_path = Path(path).expanduser().resolve()
        df = pd.read_csv(csv_path).copy()
        missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
        if missing:
            raise ValueError(f"Manifest missing required columns: {missing}")
        if df["sample_id"].duplicated().any():
            dupes = df.loc[df["sample_id"].duplicated(), "sample_id"].tolist()
            raise ValueError(f"Manifest has duplicate sample_id values: {dupes}")
        for col in REQUIRED_COLUMNS:
            if df[col].isna().any():
                raise ValueError(f"Manifest column {col!r} contains null values")
        for col in OPTIONAL_COLUMNS:
            if col not in df.columns:
                df[col] = pd.NA
        return cls(path=csv_path, dataframe=df)

    def validate_conditions(self, config) -> None:
        observed = set(self.dataframe[config.condition_key].astype(str))
        expected = {config.control_condition, config.case_condition}
        unexpected = observed - expected
        missing = expected - observed
        if unexpected:
            raise ValueError(f"Unexpected condition labels in manifest: {sorted(unexpected)}")
        if missing:
            raise ValueError(f"Manifest is missing expected condition labels: {sorted(missing)}")

    def to_records(self) -> list[dict]:
        return self.dataframe.to_dict(orient="records")


def load_manifest_adata(manifest: ManifestTable, config) -> sc.AnnData:
    """Load and concatenate all 10x H5 inputs declared in the manifest."""
    adatas = []
    for row in manifest.to_records():
        matrix_path = Path(str(row["matrix_h5"])).expanduser()
        if not matrix_path.is_absolute():
            matrix_path = (manifest.path.parent / matrix_path).resolve()
        adata = sc.read_10x_h5(str(matrix_path))
        adata.var_names_make_unique()
        adata.obs[config.sample_key] = str(row["sample_id"])
        adata.obs[config.donor_key] = str(row["donor_id"])
        adata.obs[config.condition_key] = str(row["condition"])
        adata.obs[config.batch_key] = str(row["batch"])
        if row.get("sample_name") is not None and not pd.isna(row.get("sample_name")):
            adata.obs["sample_name"] = str(row["sample_name"])
        if config.replicate_key:
            replicate_value = row.get(config.replicate_key)
            if replicate_value is not None and not pd.isna(replicate_value):
                adata.obs[config.replicate_key] = str(replicate_value)
        adatas.append(adata)

    combined = sc.concat(
        adatas,
        join="outer",
        label=config.sample_key,
        keys=[ad.obs[config.sample_key].iloc[0] for ad in adatas],
        index_unique="-",
    )
    combined.var_names_make_unique()
    return combined

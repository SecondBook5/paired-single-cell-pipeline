"""Pydantic workflow configuration models."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator


class OutputLayout(BaseModel):
    results_dir: str = "results"
    figures_dir: str = "figures"
    reports_dir: str = "reports"
    logs_dir: str = "logs"


class QCConfig(BaseModel):
    min_genes: int = 200
    min_counts: int = 500
    max_counts: int | None = None
    max_pct_mt: float = 20.0
    min_cells_per_gene: int = 10


class PreprocessConfig(BaseModel):
    target_sum: float = 10000.0
    hvg_flavor: str = "seurat_v3"
    n_top_genes: int = 3000
    scale_max_value: float = 10.0
    pca_n_comps: int = 50
    harmony_max_iter: int = 20
    neighbors_n_neighbors: int = 15
    neighbors_n_pcs: int = 30
    umap_min_dist: float = 0.3
    umap_spread: float = 1.0


class AnnotationConfig(BaseModel):
    backend: str = "celltypist"
    model: str = "Adult_Human_Skin.pkl"
    primary_leiden_resolution: float = 0.5
    leiden_resolutions: list[float] = Field(default_factory=lambda: [0.3, 0.5, 0.8, 1.2])
    cell_type_key: str = "cell_type"
    label_map_path: str | None = None

    @property
    def primary_leiden_key(self) -> str:
        return f"leiden_res{self.primary_leiden_resolution}"


class DomainConfig(BaseModel):
    enabled: list[str] = Field(
        default_factory=lambda: [
            "cell_cycle",
            "differential_expression",
            "pathway_enrichment",
            "liana",
            "magic",
            "trajectory",
            "latent",
            "regulatory",
            "cellot",
            "pseudobulk_de",
            "pathway_activity",
            "robustness",
            "integration_quality",
        ]
    )
    de_min_cells_per_condition: int = 10
    de_logfc_threshold: float = 0.5
    pathway_gene_sets: list[str] = Field(
        default_factory=lambda: [
            "GO_Biological_Process_2023",
            "KEGG_2021_Human",
            "Reactome_2022",
        ]
    )
    pathway_top_terms: int = 15
    pathway_top_genes: int = 100
    liana_expr_prop: float = 0.1
    liana_min_cells: int = 10
    magic_genes: list[str] = Field(default_factory=list)
    magic_knn: int = 5
    magic_solver: str = "approximate"
    trajectory_root_group: str | None = None
    trajectory_groupby: str | None = None
    trajectory_enable_palantir: bool = True
    trajectory_enable_cellrank: bool = True
    trajectory_max_cells: int = 5000
    latent_n_latent: int = 30
    latent_max_epochs: int = 50
    cellot_subsample_n: int = 500
    pseudobulk_min_cells_per_condition: int = 20
    pseudobulk_min_replicates_per_condition: int = 2
    pseudobulk_min_count: int = 10
    pathway_activity_library: str = "MSigDB_Hallmark_2020"
    pathway_activity_min_genes: int = 5
    robustness_top_n: int = 8
    integration_subsample_n: int = 12000
    integration_n_neighbors: int = 30
    target_groupby: str | None = None
    target_group: str | None = None
    target_marker_genes: list[str] = Field(default_factory=list)
    target_score_genes: list[str] = Field(default_factory=list)
    target_subcluster_resolution: float = 0.4
    target_subcluster_hvgs: int = 2000
    target_subcluster_neighbors: int = 20
    regulatory_tf_list: str | None = None
    regulatory_rankings: str | None = None
    regulatory_motifs: str | None = None
    regulatory_top_variable_tfs: int = 25
    regulatory_max_cells: int = 5000


class WorkflowConfig(BaseModel):
    """Validated project-level workflow configuration."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    project_name: str
    project_subtitle: str = ""
    condition_key: str = "condition"
    case_condition: str
    control_condition: str
    donor_key: str = "donor_id"
    batch_key: str = "batch"
    sample_key: str = "sample_id"
    replicate_key: str | None = "replicate_group"
    condition_colors: dict[str, str] = Field(default_factory=dict)
    qc: QCConfig = Field(default_factory=QCConfig)
    preprocess: PreprocessConfig = Field(default_factory=PreprocessConfig)
    annotation: AnnotationConfig = Field(default_factory=AnnotationConfig)
    domains: DomainConfig = Field(default_factory=DomainConfig)
    outputs: OutputLayout = Field(default_factory=OutputLayout)
    config_source: str | None = None

    @model_validator(mode="after")
    def validate_conditions(self) -> "WorkflowConfig":
        if self.case_condition == self.control_condition:
            raise ValueError("case_condition and control_condition must differ")
        if self.annotation.primary_leiden_resolution not in self.annotation.leiden_resolutions:
            raise ValueError("primary_leiden_resolution must be present in leiden_resolutions")
        return self

    @property
    def condition_order(self) -> list[str]:
        return [self.control_condition, self.case_condition]

    @classmethod
    def from_yaml(cls, path: str | Path) -> "WorkflowConfig":
        cfg_path = Path(path).expanduser().resolve()
        with cfg_path.open("r", encoding="utf-8") as handle:
            payload = yaml.safe_load(handle) or {}
        payload["config_source"] = str(cfg_path)
        return cls.model_validate(payload)

    def resolve_optional_path(self, value: str | None) -> Path | None:
        if not value:
            return None
        base = Path(self.config_source).parent if self.config_source else Path.cwd()
        return (base / value).resolve() if not Path(value).is_absolute() else Path(value).resolve()

    def to_serializable(self) -> dict[str, Any]:
        data = self.model_dump()
        data["condition_order"] = self.condition_order
        data["annotation"]["primary_leiden_key"] = self.annotation.primary_leiden_key
        return data

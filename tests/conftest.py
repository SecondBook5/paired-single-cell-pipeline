from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest
import yaml
from scipy import sparse


GENES = ["MT-CO1", "KRT14", "EPCAM", "COL1A1", "PTPRC", "MKI67"]


def write_10x_h5(path: Path, matrix: np.ndarray, gene_names: list[str], barcodes: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    csc = sparse.csc_matrix(matrix.astype(np.int32))
    with h5py.File(path, "w") as handle:
        grp = handle.create_group("matrix")
        grp.create_dataset("data", data=csc.data.astype(np.int32))
        grp.create_dataset("indices", data=csc.indices.astype(np.int64))
        grp.create_dataset("indptr", data=csc.indptr.astype(np.int64))
        grp.create_dataset("shape", data=np.asarray(csc.shape, dtype=np.int64))
        grp.create_dataset("barcodes", data=np.asarray([value.encode("utf-8") for value in barcodes]))

        features = grp.create_group("features")
        encoded_genes = np.asarray([value.encode("utf-8") for value in gene_names])
        features.create_dataset("id", data=encoded_genes)
        features.create_dataset("name", data=encoded_genes)
        features.create_dataset("feature_type", data=np.asarray([b"Gene Expression"] * len(gene_names)))
        features.create_dataset("genome", data=np.asarray([b"GRCh38"] * len(gene_names)))


@pytest.fixture()
def smoke_inputs(tmp_path: Path) -> dict[str, Path]:
    sample_dir = tmp_path / "matrices"
    matrices = {
        "N1": np.array(
            [
                [1, 1, 2],
                [8, 9, 7],
                [7, 8, 8],
                [0, 1, 0],
                [0, 0, 0],
                [1, 0, 1],
            ]
        ),
        "L1": np.array(
            [
                [1, 2, 2],
                [7, 8, 6],
                [6, 7, 6],
                [1, 0, 1],
                [0, 1, 0],
                [1, 1, 0],
            ]
        ),
        "N2": np.array(
            [
                [2, 1, 1],
                [9, 8, 8],
                [8, 7, 7],
                [0, 0, 1],
                [0, 0, 0],
                [0, 1, 0],
            ]
        ),
        "L2": np.array(
            [
                [1, 2, 1],
                [6, 7, 6],
                [6, 6, 5],
                [1, 1, 1],
                [0, 0, 1],
                [1, 0, 1],
            ]
        ),
    }

    manifest_rows = []
    for sample_id, matrix in matrices.items():
        donor_id = "P1" if sample_id.endswith("1") else "P2"
        condition = "Normal" if sample_id.startswith("N") else "LE"
        batch = "SetA" if donor_id == "P1" else "SetB"
        h5_path = sample_dir / f"{sample_id}.h5"
        barcodes = [f"{sample_id}_cell{i}" for i in range(matrix.shape[1])]
        write_10x_h5(h5_path, matrix, GENES, barcodes)
        manifest_rows.append(
            {
                "sample_id": sample_id,
                "donor_id": donor_id,
                "condition": condition,
                "batch": batch,
                "matrix_h5": str(h5_path),
            }
        )

    manifest_path = tmp_path / "manifest.csv"
    import pandas as pd

    pd.DataFrame(manifest_rows).to_csv(manifest_path, index=False)

    tf_list_path = tmp_path / "tf_list.txt"
    tf_list_path.write_text("KRT14\nEPCAM\nMKI67\n", encoding="utf-8")

    project_payload = {
        "project_name": "smoke_project",
        "project_subtitle": "Tiny paired tissue smoke test",
        "condition_key": "condition",
        "case_condition": "LE",
        "control_condition": "Normal",
        "donor_key": "donor_id",
        "batch_key": "batch",
        "sample_key": "sample_id",
        "replicate_key": None,
        "condition_colors": {"Normal": "#1B4F8A", "LE": "#C41E3A"},
        "qc": {
            "min_genes": 1,
            "min_counts": 1,
            "max_counts": None,
            "max_pct_mt": 100.0,
            "min_cells_per_gene": 1,
        },
        "preprocess": {
            "target_sum": 1000.0,
            "hvg_flavor": "cell_ranger",
            "n_top_genes": 5,
            "scale_max_value": 10.0,
            "pca_n_comps": 2,
            "harmony_max_iter": 5,
            "neighbors_n_neighbors": 2,
            "neighbors_n_pcs": 2,
            "umap_min_dist": 0.3,
            "umap_spread": 1.0,
        },
        "annotation": {
            "backend": "none",
            "model": "Adult_Human_Skin.pkl",
            "primary_leiden_resolution": 0.4,
            "leiden_resolutions": [0.2, 0.4],
            "cell_type_key": "cell_type",
            "label_map_path": None,
        },
        "domains": {
            "enabled": ["trajectory", "regulatory"],
            "trajectory_groupby": "cell_type",
            "magic_genes": ["KRT14", "EPCAM"],
            "regulatory_tf_list": str(tf_list_path),
        },
        "outputs": {
            "results_dir": "results",
            "figures_dir": "figures",
            "reports_dir": "reports",
            "logs_dir": "logs",
        },
    }
    project_path = tmp_path / "project.yaml"
    project_path.write_text(yaml.safe_dump(project_payload, sort_keys=False), encoding="utf-8")

    return {
        "project": project_path,
        "manifest": manifest_path,
        "workdir": tmp_path / "run",
    }

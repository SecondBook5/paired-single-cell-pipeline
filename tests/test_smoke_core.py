from __future__ import annotations

from pathlib import Path

from paired_sc import run_advanced_domains, run_core_pipeline


def test_core_workflow_smoke(smoke_inputs) -> None:
    outputs = run_core_pipeline(smoke_inputs["project"], smoke_inputs["manifest"], smoke_inputs["workdir"])
    expected = [
        "adata_raw_with_qc",
        "adata_filtered",
        "adata_preprocessed",
        "adata_annotated",
        "cell_type_composition_sample_level",
        "cell_type_composition_pair_level",
        "differential_abundance",
        "pseudobulk_profiles",
        "annotation_summary",
        "core_run_summary",
    ]
    for key in expected:
        assert key in outputs
        assert Path(outputs[key]).exists()


def test_advanced_domains_smoke(smoke_inputs) -> None:
    run_core_pipeline(smoke_inputs["project"], smoke_inputs["manifest"], smoke_inputs["workdir"])
    results = run_advanced_domains(
        smoke_inputs["project"],
        smoke_inputs["manifest"],
        smoke_inputs["workdir"],
        domains=["trajectory", "regulatory"],
    )
    status_map = {row["name"]: row["status"] for row in results}
    assert status_map["trajectory"] == "completed"
    assert status_map["regulatory"] == "completed"


def test_post_core_domains_smoke(smoke_inputs) -> None:
    run_core_pipeline(smoke_inputs["project"], smoke_inputs["manifest"], smoke_inputs["workdir"])
    results = run_advanced_domains(
        smoke_inputs["project"],
        smoke_inputs["manifest"],
        smoke_inputs["workdir"],
        domains=["cell_cycle", "differential_expression", "pseudobulk_de", "robustness", "integration_quality"],
    )
    status_map = {row["name"]: row["status"] for row in results}
    assert status_map["cell_cycle"] in {"completed", "skipped"}
    assert status_map["differential_expression"] in {"completed", "skipped"}
    assert status_map["pseudobulk_de"] in {"completed", "skipped"}
    assert status_map["robustness"] in {"completed", "skipped"}
    assert status_map["integration_quality"] == "completed"


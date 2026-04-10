from __future__ import annotations

from paired_sc.domains import get_domain_status


def test_domain_status_surface() -> None:
    rows = get_domain_status()
    names = {row["name"] for row in rows}
    assert names == {
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
        "target_population",
        "target_subclustering",
    }
    for row in rows:
        assert "available" in row
        assert "missing_dependencies" in row


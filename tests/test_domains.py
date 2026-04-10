from __future__ import annotations

from paired_sc.domains import get_domain_status


def test_domain_status_surface() -> None:
    rows = get_domain_status()
    names = {row["name"] for row in rows}
    assert names == {"liana", "magic", "trajectory", "latent", "regulatory"}
    for row in rows:
        assert "available" in row
        assert "missing_dependencies" in row


"""Annotation label adapter loading."""

from __future__ import annotations

import importlib.util
from pathlib import Path


def load_label_adapter(path: Path | None):
    if path is None:
        return None
    spec = importlib.util.spec_from_file_location("pts_label_adapter", str(path))
    if spec is None or spec.loader is None:
        raise ValueError(f"Could not load label adapter from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "remap_labels"):
        raise ValueError(f"Adapter {path} does not define remap_labels(labels)")
    return module.remap_labels


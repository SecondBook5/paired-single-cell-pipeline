"""Shared domain models and helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import importlib.util
from pydantic import BaseModel, ConfigDict, Field


class DomainResult(BaseModel):
    """Structured status for a callable analysis domain."""

    name: str
    status: str
    message: str
    outputs: list[str] = Field(default_factory=list)
    metadata: dict[str, Any] = Field(default_factory=dict)


class DomainContext(BaseModel):
    """State passed into each domain."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    config: Any
    manifest: Any
    paths: Any
    adata: Any

    def domain_dir(self, name: str) -> Path:
        out = self.paths.domains / name
        out.mkdir(parents=True, exist_ok=True)
        return out


def dependency_available(module_name: str) -> bool:
    return importlib.util.find_spec(module_name) is not None

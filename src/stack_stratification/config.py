"""Shared configuration for the scaffolded pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PipelinePaths:
    """Common path bundle used by CLI commands."""

    metadata: Path | None = None
    expression: Path | None = None
    manifest: Path | None = None
    embeddings: Path | None = None
    clusters: Path | None = None
    output: Path | None = None

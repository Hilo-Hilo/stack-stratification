"""Stack embedding scaffolding."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class EmbeddingJob:
    manifest_path: Path
    backend: str = "stack"


def build_embedding_job(job: EmbeddingJob) -> dict[str, object]:
    return {
        "stage": "embed",
        "status": "placeholder",
        "backend": job.backend,
        "manifest_path": str(job.manifest_path),
        "output_format": "parquet",
        "notes": [
            "Resolve Stack API or local model access.",
            "Record model identifier and embedding version.",
            "Store per-sample embedding provenance.",
        ],
    }

"""Stratification scaffolding."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class StratificationJob:
    embeddings_path: Path
    method: str = "hdbscan_or_kmeans_baseline"


def build_stratification_job(job: StratificationJob) -> dict[str, object]:
    return {
        "stage": "stratify",
        "status": "placeholder",
        "embeddings_path": str(job.embeddings_path),
        "method": job.method,
        "outputs": ["cluster_assignments", "diagnostics", "run_metadata"],
        "notes": [
            "Select the first reproducible clustering baseline.",
            "Add stability diagnostics and embedding-space QC.",
            "Benchmark against simple transcriptomics baselines.",
        ],
    }

"""Cluster annotation scaffolding."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class AnnotationJob:
    clusters_path: Path
    metadata_path: Path


def build_annotation_job(job: AnnotationJob) -> dict[str, object]:
    return {
        "stage": "annotate",
        "status": "placeholder",
        "clusters_path": str(job.clusters_path),
        "metadata_path": str(job.metadata_path),
        "outputs": ["cluster_summary", "metadata_enrichment", "interpretation_notes"],
        "notes": [
            "Summarize clinical and genomic correlates by cluster.",
            "Keep outputs descriptive rather than prescriptive.",
            "Flag any trial-enrichment-adjacent language for careful review.",
        ],
    }

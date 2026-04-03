"""Cohort ingestion scaffolding."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class IngestInputs:
    cohort: str
    metadata_path: Path
    expression_path: Path


def build_manifest(inputs: IngestInputs) -> dict[str, object]:
    return {
        "stage": "ingest",
        "status": "placeholder",
        "cohort": inputs.cohort,
        "representation": "transcriptomics",
        "metadata_path": str(inputs.metadata_path),
        "expression_path": str(inputs.expression_path),
        "notes": [
            "Validate LUAD sample identifiers.",
            "Capture input schema and provenance.",
            "Add cohort-level QC summary in a future iteration.",
        ],
    }

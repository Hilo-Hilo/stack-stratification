"""Command-line interface for the Stack-Stratification scaffold."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from stack_stratification.annotate import AnnotationJob, build_annotation_job
from stack_stratification.embeddings import EmbeddingJob, build_embedding_job
from stack_stratification.ingest import IngestInputs, build_manifest
from stack_stratification.io import write_json
from stack_stratification.stratify import StratificationJob, build_stratification_job


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="stack-stratification",
        description="LUAD stratification scaffold built around Stack embeddings.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    ingest = subparsers.add_parser("ingest", help="Create a cohort manifest.")
    ingest.add_argument("--cohort", default="luad")
    ingest.add_argument("--metadata", required=True, type=Path)
    ingest.add_argument("--expression", required=True, type=Path)
    ingest.add_argument("--output", required=True, type=Path)
    ingest.set_defaults(func=_run_ingest)

    embed = subparsers.add_parser("embed", help="Create a Stack embedding job spec.")
    embed.add_argument("--manifest", required=True, type=Path)
    embed.add_argument("--backend", default="stack")
    embed.add_argument("--output", required=True, type=Path)
    embed.set_defaults(func=_run_embed)

    stratify = subparsers.add_parser("stratify", help="Create a stratification job spec.")
    stratify.add_argument("--embeddings", required=True, type=Path)
    stratify.add_argument("--method", default="hdbscan_or_kmeans_baseline")
    stratify.add_argument("--output", required=True, type=Path)
    stratify.set_defaults(func=_run_stratify)

    annotate = subparsers.add_parser("annotate", help="Create a cluster annotation job spec.")
    annotate.add_argument("--clusters", required=True, type=Path)
    annotate.add_argument("--metadata", required=True, type=Path)
    annotate.add_argument("--output", required=True, type=Path)
    annotate.set_defaults(func=_run_annotate)

    return parser


def _run_ingest(args: argparse.Namespace) -> int:
    payload = build_manifest(
        IngestInputs(
            cohort=args.cohort,
            metadata_path=args.metadata,
            expression_path=args.expression,
        )
    )
    write_json(args.output, payload)
    return 0


def _run_embed(args: argparse.Namespace) -> int:
    payload = build_embedding_job(EmbeddingJob(manifest_path=args.manifest, backend=args.backend))
    write_json(args.output, payload)
    return 0


def _run_stratify(args: argparse.Namespace) -> int:
    payload = build_stratification_job(
        StratificationJob(embeddings_path=args.embeddings, method=args.method)
    )
    write_json(args.output, payload)
    return 0


def _run_annotate(args: argparse.Namespace) -> int:
    payload = build_annotation_job(
        AnnotationJob(clusters_path=args.clusters, metadata_path=args.metadata)
    )
    write_json(args.output, payload)
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())

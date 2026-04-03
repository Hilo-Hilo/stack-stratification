#!/usr/bin/env python3
"""Compare radiological-type signal across coarse cell compartments."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.neighbors import NearestNeighbors

matplotlib.use("Agg")
os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")

REPO_ROOT = Path(__file__).resolve().parent.parent
COMPARTMENT_ORDER = ["all", "epithelial", "immune", "stromal", "endothelial"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--adata",
        type=Path,
        default=REPO_ROOT / "runs/luad_stage_benchmark_300/gse189357_stage_subset.h5ad",
    )
    parser.add_argument(
        "--stack-embeddings",
        type=Path,
        default=REPO_ROOT / "runs/luad_stage_benchmark_300/gse189357_stage_subset_stack_embeddings.h5ad",
    )
    parser.add_argument(
        "--compartments",
        type=Path,
        default=REPO_ROOT / "results/luad_stage_compartments/cell_compartments.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "results/luad_radiology_compartments",
    )
    parser.add_argument("--seed", type=int, default=7)
    return parser.parse_args()


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def build_pca_embedding(adata: ad.AnnData) -> np.ndarray:
    work = adata.copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    sc.pp.highly_variable_genes(work, n_top_genes=2000, flavor="seurat")
    work = work[:, work.var["highly_variable"].values].copy()
    sc.pp.pca(work, n_comps=50)
    return np.asarray(work.obsm["X_pca"])


def nearest_neighbor_accuracy(features: np.ndarray, labels: np.ndarray) -> float:
    nn = NearestNeighbors(n_neighbors=2, metric="euclidean")
    nn.fit(features)
    indices = nn.kneighbors(return_distance=False)
    return float(np.mean(labels[indices[:, 1]] == labels))


def evaluate(
    obs: pd.DataFrame,
    embeddings: np.ndarray,
    compartments: pd.Series,
    representation: str,
    seed: int,
) -> pd.DataFrame:
    label_map = {label: idx for idx, label in enumerate(sorted(obs["radiological_type"].unique()))}
    rows: list[dict[str, object]] = []

    for compartment in COMPARTMENT_ORDER:
        mask = np.ones(obs.shape[0], dtype=bool) if compartment == "all" else (
            compartments.values == compartment
        )
        subset_obs = obs.loc[mask].copy()
        if subset_obs.empty:
            continue
        subset_embeddings = embeddings[mask]

        centroid_df = (
            pd.DataFrame(subset_embeddings)
            .assign(
                sample_id=subset_obs["sample_id"].values,
                radiological_type=subset_obs["radiological_type"].values,
            )
            .groupby(["sample_id", "radiological_type"], as_index=False, observed=True)
            .mean()
        )
        if centroid_df.shape[0] < 4:
            continue

        centroid_matrix = centroid_df.drop(columns=["sample_id", "radiological_type"]).to_numpy()
        labels = centroid_df["radiological_type"].map(label_map).to_numpy()
        clusters = KMeans(n_clusters=min(3, centroid_df.shape[0]), n_init=50, random_state=seed).fit_predict(
            centroid_matrix
        )

        rows.append(
            {
                "representation": representation,
                "compartment": compartment,
                "n_cells": int(mask.sum()),
                "n_samples": int(centroid_df.shape[0]),
                "ari": float(adjusted_rand_score(labels, clusters)),
                "nmi": float(normalized_mutual_info_score(labels, clusters)),
                "nn_accuracy": float(nearest_neighbor_accuracy(centroid_matrix, labels)),
                "silhouette": float(silhouette_score(centroid_matrix, labels)),
            }
        )
    return pd.DataFrame(rows)


def render_plot(metrics: pd.DataFrame, output_path: Path) -> None:
    plot_df = metrics.melt(
        id_vars=["representation", "compartment", "n_cells", "n_samples"],
        value_vars=["ari", "nmi", "nn_accuracy", "silhouette"],
        var_name="metric",
        value_name="value",
    )
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for ax, metric in zip(axes.flat, ["ari", "nmi", "nn_accuracy", "silhouette"]):
        sns.barplot(
            data=plot_df[plot_df["metric"] == metric],
            x="compartment",
            y="value",
            hue="representation",
            order=COMPARTMENT_ORDER,
            ax=ax,
        )
        ax.set_title(metric)
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=45)
        ax.legend(frameon=False)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    ensure_dir(args.output_dir)
    ensure_dir(args.output_dir / "figures")

    adata = ad.read_h5ad(args.adata)
    stack_adata = ad.read_h5ad(args.stack_embeddings)
    compartments = pd.read_csv(args.compartments).set_index(stack_adata.obs_names)["compartment"]

    stack_metrics = evaluate(
        stack_adata.obs,
        np.asarray(stack_adata.X),
        compartments,
        representation="stack",
        seed=args.seed,
    )
    pca_metrics = evaluate(
        stack_adata.obs,
        build_pca_embedding(adata),
        compartments,
        representation="pca",
        seed=args.seed,
    )
    metrics = pd.concat([stack_metrics, pca_metrics], ignore_index=True)

    metrics.to_csv(args.output_dir / "radiology_compartment_metrics.csv", index=False)
    (args.output_dir / "radiology_compartment_metrics.json").write_text(
        json.dumps(metrics.to_dict(orient="records"), indent=2) + "\n",
        encoding="utf-8",
    )
    render_plot(metrics, args.output_dir / "figures" / "radiology_compartment_comparison.png")


if __name__ == "__main__":
    main()

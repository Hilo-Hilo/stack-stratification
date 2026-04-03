#!/usr/bin/env python3
"""Benchmark Stack against a PCA baseline for LUAD radiological type."""

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
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.neighbors import NearestNeighbors

matplotlib.use("Agg")
os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")

REPO_ROOT = Path(__file__).resolve().parent.parent


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
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "results/luad_radiology_benchmark",
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


def centroid_metrics(obs: pd.DataFrame, embeddings: np.ndarray, seed: int) -> tuple[pd.DataFrame, dict[str, float]]:
    centroid_df = (
        pd.DataFrame(embeddings)
        .assign(
            sample_id=obs["sample_id"].values,
            radiological_type=obs["radiological_type"].values,
        )
        .groupby(["sample_id", "radiological_type"], as_index=False, observed=True)
        .mean()
    )
    centroid_matrix = centroid_df.drop(columns=["sample_id", "radiological_type"]).to_numpy()
    label_order = sorted(centroid_df["radiological_type"].unique())
    label_map = {label: idx for idx, label in enumerate(label_order)}
    labels = centroid_df["radiological_type"].map(label_map).to_numpy()
    clusters = KMeans(n_clusters=3, n_init=50, random_state=seed).fit_predict(centroid_matrix)

    metrics = {
        "ari": float(adjusted_rand_score(labels, clusters)),
        "nmi": float(normalized_mutual_info_score(labels, clusters)),
        "nn_accuracy": float(nearest_neighbor_accuracy(centroid_matrix, labels)),
        "silhouette": float(silhouette_score(centroid_matrix, labels)),
    }
    return centroid_df, metrics


def render_centroid_plot(stack_centroids: pd.DataFrame, pca_centroids: pd.DataFrame, output_path: Path) -> None:
    stack_coords = PCA(n_components=2, random_state=7).fit_transform(
        stack_centroids.drop(columns=["sample_id", "radiological_type"]).to_numpy()
    )
    pca_coords = PCA(n_components=2, random_state=7).fit_transform(
        pca_centroids.drop(columns=["sample_id", "radiological_type"]).to_numpy()
    )

    stack_plot = stack_centroids[["sample_id", "radiological_type"]].copy()
    stack_plot["x"] = stack_coords[:, 0]
    stack_plot["y"] = stack_coords[:, 1]

    pca_plot = pca_centroids[["sample_id", "radiological_type"]].copy()
    pca_plot["x"] = pca_coords[:, 0]
    pca_plot["y"] = pca_coords[:, 1]

    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    for ax, plot_df, title in [
        (axes[0], stack_plot, "Stack sample centroids"),
        (axes[1], pca_plot, "PCA sample centroids"),
    ]:
        sns.scatterplot(data=plot_df, x="x", y="y", hue="radiological_type", style="sample_id", s=160, ax=ax)
        for row in plot_df.itertuples():
            ax.text(row.x + 0.02, row.y + 0.02, row.sample_id, fontsize=9)
        ax.set_title(title)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.grid(alpha=0.15)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def render_metric_plot(metrics_df: pd.DataFrame, output_path: Path) -> None:
    plot_df = metrics_df.melt(id_vars=["representation"], var_name="metric", value_name="value")
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    sns.barplot(data=plot_df, x="metric", y="value", hue="representation", ax=ax)
    ax.set_title("Radiological-type benchmark: Stack versus PCA")
    ax.set_xlabel("")
    ax.set_ylabel("Score")
    ax.legend(frameon=False)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    ensure_dir(args.output_dir)
    ensure_dir(args.output_dir / "figures")

    adata = ad.read_h5ad(args.adata)
    stack_adata = ad.read_h5ad(args.stack_embeddings)
    pca_embedding = build_pca_embedding(adata)

    stack_centroids, stack_metrics = centroid_metrics(stack_adata.obs, np.asarray(stack_adata.X), args.seed)
    pca_centroids, pca_metrics = centroid_metrics(stack_adata.obs, pca_embedding, args.seed)

    metrics_df = pd.DataFrame(
        [
            {"representation": "stack", **stack_metrics},
            {"representation": "pca", **pca_metrics},
        ]
    )

    stack_centroids.to_csv(args.output_dir / "stack_sample_centroids.csv", index=False)
    pca_centroids.to_csv(args.output_dir / "pca_sample_centroids.csv", index=False)
    metrics_df.to_csv(args.output_dir / "radiology_metrics.csv", index=False)
    (args.output_dir / "radiology_metrics.json").write_text(
        json.dumps(metrics_df.to_dict(orient="records"), indent=2) + "\n",
        encoding="utf-8",
    )

    render_centroid_plot(stack_centroids, pca_centroids, args.output_dir / "figures" / "radiology_centroids.png")
    render_metric_plot(metrics_df, args.output_dir / "figures" / "radiology_metric_comparison.png")


if __name__ == "__main__":
    main()

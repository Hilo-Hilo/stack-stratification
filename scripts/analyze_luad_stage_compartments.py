#!/usr/bin/env python3
"""Compare LUAD stage concordance across coarse cell compartments."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.neighbors import NearestNeighbors

matplotlib.use("Agg")
os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")

REPO_ROOT = Path(__file__).resolve().parent.parent
HISTOLOGY_CODE = {"AIS": 0, "MIA": 1, "IAC": 2}
COMPARTMENT_MARKERS = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "MSLN"],
    "immune": ["PTPRC", "LST1", "TYROBP", "C1QB", "HLA-DRA"],
    "stromal": ["COL1A1", "COL1A2", "DCN", "TAGLN", "RGS5"],
    "endothelial": ["PECAM1", "VWF", "KDR", "EMCN", "RAMP2"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--adata",
        type=Path,
        default=REPO_ROOT / "runs/luad_stage_benchmark_300/gse189357_stage_subset.h5ad",
    )
    parser.add_argument(
        "--embeddings",
        type=Path,
        default=REPO_ROOT / "runs/luad_stage_benchmark_300/gse189357_stage_subset_stack_embeddings.h5ad",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "results/luad_stage_compartments",
    )
    parser.add_argument("--seed", type=int, default=7)
    return parser.parse_args()


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def compute_marker_scores(adata: ad.AnnData) -> pd.DataFrame:
    full = adata.X
    if not isinstance(full, np.ndarray):
        library_size = np.asarray(full.sum(axis=1)).ravel()
    else:
        library_size = full.sum(axis=1)
    library_size[library_size == 0] = 1

    scores: dict[str, np.ndarray] = {}
    for compartment, genes in COMPARTMENT_MARKERS.items():
        available = [gene for gene in genes if gene in adata.var_names]
        matrix = adata[:, available].X
        if not isinstance(matrix, np.ndarray):
            matrix = matrix.toarray()
        scores[compartment] = np.log1p((matrix / library_size[:, None]) * 1e4).mean(axis=1)
    return pd.DataFrame(scores, index=adata.obs_names)


def nearest_neighbor_accuracy(features: np.ndarray, labels: np.ndarray) -> float:
    nn = NearestNeighbors(n_neighbors=2, metric="euclidean")
    nn.fit(features)
    indices = nn.kneighbors(return_distance=False)
    return float(np.mean(labels[indices[:, 1]] == labels))


def compartment_metrics(
    embedding_adata: ad.AnnData,
    compartments: pd.Series,
    seed: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    ordered = ["all", *COMPARTMENT_MARKERS.keys()]
    for compartment in ordered:
        mask = np.ones(embedding_adata.n_obs, dtype=bool) if compartment == "all" else (
            compartments.values == compartment
        )
        subset_obs = embedding_adata.obs.loc[mask].copy()
        if subset_obs.empty:
            continue

        subset_embeddings = np.asarray(embedding_adata.X[mask])
        centroid_df = (
            pd.DataFrame(
                subset_embeddings,
                columns=[f"latent_{idx}" for idx in range(subset_embeddings.shape[1])],
            )
            .assign(
                sample_id=subset_obs["sample_id"].values,
                histology_type=subset_obs["histology_type"].values,
            )
            .groupby(["sample_id", "histology_type"], as_index=False, observed=True)
            .mean()
        )
        if centroid_df.shape[0] < 4:
            continue

        centroid_matrix = centroid_df.filter(regex=r"^latent_").to_numpy()
        labels = centroid_df["histology_type"].map(HISTOLOGY_CODE).to_numpy()
        kmeans = KMeans(n_clusters=min(3, centroid_df.shape[0]), n_init=50, random_state=seed)
        cluster_labels = kmeans.fit_predict(centroid_matrix)

        silhouette = np.nan
        if len(np.unique(labels)) > 1 and centroid_df.shape[0] > len(np.unique(labels)):
            silhouette = silhouette_score(centroid_matrix, labels, metric="euclidean")

        rows.append(
            {
                "compartment": compartment,
                "n_cells": int(mask.sum()),
                "n_samples": int(centroid_df.shape[0]),
                "ari": float(adjusted_rand_score(labels, cluster_labels)),
                "nmi": float(normalized_mutual_info_score(labels, cluster_labels)),
                "nn_accuracy": float(nearest_neighbor_accuracy(centroid_matrix, labels)),
                "silhouette": float(silhouette),
            }
        )

    return pd.DataFrame(rows)


def render_composition_plot(summary: pd.DataFrame, output_path: Path) -> None:
    ordered_samples = (
        summary[["sample_id", "histology_type"]]
        .drop_duplicates()
        .sort_values(["histology_type", "sample_id"])
        ["sample_id"]
        .tolist()
    )
    pivot = (
        summary.pivot(index="sample_id", columns="compartment", values="fraction")
        .fillna(0)
        .reindex(ordered_samples)
    )

    fig, ax = plt.subplots(figsize=(10, 5.5), constrained_layout=True)
    bottom = np.zeros(len(pivot))
    for compartment in COMPARTMENT_MARKERS:
        values = pivot.get(compartment, pd.Series(0, index=pivot.index)).to_numpy()
        ax.bar(pivot.index, values, bottom=bottom, label=compartment)
        bottom += values
    ax.set_title("Marker-based compartment composition by sample")
    ax.set_ylabel("Fraction of cells")
    ax.set_xlabel("Sample")
    ax.tick_params(axis="x", rotation=45)
    ax.legend(frameon=False)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def render_metric_plot(metrics: pd.DataFrame, output_path: Path) -> None:
    plot_df = metrics.melt(
        id_vars=["compartment", "n_cells", "n_samples"],
        value_vars=["ari", "nmi", "nn_accuracy", "silhouette"],
        var_name="metric",
        value_name="value",
    )
    fig, ax = plt.subplots(figsize=(10, 5.5), constrained_layout=True)
    sns.barplot(data=plot_df, x="metric", y="value", hue="compartment", ax=ax)
    ax.set_title("Stage concordance metrics by compartment")
    ax.set_xlabel("")
    ax.set_ylabel("Score")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    ensure_dir(args.output_dir)
    ensure_dir(args.output_dir / "figures")

    adata = ad.read_h5ad(args.adata)
    embedding_adata = ad.read_h5ad(args.embeddings)

    score_df = compute_marker_scores(adata)
    compartments = score_df.idxmax(axis=1).rename("compartment")

    cell_table = embedding_adata.obs.copy()
    cell_table["compartment"] = compartments.values
    for compartment in COMPARTMENT_MARKERS:
        cell_table[f"{compartment}_score"] = score_df[compartment].values

    summary = (
        cell_table.groupby(["sample_id", "histology_type", "compartment"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    totals = summary.groupby("sample_id", observed=True)["n_cells"].transform("sum")
    summary["fraction"] = summary["n_cells"] / totals

    metrics = compartment_metrics(embedding_adata, compartments, seed=args.seed)

    cell_table.to_csv(args.output_dir / "cell_compartments.csv", index=False)
    summary.to_csv(args.output_dir / "sample_compartment_summary.csv", index=False)
    metrics.to_csv(args.output_dir / "compartment_metrics.csv", index=False)
    (args.output_dir / "compartment_metrics.json").write_text(
        json.dumps(metrics.to_dict(orient="records"), indent=2) + "\n",
        encoding="utf-8",
    )

    render_composition_plot(summary, args.output_dir / "figures" / "compartment_composition_by_sample.png")
    render_metric_plot(metrics, args.output_dir / "figures" / "compartment_metric_comparison.png")


if __name__ == "__main__":
    main()

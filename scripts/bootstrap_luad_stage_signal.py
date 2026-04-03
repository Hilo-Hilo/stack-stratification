#!/usr/bin/env python3
"""Bootstrap stage-concordance metrics for Stack and PCA representations."""

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
HISTOLOGY_CODE = {"AIS": 0, "MIA": 1, "IAC": 2}
COMPARTMENTS = ["all", "epithelial", "immune"]


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
        default=REPO_ROOT / "results/luad_stage_bootstrap",
    )
    parser.add_argument("--n-boot", type=int, default=200)
    parser.add_argument("--seed", type=int, default=7)
    return parser.parse_args()


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def nearest_neighbor_accuracy(features: np.ndarray, labels: np.ndarray) -> float:
    nn = NearestNeighbors(n_neighbors=2, metric="euclidean")
    nn.fit(features)
    indices = nn.kneighbors(return_distance=False)
    return float(np.mean(labels[indices[:, 1]] == labels))


def build_pca_embedding(adata: ad.AnnData) -> np.ndarray:
    work = adata.copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    sc.pp.highly_variable_genes(work, n_top_genes=2000, flavor="seurat")
    work = work[:, work.var["highly_variable"].values].copy()
    sc.pp.pca(work, n_comps=50)
    return np.asarray(work.obsm["X_pca"])


def bootstrap_representation(
    obs: pd.DataFrame,
    embeddings: np.ndarray,
    compartments: pd.Series,
    representation: str,
    n_boot: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []

    for compartment in COMPARTMENTS:
        mask = np.ones(obs.shape[0], dtype=bool) if compartment == "all" else (
            compartments.values == compartment
        )
        subset_obs = obs.loc[mask].copy()
        subset_embeddings = embeddings[mask]
        if subset_obs.empty:
            continue

        sample_groups = {}
        for sample_id, frame in subset_obs.groupby("sample_id", observed=True):
            sample_groups[sample_id] = frame.index.to_numpy()

        sample_labels = (
            subset_obs[["sample_id", "histology_type"]]
            .drop_duplicates()
            .sort_values("sample_id")
            .set_index("sample_id")["histology_type"]
        )
        ordered_samples = sample_labels.index.tolist()
        label_codes = sample_labels.map(HISTOLOGY_CODE).to_numpy()

        index_lookup = pd.Series(np.arange(subset_obs.shape[0]), index=subset_obs.index)

        for bootstrap_id in range(n_boot):
            centroids = []
            for sample_id in ordered_samples:
                sample_indices = index_lookup.loc[sample_groups[sample_id]].to_numpy()
                resampled = rng.choice(sample_indices, size=len(sample_indices), replace=True)
                centroids.append(subset_embeddings[resampled].mean(axis=0))
            centroid_matrix = np.vstack(centroids)
            clusters = KMeans(n_clusters=3, n_init=25, random_state=seed + bootstrap_id).fit_predict(
                centroid_matrix
            )
            rows.append(
                {
                    "representation": representation,
                    "compartment": compartment,
                    "bootstrap_id": bootstrap_id,
                    "ari": float(adjusted_rand_score(label_codes, clusters)),
                    "nmi": float(normalized_mutual_info_score(label_codes, clusters)),
                    "nn_accuracy": float(nearest_neighbor_accuracy(centroid_matrix, label_codes)),
                    "silhouette": float(silhouette_score(centroid_matrix, label_codes)),
                }
            )

    return pd.DataFrame(rows)


def render_bootstrap_plot(results: pd.DataFrame, output_path: Path) -> None:
    plot_df = results.melt(
        id_vars=["representation", "compartment", "bootstrap_id"],
        value_vars=["ari", "nmi", "nn_accuracy", "silhouette"],
        var_name="metric",
        value_name="value",
    )
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for ax, metric in zip(axes.flat, ["ari", "nmi", "nn_accuracy", "silhouette"]):
        sns.boxplot(
            data=plot_df[plot_df["metric"] == metric],
            x="compartment",
            y="value",
            hue="representation",
            order=COMPARTMENTS,
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
    compartment_series = pd.read_csv(args.compartments).set_index(stack_adata.obs_names)["compartment"]

    stack_results = bootstrap_representation(
        stack_adata.obs,
        np.asarray(stack_adata.X),
        compartment_series,
        representation="stack",
        n_boot=args.n_boot,
        seed=args.seed,
    )
    pca_results = bootstrap_representation(
        stack_adata.obs,
        build_pca_embedding(adata),
        compartment_series,
        representation="pca",
        n_boot=args.n_boot,
        seed=args.seed,
    )

    results = pd.concat([stack_results, pca_results], ignore_index=True)
    summary = (
        results.groupby(["representation", "compartment"], observed=True)[
            ["ari", "nmi", "nn_accuracy", "silhouette"]
        ]
        .agg(["median", "mean"])
    )
    summary.columns = [f"{metric}_{stat}" for metric, stat in summary.columns]
    summary = summary.reset_index()

    results.to_csv(args.output_dir / "bootstrap_metrics.csv", index=False)
    summary.to_csv(args.output_dir / "bootstrap_metric_summary.csv")
    (args.output_dir / "bootstrap_metric_summary.json").write_text(
        json.dumps(
            summary.to_dict(orient="records"),
            indent=2,
            default=lambda x: float(x) if isinstance(x, np.floating) else x,
        )
        + "\n",
        encoding="utf-8",
    )
    render_bootstrap_plot(results, args.output_dir / "figures" / "bootstrap_metric_distributions.png")


if __name__ == "__main__":
    main()

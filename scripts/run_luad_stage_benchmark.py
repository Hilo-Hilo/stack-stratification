#!/usr/bin/env python3
"""Run a first-pass LUAD pathology-stage benchmark with Stack embeddings."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import sys
import tarfile
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import h5py
from scipy.io import mmread
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.neighbors import NearestNeighbors

matplotlib.use("Agg")
os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")

REPO_ROOT = Path(__file__).resolve().parent.parent


SAMPLE_METADATA = {
    "TD1": {"histology_type": "IAC", "radiological_type": "SN", "gender": "Female"},
    "TD2": {"histology_type": "IAC", "radiological_type": "SN", "gender": "Male"},
    "TD3": {"histology_type": "MIA", "radiological_type": "pGGN", "gender": "Male"},
    "TD4": {"histology_type": "MIA", "radiological_type": "pGGN", "gender": "Female"},
    "TD5": {"histology_type": "AIS", "radiological_type": "SSN", "gender": "Female"},
    "TD6": {"histology_type": "MIA", "radiological_type": "SSN", "gender": "Female"},
    "TD7": {"histology_type": "AIS", "radiological_type": "SSN", "gender": "Female"},
    "TD8": {"histology_type": "AIS", "radiological_type": "pGGN", "gender": "Male"},
    "TD9": {"histology_type": "IAC", "radiological_type": "SSN", "gender": "Female"},
}

HISTOLOGY_ORDER = ["AIS", "MIA", "IAC"]
HISTOLOGY_CODE = {label: idx for idx, label in enumerate(HISTOLOGY_ORDER)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--geo-tar",
        type=Path,
        default=REPO_ROOT / "data/raw/gse189357/GSE189357_RAW.tar",
        help="Path to the downloaded GEO tarball.",
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=REPO_ROOT / "runs/luad_stage_benchmark",
        help="Directory for intermediate artifacts, logs, and results.",
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=REPO_ROOT / "models/Stack-Large/bc_large.ckpt",
        help="Stack checkpoint path.",
    )
    parser.add_argument(
        "--genelist",
        type=Path,
        default=REPO_ROOT / "models/Stack-Large/basecount_1000per_15000max.pkl",
        help="Stack gene list path.",
    )
    parser.add_argument(
        "--max-cells-per-sample",
        type=int,
        default=300,
        help="Cap cells per sample for a controlled CPU run.",
    )
    parser.add_argument(
        "--min-genes",
        type=int,
        default=200,
        help="Filter cells with fewer detected genes than this threshold.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=8,
        help="Batch size passed to stack-embedding.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=7,
        help="Random seed used for subsampling and clustering.",
    )
    parser.add_argument(
        "--skip-embedding",
        action="store_true",
        help="Skip Stack inference and only run downstream analysis on existing artifacts.",
    )
    return parser.parse_args()


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def h5ad_has_groups(path: Path, required_groups: tuple[str, ...]) -> bool:
    if not path.exists():
        return False
    try:
        with h5py.File(path, "r") as handle:
            return all(group in handle for group in required_groups)
    except OSError:
        return False


def unpack_geo_bundle(tar_path: Path, extract_dir: Path) -> None:
    sentinel = extract_dir / ".complete"
    if sentinel.exists():
        return
    ensure_dir(extract_dir)
    with tarfile.open(tar_path) as archive:
        archive.extractall(extract_dir)
    sentinel.write_text("ok\n", encoding="utf-8")


def read_10x_triplet(prefix: str, source_dir: Path) -> ad.AnnData:
    matrix_path = source_dir / f"{prefix}_matrix.mtx.gz"
    barcodes_path = source_dir / f"{prefix}_barcodes.tsv.gz"
    features_path = source_dir / f"{prefix}_features.tsv.gz"

    with gzip.open(matrix_path, "rb") as handle:
        matrix = mmread(handle).tocsr().transpose().tocsr()
    barcodes = pd.read_csv(barcodes_path, header=None, sep="\t")
    features = pd.read_csv(features_path, header=None, sep="\t")

    obs = pd.DataFrame(index=barcodes.iloc[:, 0].astype(str).values)
    obs.index = pd.Index(obs.index.astype(str))
    obs.index.name = None
    var = pd.DataFrame(
        {
            "gene_id": features.iloc[:, 0].astype(str).values,
            "feature_name": features.iloc[:, 1].astype(str).values,
            "feature_type": features.iloc[:, 2].astype(str).values,
        },
        index=features.iloc[:, 1].astype(str).values,
    )
    var.index = pd.Index(var.index.astype(str))
    var.index.name = None

    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.var_names_make_unique()
    return adata


def filter_and_subsample(adata: ad.AnnData, sample_id: str, max_cells: int, min_genes: int, seed: int) -> ad.AnnData:
    gene_counts = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    keep_mask = gene_counts >= min_genes
    adata = adata[keep_mask].copy()

    if adata.n_obs > max_cells:
        rng = np.random.default_rng(seed)
        selected = np.sort(rng.choice(adata.n_obs, size=max_cells, replace=False))
        adata = adata[selected].copy()

    metadata = SAMPLE_METADATA[sample_id]
    adata.obs["sample_id"] = sample_id
    adata.obs["histology_type"] = metadata["histology_type"]
    adata.obs["radiological_type"] = metadata["radiological_type"]
    adata.obs["gender"] = metadata["gender"]
    adata.obs["organism"] = "Homo sapiens"
    adata.obs["cohort"] = "GSE189357"
    return adata


def build_combined_adata(source_dir: Path, max_cells_per_sample: int, min_genes: int, seed: int) -> ad.AnnData:
    adatas: list[ad.AnnData] = []
    prefixes = sorted({path.name.replace("_barcodes.tsv.gz", "") for path in source_dir.glob("*_barcodes.tsv.gz")})
    for prefix in prefixes:
        sample_id = prefix.split("_")[1]
        adata = read_10x_triplet(prefix, source_dir)
        adata = filter_and_subsample(
            adata,
            sample_id=sample_id,
            max_cells=max_cells_per_sample,
            min_genes=min_genes,
            seed=seed,
        )
        adatas.append(adata)

    combined = ad.concat(adatas, join="inner", merge="same")
    combined.layers["counts"] = combined.X.copy()
    return combined


def run_stack_embedding(args: argparse.Namespace, combined_h5ad: Path, embeddings_h5ad: Path) -> None:
    cmd = [
        str(Path(sys.prefix) / "bin/stack-embedding"),
        "--checkpoint",
        str(args.checkpoint),
        "--adata",
        str(combined_h5ad),
        "--genelist",
        str(args.genelist),
        "--output",
        str(embeddings_h5ad),
        "--batch-size",
        str(args.batch_size),
        "--num-workers",
        "0",
        "--device",
        "cpu",
        "--show-progress",
        "--obs-source",
        str(combined_h5ad),
    ]
    env = os.environ.copy()
    env.setdefault("TMPDIR", str((args.workdir / "tmp").resolve()))
    env.setdefault("MPLCONFIGDIR", str((args.workdir / "matplotlib").resolve()))
    ensure_dir(Path(env["TMPDIR"]))
    ensure_dir(Path(env["MPLCONFIGDIR"]))
    subprocess.run(cmd, check=True, env=env)


def nearest_neighbor_accuracy(embeddings: np.ndarray, labels: np.ndarray) -> float:
    nn = NearestNeighbors(n_neighbors=2, metric="euclidean")
    nn.fit(embeddings)
    indices = nn.kneighbors(return_distance=False)
    preds = labels[indices[:, 1]]
    return float(np.mean(preds == labels))


def render_cell_embedding_figure(embedding_2d: np.ndarray, obs: pd.DataFrame, output_path: Path) -> None:
    plot_df = obs.copy()
    plot_df["dim1"] = embedding_2d[:, 0]
    plot_df["dim2"] = embedding_2d[:, 1]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)
    sns.scatterplot(
        data=plot_df,
        x="dim1",
        y="dim2",
        hue="histology_type",
        hue_order=HISTOLOGY_ORDER,
        s=18,
        linewidth=0,
        ax=axes[0],
    )
    axes[0].set_title("Stack cell embeddings colored by pathology label")
    axes[0].set_xlabel("Embedding axis 1")
    axes[0].set_ylabel("Embedding axis 2")

    sns.scatterplot(
        data=plot_df,
        x="dim1",
        y="dim2",
        hue="sample_id",
        s=18,
        linewidth=0,
        ax=axes[1],
    )
    axes[1].set_title("Same embedding colored by sample")
    axes[1].set_xlabel("Embedding axis 1")
    axes[1].set_ylabel("Embedding axis 2")

    for ax in axes:
        ax.grid(alpha=0.15)

    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def render_sample_centroid_figure(centroids: pd.DataFrame, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 5.5), constrained_layout=True)
    sns.scatterplot(
        data=centroids,
        x="pc1",
        y="pc2",
        hue="histology_type",
        hue_order=HISTOLOGY_ORDER,
        style="sample_id",
        s=160,
        ax=ax,
    )
    for row in centroids.itertuples():
        ax.text(row.pc1 + 0.02, row.pc2 + 0.02, row.sample_id, fontsize=9)
    ax.set_title("Sample centroids in Stack embedding space")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.grid(alpha=0.15)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def render_distance_heatmap(centroids: pd.DataFrame, output_path: Path) -> None:
    ordered = centroids.sort_values(["histology_type", "sample_id"]).reset_index(drop=True)
    sample_labels = [f"{sid} ({label})" for sid, label in zip(ordered["sample_id"], ordered["histology_type"])]
    centroid_matrix = ordered.filter(regex=r"^latent_").to_numpy()
    distances = euclidean_distances(centroid_matrix)

    fig, ax = plt.subplots(figsize=(7.5, 6.5), constrained_layout=True)
    sns.heatmap(
        distances,
        xticklabels=sample_labels,
        yticklabels=sample_labels,
        cmap="mako",
        square=True,
        ax=ax,
    )
    ax.set_title("Sample-to-sample distances in Stack space")
    ax.tick_params(axis="x", rotation=45)
    ax.tick_params(axis="y", rotation=0)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def analyze_embeddings(embeddings_h5ad: Path, results_dir: Path) -> dict[str, object]:
    embedding_adata = ad.read_h5ad(embeddings_h5ad)
    cell_embeddings = np.asarray(embedding_adata.X)
    obs = embedding_adata.obs.copy()

    perplexity = max(5, min(30, (cell_embeddings.shape[0] - 1) // 3))
    tsne = TSNE(n_components=2, random_state=7, init="pca", learning_rate="auto", perplexity=perplexity)
    cell_embedding_2d = tsne.fit_transform(cell_embeddings)

    render_cell_embedding_figure(
        cell_embedding_2d,
        obs,
        results_dir / "figures" / "cell_embeddings_by_label.png",
    )

    centroid_df = (
        pd.DataFrame(cell_embeddings, columns=[f"latent_{idx}" for idx in range(cell_embeddings.shape[1])])
        .assign(sample_id=obs["sample_id"].values, histology_type=obs["histology_type"].values)
        .groupby(["sample_id", "histology_type"], as_index=False, observed=True)
        .mean()
    )

    centroid_matrix = centroid_df.filter(regex=r"^latent_").to_numpy()
    centroid_pca = PCA(n_components=2, random_state=7).fit_transform(centroid_matrix)
    centroid_df["pc1"] = centroid_pca[:, 0]
    centroid_df["pc2"] = centroid_pca[:, 1]

    kmeans = KMeans(n_clusters=3, n_init=50, random_state=7)
    centroid_df["cluster"] = kmeans.fit_predict(centroid_matrix)

    label_codes = centroid_df["histology_type"].map(HISTOLOGY_CODE).to_numpy()
    ari = adjusted_rand_score(label_codes, centroid_df["cluster"].to_numpy())
    nmi = normalized_mutual_info_score(label_codes, centroid_df["cluster"].to_numpy())
    nn_acc = nearest_neighbor_accuracy(centroid_matrix, label_codes)
    silhouette = float(silhouette_score(centroid_matrix, label_codes, metric="euclidean"))

    render_sample_centroid_figure(centroid_df, results_dir / "figures" / "sample_centroids.png")
    render_distance_heatmap(centroid_df, results_dir / "figures" / "sample_distance_heatmap.png")

    centroid_df.to_csv(results_dir / "sample_centroids.csv", index=False)

    metrics = {
        "n_cells": int(cell_embeddings.shape[0]),
        "embedding_dim": int(cell_embeddings.shape[1]),
        "n_samples": int(centroid_df.shape[0]),
        "ari_sample_clusters_vs_histology": float(ari),
        "nmi_sample_clusters_vs_histology": float(nmi),
        "nearest_neighbor_histology_accuracy": float(nn_acc),
        "silhouette_histology_on_sample_centroids": float(silhouette),
        "histology_counts": centroid_df["histology_type"].value_counts().sort_index().to_dict(),
    }

    (results_dir / "metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return metrics


def main() -> None:
    args = parse_args()
    ensure_dir(args.workdir)
    ensure_dir(args.workdir / "figures")
    extract_dir = args.workdir / "gse189357_extracted"
    combined_h5ad = args.workdir / "gse189357_stage_subset.h5ad"
    embeddings_h5ad = args.workdir / "gse189357_stage_subset_stack_embeddings.h5ad"

    unpack_geo_bundle(args.geo_tar, extract_dir)

    if not h5ad_has_groups(combined_h5ad, ("X", "obs", "var")):
        if combined_h5ad.exists():
            combined_h5ad.unlink()
        combined = build_combined_adata(
            extract_dir,
            max_cells_per_sample=args.max_cells_per_sample,
            min_genes=args.min_genes,
            seed=args.seed,
        )
        combined.write_h5ad(combined_h5ad)

    if not args.skip_embedding and not h5ad_has_groups(embeddings_h5ad, ("X", "obs", "var")):
        if embeddings_h5ad.exists():
            embeddings_h5ad.unlink()
        run_stack_embedding(args, combined_h5ad, embeddings_h5ad)

    metrics = analyze_embeddings(embeddings_h5ad, args.workdir)
    print(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Consolidate current LUAD benchmark results into summary tables and figures."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

matplotlib.use("Agg")
os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")

REPO_ROOT = Path(__file__).resolve().parent.parent
METRICS = ["ari", "nmi", "nn_accuracy", "silhouette"]
METRIC_LABELS = {
    "ari": "ARI",
    "nmi": "NMI",
    "nn_accuracy": "1-NN",
    "silhouette": "Silhouette",
}
BENCHMARK_LABELS = {
    "stage": "Pathology stage",
    "radiology": "Radiology",
}
COMPARTMENT_ORDER = ["all", "immune", "epithelial", "stromal", "endothelial"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage-metrics",
        type=Path,
        default=REPO_ROOT / "results/luad_stage_benchmark_300/metrics.json",
    )
    parser.add_argument(
        "--stage-representation",
        type=Path,
        default=REPO_ROOT / "results/luad_stage_representation_comparison/representation_metrics.csv",
    )
    parser.add_argument(
        "--stage-bootstrap",
        type=Path,
        default=REPO_ROOT / "results/luad_stage_bootstrap/bootstrap_metric_summary.csv",
    )
    parser.add_argument(
        "--radiology-metrics",
        type=Path,
        default=REPO_ROOT / "results/luad_radiology_benchmark/radiology_metrics.csv",
    )
    parser.add_argument(
        "--radiology-bootstrap",
        type=Path,
        default=REPO_ROOT / "results/luad_radiology_bootstrap/radiology_bootstrap_summary.csv",
    )
    parser.add_argument(
        "--radiology-compartments",
        type=Path,
        default=REPO_ROOT / "results/luad_radiology_compartments/radiology_compartment_metrics.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "results/luad_benchmark_summary",
    )
    return parser.parse_args()


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def load_stage_point_metrics(metrics_path: Path, representation_path: Path) -> pd.DataFrame:
    payload = json.loads(metrics_path.read_text(encoding="utf-8"))
    rows = [
        {
            "benchmark": "stage",
            "representation": "stack",
            "estimate": "point",
            "compartment": "all",
            "ari": payload["ari_sample_clusters_vs_histology"],
            "nmi": payload["nmi_sample_clusters_vs_histology"],
            "nn_accuracy": payload["nearest_neighbor_histology_accuracy"],
            "silhouette": payload["silhouette_histology_on_sample_centroids"],
        }
    ]

    representation_df = pd.read_csv(representation_path)
    pca_all = representation_df[
        (representation_df["representation"] == "pca") & (representation_df["compartment"] == "all")
    ].iloc[0]
    rows.append(
        {
            "benchmark": "stage",
            "representation": "pca",
            "estimate": "point",
            "compartment": "all",
            "ari": pca_all["ari"],
            "nmi": pca_all["nmi"],
            "nn_accuracy": pca_all["nn_accuracy"],
            "silhouette": pca_all["silhouette"],
        }
    )
    return pd.DataFrame(rows)


def load_stage_bootstrap(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df[df["compartment"] == "all"].copy()
    df["benchmark"] = "stage"
    df["estimate"] = "bootstrap_median"
    df["compartment"] = "all"
    rename_map = {f"{metric}_median": metric for metric in METRICS}
    return df.rename(columns=rename_map)[["benchmark", "representation", "estimate", "compartment", *METRICS]]


def load_radiology_point(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path).copy()
    df["benchmark"] = "radiology"
    df["estimate"] = "point"
    df["compartment"] = "all"
    return df[["benchmark", "representation", "estimate", "compartment", *METRICS]]


def load_radiology_bootstrap(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path).copy()
    df["benchmark"] = "radiology"
    df["estimate"] = "bootstrap_median"
    df["compartment"] = "all"
    rename_map = {f"{metric}_median": metric for metric in METRICS}
    return df.rename(columns=rename_map)[["benchmark", "representation", "estimate", "compartment", *METRICS]]


def load_compartment_tables(stage_representation_path: Path, radiology_path: Path) -> pd.DataFrame:
    stage_df = pd.read_csv(stage_representation_path).copy()
    stage_df["benchmark"] = "stage"
    radiology_df = pd.read_csv(radiology_path).copy()
    radiology_df["benchmark"] = "radiology"
    combined = pd.concat([stage_df, radiology_df], ignore_index=True)
    return combined[["benchmark", "representation", "compartment", *METRICS, "n_cells", "n_samples"]]


def compute_gain_table(summary_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for benchmark in sorted(summary_df["benchmark"].unique()):
        for estimate in sorted(summary_df["estimate"].unique()):
            subset = summary_df[
                (summary_df["benchmark"] == benchmark)
                & (summary_df["estimate"] == estimate)
                & (summary_df["compartment"] == "all")
            ]
            if set(subset["representation"]) != {"pca", "stack"}:
                continue
            pivot = subset.set_index("representation")[METRICS]
            rows.append(
                {
                    "benchmark": benchmark,
                    "estimate": estimate,
                    **{metric: float(pivot.loc["stack", metric] - pivot.loc["pca", metric]) for metric in METRICS},
                }
            )
    return pd.DataFrame(rows)


def render_overview(summary_df: pd.DataFrame, gain_df: pd.DataFrame, compartments_df: pd.DataFrame, output_path: Path) -> None:
    summary_plot = summary_df[summary_df["compartment"] == "all"].melt(
        id_vars=["benchmark", "representation", "estimate"],
        value_vars=["ari", "nmi"],
        var_name="metric",
        value_name="value",
    )
    summary_plot["panel"] = summary_plot["estimate"].map(
        {"point": "Point estimate", "bootstrap_median": "Bootstrap median"}
    )
    summary_plot["benchmark_label"] = summary_plot["benchmark"].map(BENCHMARK_LABELS)
    summary_plot["metric_label"] = summary_plot["metric"].map(METRIC_LABELS)
    summary_plot["benchmark_metric"] = (
        summary_plot["benchmark_label"] + "\n" + summary_plot["metric_label"]
    )

    gain_plot = gain_df.melt(
        id_vars=["benchmark", "estimate"],
        value_vars=METRICS,
        var_name="metric",
        value_name="delta",
    )
    gain_heatmap = gain_plot.pivot(
        index="metric",
        columns=["benchmark", "estimate"],
        values="delta",
    ).reindex(index=METRICS)
    gain_heatmap.index = [METRIC_LABELS[metric] for metric in gain_heatmap.index]
    gain_heatmap.columns = [
        f"{BENCHMARK_LABELS[benchmark]}\n{'Point' if estimate == 'point' else 'Bootstrap'}"
        for benchmark, estimate in gain_heatmap.columns
    ]

    stack_compartments = compartments_df[compartments_df["representation"] == "stack"].copy()
    stack_heatmap = stack_compartments.pivot(
        index="compartment",
        columns="benchmark",
        values="nmi",
    ).reindex(COMPARTMENT_ORDER)
    stack_heatmap = stack_heatmap.dropna(how="all")
    stack_heatmap.index = [label.title() for label in stack_heatmap.index]
    stack_heatmap.columns = [BENCHMARK_LABELS[column] for column in stack_heatmap.columns]

    pca_compartments = compartments_df[compartments_df["representation"] == "pca"].copy()
    pca_heatmap = pca_compartments.pivot(
        index="compartment",
        columns="benchmark",
        values="nmi",
    ).reindex(COMPARTMENT_ORDER)
    pca_heatmap = pca_heatmap.dropna(how="all")
    pca_heatmap.index = [label.title() for label in pca_heatmap.index]
    pca_heatmap.columns = [BENCHMARK_LABELS[column] for column in pca_heatmap.columns]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)

    for ax, panel in zip(axes[0], ["Point estimate", "Bootstrap median"]):
        sns.barplot(
            data=summary_plot[summary_plot["panel"] == panel],
            x="benchmark_metric",
            y="value",
            hue="representation",
            ax=ax,
        )
        ax.set_title(panel)
        ax.set_xlabel("")
        ax.set_ylabel("Score")
        ax.tick_params(axis="x", rotation=0)
        ax.legend(frameon=False, title="")
        ax.grid(axis="y", alpha=0.15)

    sns.heatmap(gain_heatmap, annot=True, fmt=".3f", cmap="RdYlGn", center=0.0, ax=axes[1, 0])
    axes[1, 0].set_title("Stack minus PCA gain")
    axes[1, 0].set_xlabel("")
    axes[1, 0].set_ylabel("")

    sns.heatmap(stack_heatmap, annot=True, fmt=".3f", cmap="crest", vmin=0.0, vmax=0.65, ax=axes[1, 1])
    axes[1, 1].set_title("Stack compartment NMI")
    axes[1, 1].set_xlabel("")
    axes[1, 1].set_ylabel("")

    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
    sns.heatmap(pca_heatmap, annot=True, fmt=".3f", cmap="crest", vmin=0.0, vmax=0.65, ax=ax)
    ax.set_title("PCA compartment NMI")
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.savefig(output_path.with_name("benchmark_summary_pca_compartments.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    ensure_dir(args.output_dir)
    ensure_dir(args.output_dir / "figures")

    summary_df = pd.concat(
        [
            load_stage_point_metrics(args.stage_metrics, args.stage_representation),
            load_stage_bootstrap(args.stage_bootstrap),
            load_radiology_point(args.radiology_metrics),
            load_radiology_bootstrap(args.radiology_bootstrap),
        ],
        ignore_index=True,
    )
    gain_df = compute_gain_table(summary_df)
    compartments_df = load_compartment_tables(args.stage_representation, args.radiology_compartments)

    summary_df.to_csv(args.output_dir / "benchmark_summary.csv", index=False)
    gain_df.to_csv(args.output_dir / "benchmark_gain_summary.csv", index=False)
    compartments_df.to_csv(args.output_dir / "compartment_summary.csv", index=False)
    (args.output_dir / "benchmark_summary.json").write_text(
        json.dumps(
            {
                "summary": summary_df.to_dict(orient="records"),
                "gains": gain_df.to_dict(orient="records"),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    render_overview(
        summary_df,
        gain_df,
        compartments_df,
        args.output_dir / "figures" / "benchmark_summary_overview.png",
    )


if __name__ == "__main__":
    main()

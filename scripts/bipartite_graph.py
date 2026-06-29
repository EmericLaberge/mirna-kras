#!/usr/bin/env python3
"""
Bipartite graph visualization for miRNA–transcript binding affinities.

Reads the binding matrix from DuckDB, filters by strength, and renders
a bipartite graph (miRNAs left, genes right) with edges weighted by
mean |binding_dG|.

Usage:
    .venv/bin/python3 scripts/bipartite_graph.py \\
        --threshold 16 --top-mirnas 80 --top-edges 300 \\
        --output figures/bipartite_graph.png

    .venv/bin/python3 scripts/bipartite_graph.py \\
        --use-residual --threshold 0.5 --top-mirnas 80 --top-edges 300 \\
        --output figures/bipartite_graph_residual.png
"""

import warnings
warnings.filterwarnings("ignore")

import os
import argparse

import duckdb
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll
from matplotlib.patches import FancyArrowPatch


# ────────────────────────────────────────────────────────────────
# Constants
# ────────────────────────────────────────────────────────────────
TRANSCRIPT_GENE_MAP = {
    "ENST00000215832": "RAF1",
    "ENST00000226574": "TP53",
    "ENST00000227507": "ERBB2",
    "ENST00000251849": "PTEN",
    "ENST00000256078": "BRAF",
    "ENST00000259523": "MAP2K1",
    "ENST00000262948": "KRAS",
    "ENST00000263025": "NRAS",
    "ENST00000263826": "HRAS",
    "ENST00000263967": "MAP2K2",
    "ENST00000267163": "MAPK3",
    "ENST00000269305": "MAPK1",
    "ENST00000275493": "RB1",
    "ENST00000307102": "DUSP6",
    "ENST00000308639": "MTOR",
    "ENST00000311278": "PIK3CB",
    "ENST00000315204": "APC",
    "ENST00000318493": "NF1",
    "ENST00000316629": "PIK3CA",
    "ENST00000331710": "AKT1",
    "ENST00000349310": "AKT2",
    "ENST00000361445": "EGFR",
    "ENST00000366790": "MYC",
    "ENST00000371953": "SPRY2",
    "ENST00000377687": "CDKN1A",
}

KRAS_PATHWAY = {"KRAS", "NRAS", "HRAS", "PIK3CA", "PIK3CB", "MTOR",
                "BRAF", "RAF1", "MAP2K1", "MAP2K2", "MAPK1"}
TUMOR_SUPPRESSORS = {"TP53", "PTEN", "RB1"}
MAPK_PATHWAY = {"EGFR", "ERBB2"}

PATHWAY_COLORS = {
    "KRas": "#e74c3c",
    "Tumor suppressor": "#3498db",
    "MAPK": "#2ecc71",
    "Other": "#95a5a6",
}


def gene_to_pathway(gene: str) -> str:
    """Classify a gene symbol into one of KRas / Tumor suppressor / MAPK / Other.

    Args:
        gene: HGNC gene symbol (e.g. ``KRAS``).

    Returns:
        The pathway label.
    """
    if gene in KRAS_PATHWAY:
        return "KRas"
    if gene in TUMOR_SUPPRESSORS:
        return "Tumor suppressor"
    if gene in MAPK_PATHWAY:
        return "MAPK"
    return "Other"


# ────────────────────────────────────────────────────────────────
# Data loading
# ────────────────────────────────────────────────────────────────
def load_binding_matrix(
    db_path: str = "rimap_results.duckdb",
    exclude_seedless: bool = True,
) -> tuple[np.ndarray, np.ndarray, list[str], list[str]]:
    """Load binding matrix from DuckDB.

    Returns:
        V_raw: raw binding matrix (n_mirna × n_transcript)
        V_residual: residual matrix (col-mean subtracted, clipped, zero-filled)
        mirna_ids: list of miRNA identifiers
        transcript_ids: list of transcript identifiers
    """
    con = duckdb.connect(db_path, read_only=True)
    df = con.execute("SELECT * FROM results").fetchdf()
    con.close()

    n_total = len(df)
    print(f"  Total interactions: {n_total:,}")

    if exclude_seedless:
        df = df[df["seed_type"] != "seedless"]
        print(f"  After removing seedless: {len(df):,}")

    agg = (
        df.groupby(["mirna_id", "transcript_id"])
        .agg(weight=("binding_dG", lambda x: abs(x).mean()))
        .reset_index()
    )

    matrix_df = agg.pivot(
        index="mirna_id", columns="transcript_id", values="weight"
    ).fillna(0)

    V_raw = matrix_df.values.astype(np.float64)
    mirna_ids = matrix_df.index.tolist()
    transcript_ids = matrix_df.columns.tolist()

    print(f"  Matrix shape: {V_raw.shape[0]} miRNAs × {V_raw.shape[1]} transcripts")
    print(f"  Value range: [{V_raw.min():.3f}, {V_raw.max():.3f}]")

    col_means = V_raw.mean(axis=0)
    V_residual = np.clip(V_raw - col_means, 0, None)
    V_residual[V_residual == 0] = 1e-6

    return V_raw, V_residual, mirna_ids, transcript_ids


# ────────────────────────────────────────────────────────────────
# Edge extraction and filtering
# ────────────────────────────────────────────────────────────────
def extract_edges(
    matrix: np.ndarray,
    mirna_ids: list[str],
    transcript_ids: list[str],
    threshold: float,
    top_edges: int,
    top_mirnas: int,
    min_degree: int,
    power: float,
) -> pd.DataFrame:
    """Extract edges from binding matrix and apply filters.

    Returns a DataFrame with columns: mirna, gene, weight, pathway.
    """
    gene_names = [TRANSCRIPT_GENE_MAP.get(tid, tid) for tid in transcript_ids]

    rows: list[dict] = []
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            w = matrix[i, j]
            if w >= threshold:
                rows.append({
                    "mirna": mirna_ids[i],
                    "gene": gene_names[j],
                    "weight": w,
                })

    if not rows:
        print("  WARNING: No edges above threshold. Try lowering --threshold.")
        return pd.DataFrame(columns=["mirna", "gene", "weight", "pathway"])

    edges = pd.DataFrame(rows)
    edges["weight"] = np.power(edges["weight"].values, power)
    edges["pathway"] = edges["gene"].apply(gene_to_pathway)

    print(f"  Edges after threshold (≥{threshold}): {len(edges):,}")

    if top_edges > 0 and len(edges) > top_edges:
        edges = edges.nlargest(top_edges, "weight")
        print(f"  Edges after top-edges filter (top {top_edges}): {len(edges):,}")

    # Filter by min_degree — iterate until stable
    if min_degree > 1:
        for _ in range(5):
            mirna_counts = edges["mirna"].value_counts()
            gene_counts = edges["gene"].value_counts()
            valid_mirnas = mirna_counts[mirna_counts >= min_degree].index
            valid_genes = gene_counts[gene_counts >= min_degree].index
            before = len(edges)
            edges = edges[edges["mirna"].isin(valid_mirnas) & edges["gene"].isin(valid_genes)]
            if len(edges) == before:
                break
        print(f"  Edges after min-degree filter (≥{min_degree}): {len(edges):,}")

    if top_mirnas > 0:
        mirna_total = edges.groupby("mirna")["weight"].sum()
        top_m = mirna_total.nlargest(top_mirnas).index
        edges = edges[edges["mirna"].isin(top_m)]
        print(f"  Edges after top-mirnas filter (top {top_mirnas}): {len(edges):,}")

    return edges.reset_index(drop=True)


# ────────────────────────────────────────────────────────────────
# Graph rendering
# ────────────────────────────────────────────────────────────────
def render_bipartite(
    edges: pd.DataFrame,
    output: str,
    edge_scale: float = 3.0,
    dpi: int = 200,
    width: float = 20,
    height: float = 14,
    label_fontsize: int = 7,
    gene_fontsize: int = 9,
    threshold: float = 15.0,
    use_residual: bool = False,
    power: float = 1.0,
) -> None:
    """Render bipartite graph to PNG using matplotlib only."""
    if edges.empty:
        print("  No edges to plot.")
        return

    unique_mirnas = sorted(edges["mirna"].unique())
    unique_genes = sorted(edges["gene"].unique(),
                          key=lambda g: (gene_to_pathway(g), g))

    n_mirnas = len(unique_mirnas)
    n_genes = len(unique_genes)

    mirna_y = {m: i for i, m in enumerate(unique_mirnas)}
    gene_y = {g: i for i, g in enumerate(unique_genes)}

    n_max = max(n_mirnas, n_genes)
    mirna_scale = n_max / n_mirnas if n_mirnas > 0 else 1.0
    gene_scale = n_max / n_genes if n_genes > 0 else 1.0

    fig, ax = plt.subplots(figsize=(width, height))

    x_left = 0.0
    x_right = 1.0

    w_min = edges["weight"].min()
    w_max = edges["weight"].max()
    w_range = w_max - w_min if w_max > w_min else 1.0

    edges_sorted = edges.sort_values("weight", ascending=True)

    for _, row in edges_sorted.iterrows():
        y_left = mirna_y[row["mirna"]] * mirna_scale
        y_right = gene_y[row["gene"]] * gene_scale
        norm_w = (row["weight"] - w_min) / w_range
        lw = 0.3 + norm_w * edge_scale
        alpha = 0.08 + norm_w * 0.72

        pathway = row["pathway"]
        color = PATHWAY_COLORS.get(pathway, "#95a5a6")

        mid_x = 0.5
        ctrl_x1 = x_left + 0.3
        ctrl_x2 = x_right - 0.3
        xs = [x_left, ctrl_x1, ctrl_x2, x_right]
        ys = [y_left, y_left, y_right, y_right]

        from matplotlib.path import Path
        import matplotlib.patches as mpatches
        verts = [(xs[0], ys[0]),
                 (xs[1], ys[1]),
                 (xs[2], ys[2]),
                 (xs[3], ys[3])]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        path = Path(verts, codes)
        patch = mpatches.PathPatch(
            path,
            facecolor="none",
            edgecolor=color,
            alpha=alpha,
            linewidth=lw,
        )
        ax.add_patch(patch)

    for m in unique_mirnas:
        y = mirna_y[m] * mirna_scale
        ax.plot(x_left, y, "o", color="#7f8c8d", markersize=3, zorder=5)
        ax.text(x_left - 0.01, y, m, ha="right", va="center",
                fontsize=label_fontsize, color="#2c3e50", family="monospace")

    for g in unique_genes:
        y = gene_y[g] * gene_scale
        pathway = gene_to_pathway(g)
        color = PATHWAY_COLORS[pathway]
        ax.plot(x_right, y, "o", color=color, markersize=7, zorder=5)
        ax.text(x_right + 0.01, y, g, ha="left", va="center",
                fontsize=gene_fontsize, fontweight="bold", color="#2c3e50")

    ax.text(x_left, n_max + 0.5, "miRNAs", ha="center", va="bottom",
            fontsize=11, fontweight="bold", color="#2c3e50")
    ax.text(x_right, n_max + 0.5, "Target Genes", ha="center", va="bottom",
            fontsize=11, fontweight="bold", color="#2c3e50")

    legend_handles = []
    for pathway, color in PATHWAY_COLORS.items():
        count = len(unique_genes)
        legend_handles.append(
            plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color,
                       markersize=8, label=f"{pathway}")
        )
    ax.legend(handles=legend_handles, loc="upper center", ncol=4,
              fontsize=9, frameon=True, fancybox=True,
              bbox_to_anchor=(0.5, -0.02))

    matrix_label = "residual" if use_residual else "raw"
    title = (
        f"miRNA–Gene Binding Bipartite Graph "
        f"({matrix_label} matrix, threshold={threshold}, "
        f"power={power})\n"
        f"{n_mirnas} miRNAs × {n_genes} genes — {len(edges)} edges"
    )
    ax.set_title(title, fontsize=13, fontweight="bold", pad=15)

    ax.set_xlim(-0.35, 1.35)
    ax.set_ylim(-1, n_max + 1.5)
    ax.axis("off")

    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
    fig.savefig(output, dpi=dpi, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"  ✓ Saved: {output}")


# ────────────────────────────────────────────────────────────────
# Summary table
# ────────────────────────────────────────────────────────────────
def print_summary(edges: pd.DataFrame) -> None:
    """Print graph statistics and top-10 edges table."""
    if edges.empty:
        print("  No edges to summarize.")
        return

    unique_mirnas = edges["mirna"].nunique()
    unique_genes = edges["gene"].nunique()
    n_edges = len(edges)
    max_possible = unique_mirnas * unique_genes
    density = n_edges / max_possible if max_possible > 0 else 0.0

    print(f"\n  Graph summary:")
    print(f"    miRNA nodes:   {unique_mirnas}")
    print(f"    Gene nodes:    {unique_genes}")
    print(f"    Edges:         {n_edges}")
    print(f"    Density:       {density:.3%}")
    print(f"    Weight range:  [{edges['weight'].min():.3f}, {edges['weight'].max():.3f}]")

    mirna_deg = edges.groupby("mirna")["gene"].nunique()
    gene_deg = edges.groupby("gene")["mirna"].nunique()
    print(f"    miRNA degree:  mean={mirna_deg.mean():.1f}, max={mirna_deg.max()}")
    print(f"    Gene degree:   mean={gene_deg.mean():.1f}, max={gene_deg.max()}")

    print(f"\n  Top 10 edges by weight:")
    print(f"  {'miRNA':<22s} {'Gene':<10s} {'Weight':>8s}  {'Pathway':<20s}")
    print(f"  {'─'*22} {'─'*10} {'─'*8}  {'─'*20}")
    top10 = edges.nlargest(10, "weight")
    for _, row in top10.iterrows():
        print(f"  {row['mirna']:<22s} {row['gene']:<10s} "
              f"{row['weight']:>8.2f}  {row['pathway']:<20s}")

    # Per-gene summary
    print(f"\n  Gene connectivity (sorted by edge count):")
    print(f"  {'Gene':<10s} {'Edges':>6s} {'Max weight':>11s} {'Pathway':<20s}")
    print(f"  {'─'*10} {'─'*6} {'─'*11} {'─'*20}")
    gene_stats = (edges.groupby(["gene", "pathway"])
                  .agg(n_edges=("weight", "count"),
                       max_weight=("weight", "max"))
                  .reset_index()
                  .sort_values("n_edges", ascending=False))
    for _, row in gene_stats.iterrows():
        print(f"  {row['gene']:<10s} {row['n_edges']:>6d} "
              f"{row['max_weight']:>11.2f} {row['pathway']:<20s}")


# ────────────────────────────────────────────────────────────────
# CLI
# ────────────────────────────────────────────────────────────────
def parse_args() -> argparse.Namespace:
    """Parse CLI args for the bipartite graph renderer.

    Returns:
        The parsed ``argparse.Namespace`` (threshold, top-mirnas, top-edges,
        output path, residual flag, power weighting, etc.).
    """
    parser = argparse.ArgumentParser(
        description="Bipartite graph visualization for miRNA–transcript binding"
    )
    parser.add_argument("--threshold", type=float, default=15.0,
                        help="Minimum binding strength to show an edge")
    parser.add_argument("--top-mirnas", type=int, default=50,
                        help="Show only top N miRNAs by total binding weight")
    parser.add_argument("--top-edges", type=int, default=200,
                        help="Show only top N edges by weight")
    parser.add_argument("--output", type=str, default="figures/bipartite_graph.png",
                        help="Output file path")
    parser.add_argument("--use-residual", action="store_true",
                        help="Use residual matrix (col-mean subtracted, clipped)")
    parser.add_argument("--power", type=float, default=1.0,
                        help="Power weighting for binding values")
    parser.add_argument("--exclude-seedless", action="store_true", default=True,
                        help="Exclude seedless interactions (default: True)")
    parser.add_argument("--min-degree", type=int, default=1,
                        help="Minimum edges a node must have to be shown")
    parser.add_argument("--edge-scale", type=float, default=3.0,
                        help="Scaling factor for edge width")
    parser.add_argument("--dpi", type=int, default=200,
                        help="Output DPI")
    parser.add_argument("--width", type=float, default=20,
                        help="Figure width in inches")
    parser.add_argument("--height", type=float, default=14,
                        help="Figure height in inches")
    parser.add_argument("--label-fontsize", type=int, default=7,
                        help="Font size for miRNA labels")
    parser.add_argument("--gene-fontsize", type=int, default=9,
                        help="Font size for gene labels")
    return parser.parse_args()


# ────────────────────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────────────────────
def main() -> None:
    """CLI entry: load matrix, filter, draw bipartite graph (raw or residual), save PNG."""
    args = parse_args()

    separator = "=" * 60
    print(f"\n  {separator}")
    print(f"  Bipartite Graph: miRNA–Transcript Binding")
    print(f"  {separator}")

    # 1. Load data
    print(f"\n  1. Loading binding matrix from DuckDB")
    V_raw, V_residual, mirna_ids, transcript_ids = load_binding_matrix(
        exclude_seedless=args.exclude_seedless,
    )

    # 2. Select matrix
    if args.use_residual:
        matrix = V_residual
        print(f"  Using residual matrix")
    else:
        matrix = V_raw
        print(f"  Using raw matrix")

    # 3. Extract and filter edges
    print(f"\n  2. Filtering edges")
    edges = extract_edges(
        matrix=matrix,
        mirna_ids=mirna_ids,
        transcript_ids=transcript_ids,
        threshold=args.threshold,
        top_edges=args.top_edges,
        top_mirnas=args.top_mirnas,
        min_degree=args.min_degree,
        power=args.power,
    )

    # 4. Render graph
    print(f"\n  3. Rendering bipartite graph")
    render_bipartite(
        edges=edges,
        output=args.output,
        edge_scale=args.edge_scale,
        dpi=args.dpi,
        width=args.width,
        height=args.height,
        label_fontsize=args.label_fontsize,
        gene_fontsize=args.gene_fontsize,
        threshold=args.threshold,
        use_residual=args.use_residual,
        power=args.power,
    )

    # 5. Summary
    print(f"\n  4. Summary")
    print_summary(edges)

    print(f"\n  {separator}")
    print(f"  Done.")
    print(f"  {separator}\n")


if __name__ == "__main__":
    main()

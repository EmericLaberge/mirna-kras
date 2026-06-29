#!/usr/bin/env python3
"""
Bipartite network visualization filtered by NMF cluster.

For each NMF signature (k=2), select the miRNAs where that signature is
dominant, compute residual binding, and draw a bipartite graph showing the
strongest differential bindings.

Usage:
    .venv/bin/python3 scripts/bipartite_by_cluster.py \\
        --top-genes 3 --threshold 0.5 --max-mirnas 100 \\
        --output figures/bipartite_by_cluster.png
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
import matplotlib.patches as mpatches
from matplotlib.path import Path


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

KRAS_PATHWAY = {
    "KRAS", "NRAS", "HRAS", "PIK3CA", "PIK3CB", "MTOR",
    "BRAF", "RAF1", "MAP2K1", "MAP2K2", "MAPK1",
}
TUMOR_SUPPRESSORS = {"TP53", "PTEN", "RB1"}
MAPK_PATHWAY = {"EGFR", "ERBB2"}

PATHWAY_COLORS = {
    "KRas": "#e74c3c",
    "Tumor suppressor": "#3498db",
    "MAPK": "#2ecc71",
    "Other": "#95a5a6",
}

SIGNATURE_LABELS = {
    0: "Signature 1 (RAS–MAPK)",
    1: "Signature 2 (PI3K–AKT)",
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
def load_data(
    db_path: str = "rimap_results.duckdb",
    w_path: str = "figures_nmf/W_matrix.npy",
    h_path: str = "figures_nmf/H_matrix.npy",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[str], list[str]]:
    """Load binding matrix + NMF factors.

    Returns (V_raw, V_residual, W, H, mirna_ids, transcript_ids).
    """
    con = duckdb.connect(db_path, read_only=True)
    df = con.execute("SELECT * FROM results").fetchdf()
    con.close()

    df = df[df["seed_type"] != "seedless"]
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

    # Residual: subtract column means, clip to >= 0
    col_means = V_raw.mean(axis=0)
    V_residual = np.clip(V_raw - col_means, 0, None)

    W = np.load(w_path)
    H = np.load(h_path)

    print(f"  Matrix: {V_raw.shape[0]} miRNAs × {V_raw.shape[1]} transcripts")
    print(f"  W: {W.shape}  H: {H.shape}")
    return V_raw, V_residual, W, H, mirna_ids, transcript_ids


# ────────────────────────────────────────────────────────────────
# Edge extraction for one signature
# ────────────────────────────────────────────────────────────────
def extract_edges_for_signature(
    V_residual: np.ndarray,
    W: np.ndarray,
    mirna_ids: list[str],
    transcript_ids: list[str],
    sig_idx: int,
    top_genes: int = 3,
    threshold: float = 0.5,
    max_mirnas: int = 100,
) -> pd.DataFrame:
    """Build edge list for one NMF signature's dominant miRNAs.

    Steps:
      1. Select miRNAs whose dominant signature == sig_idx.
      2. For each miRNA, keep top K genes by residual binding.
      3. Filter edges below threshold.
      4. Keep only genes that have at least one edge above threshold.
      5. Limit to top max_mirnas by total residual binding.
    """
    gene_names = [TRANSCRIPT_GENE_MAP.get(tid, tid) for tid in transcript_ids]
    dominant = W.argmax(axis=1)
    mask = dominant == sig_idx
    mirna_indices = np.where(mask)[0]

    if len(mirna_indices) == 0:
        return pd.DataFrame(columns=["mirna", "gene", "weight", "pathway"])

    # Sub-select to top max_mirnas by total residual binding
    total_binding = V_residual[mirna_indices].sum(axis=1)
    top_order = np.argsort(total_binding)[::-1][:max_mirnas]
    mirna_indices = mirna_indices[top_order]

    # For each miRNA, keep top K genes
    rows = []
    for mi in mirna_indices:
        vals = V_residual[mi]
        top_k_idx = np.argsort(vals)[::-1][:top_genes]
        for j in top_k_idx:
            w = vals[j]
            if w >= threshold:
                rows.append({
                    "mirna": mirna_ids[mi],
                    "gene": gene_names[j],
                    "weight": float(w),
                })

    if not rows:
        return pd.DataFrame(columns=["mirna", "gene", "weight", "pathway"])

    edges = pd.DataFrame(rows)
    edges["pathway"] = edges["gene"].apply(gene_to_pathway)

    # Keep only genes with at least one edge above threshold (already enforced)
    # Re-filter: genes must have at least one edge
    gene_counts = edges["gene"].value_counts()
    valid_genes = gene_counts[gene_counts >= 1].index
    edges = edges[edges["gene"].isin(valid_genes)]

    return edges.reset_index(drop=True)


# ────────────────────────────────────────────────────────────────
# Rendering one subplot
# ────────────────────────────────────────────────────────────────
def render_subplot(
    ax: plt.Axes,
    edges: pd.DataFrame,
    sig_label: str,
    edge_scale: float = 4.0,
    label_top: int = 10,
) -> None:
    """Draw a bipartite graph on the given Axes."""
    if edges.empty:
        ax.text(0.5, 0.5, "No edges", transform=ax.transAxes,
                ha="center", va="center", fontsize=14, color="gray")
        ax.set_title(sig_label, fontsize=13, fontweight="bold", pad=12)
        ax.axis("off")
        return

    unique_mirnas = edges.groupby("mirna")["weight"].sum().sort_values(ascending=False).index.tolist()
    unique_genes = sorted(edges["gene"].unique(),
                          key=lambda g: (gene_to_pathway(g), g))

    n_mirnas = len(unique_mirnas)
    n_genes = len(unique_genes)

    mirna_y = {m: i for i, m in enumerate(unique_mirnas)}
    gene_y = {g: i for i, g in enumerate(unique_genes)}

    n_max = max(n_mirnas, n_genes)
    mirna_scale = n_max / n_mirnas if n_mirnas > 0 else 1.0
    gene_scale = n_max / n_genes if n_genes > 0 else 1.0

    x_left = 0.0
    x_right = 1.0

    w_min = edges["weight"].min()
    w_max = edges["weight"].max()
    w_range = w_max - w_min if w_max > w_min else 1.0

    # Draw edges (sorted so stronger on top)
    edges_sorted = edges.sort_values("weight", ascending=True)
    for _, row in edges_sorted.iterrows():
        y_left = mirna_y[row["mirna"]] * mirna_scale
        y_right = gene_y[row["gene"]] * gene_scale
        norm_w = (row["weight"] - w_min) / w_range
        lw = 0.3 + norm_w * edge_scale
        alpha = 0.06 + norm_w * 0.65

        pathway = row["pathway"]
        color = PATHWAY_COLORS.get(pathway, "#95a5a6")

        # Cubic Bezier
        verts = [
            (x_left, y_left),
            (x_left + 0.35, y_left),
            (x_right - 0.35, y_right),
            (x_right, y_right),
        ]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        path = Path(verts, codes)
        patch = mpatches.PathPatch(
            path, facecolor="none", edgecolor=color,
            alpha=alpha, linewidth=lw,
        )
        ax.add_patch(patch)

    # Draw miRNA nodes
    for m in unique_mirnas:
        y = mirna_y[m] * mirna_scale
        ax.plot(x_left, y, "o", color="#7f8c8d", markersize=2, zorder=5)

    # Label top N miRNAs
    top_mirnas = unique_mirnas[:label_top]
    for m in top_mirnas:
        y = mirna_y[m] * mirna_scale
        ax.text(
            x_left - 0.02, y, m, ha="right", va="center",
            fontsize=4.5, color="#2c3e50", family="monospace",
        )

    # Draw gene nodes
    for g in unique_genes:
        y = gene_y[g] * gene_scale
        pathway = gene_to_pathway(g)
        color = PATHWAY_COLORS[pathway]
        ax.plot(x_right, y, "o", color=color, markersize=8, zorder=5)
        ax.text(
            x_right + 0.02, y, g, ha="left", va="center",
            fontsize=7, fontweight="bold", color="#2c3e50",
        )

    # Column headers
    ax.text(x_left, n_max + 0.8, "miRNAs", ha="center", va="bottom",
            fontsize=10, fontweight="bold", color="#2c3e50")
    ax.text(x_right, n_max + 0.8, "Genes", ha="center", va="bottom",
            fontsize=10, fontweight="bold", color="#2c3e50")

    ax.set_title(
        f"{sig_label}\n{n_mirnas} miRNAs × {n_genes} genes — {len(edges)} edges",
        fontsize=12, fontweight="bold", pad=12,
    )
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylim(-1, n_max + 2)
    ax.axis("off")


# ────────────────────────────────────────────────────────────────
# CLI
# ────────────────────────────────────────────────────────────────
def parse_args() -> argparse.Namespace:
    """Parse CLI args: --top-genes, --threshold, --max-mirnas, --output, --dpi, --edge-scale, --label-top.

    Returns:
        The parsed ``argparse.Namespace``.
    """
    parser = argparse.ArgumentParser(
        description="Bipartite network visualization filtered by NMF cluster"
    )
    parser.add_argument("--top-genes", type=int, default=3,
                        help="Top K genes per miRNA to keep")
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="Minimum residual binding to show edge")
    parser.add_argument("--max-mirnas", type=int, default=100,
                        help="Max miRNAs to show (highest total residual binding)")
    parser.add_argument("--output", type=str, default="figures/bipartite_by_cluster.png",
                        help="Output file path")
    parser.add_argument("--dpi", type=int, default=200,
                        help="Output DPI")
    parser.add_argument("--edge-scale", type=float, default=4.0,
                        help="Scaling factor for edge width")
    parser.add_argument("--label-top", type=int, default=10,
                        help="Number of top miRNAs to label")
    return parser.parse_args()


# ────────────────────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────────────────────
def main() -> None:
    """CLI entry: load V + W + H, draw per-signature bipartite panels, save PNG."""
    args = parse_args()

    sep = "=" * 60
    print(f"\n  {sep}")
    print(f"  Bipartite Network by NMF Cluster")
    print(f"  {sep}")

    # 1. Load data
    print(f"\n  1. Loading data")
    V_raw, V_residual, W, H, mirna_ids, transcript_ids = load_data()

    dominant = W.argmax(axis=1)
    for i in range(W.shape[1]):
        n = (dominant == i).sum()
        print(f"    Signature {i+1}: {n} miRNAs")

    # 2. Extract edges per signature
    print(f"\n  2. Extracting edges (top-genes={args.top_genes}, "
          f"threshold={args.threshold}, max-mirnas={args.max_mirnas})")

    all_edges = {}
    for sig_idx in range(W.shape[1]):
        edges = extract_edges_for_signature(
            V_residual, W, mirna_ids, transcript_ids,
            sig_idx=sig_idx,
            top_genes=args.top_genes,
            threshold=args.threshold,
            max_mirnas=args.max_mirnas,
        )
        n_mirnas = edges["mirna"].nunique() if not edges.empty else 0
        n_genes = edges["gene"].nunique() if not edges.empty else 0
        print(f"    Sig {sig_idx+1}: {len(edges)} edges, "
              f"{n_mirnas} miRNAs, {n_genes} genes")
        all_edges[sig_idx] = edges

    # 3. Render
    print(f"\n  3. Rendering figure")
    fig, axes = plt.subplots(1, 2, figsize=(28, 16))

    for sig_idx in range(W.shape[1]):
        label = SIGNATURE_LABELS.get(sig_idx, f"Signature {sig_idx+1}")
        render_subplot(
            axes[sig_idx],
            all_edges[sig_idx],
            sig_label=label,
            edge_scale=args.edge_scale,
            label_top=args.label_top,
        )

    # Legend (shared)
    legend_handles = []
    for pathway, color in PATHWAY_COLORS.items():
        legend_handles.append(
            plt.Line2D([0], [0], marker="o", color="w",
                       markerfacecolor=color, markersize=8, label=pathway)
        )
    fig.legend(
        handles=legend_handles, loc="lower center", ncol=4,
        fontsize=10, frameon=True, fancybox=True, bbox_to_anchor=(0.5, 0.01),
    )

    fig.suptitle(
        "miRNA–Gene Bipartite Networks by NMF Signature",
        fontsize=16, fontweight="bold", y=0.98,
    )
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"  ✓ Saved: {args.output}")

    # 4. Summary
    print(f"\n  4. Summary")
    for sig_idx in range(W.shape[1]):
        edges = all_edges[sig_idx]
        label = SIGNATURE_LABELS.get(sig_idx, f"Signature {sig_idx+1}")
        if edges.empty:
            print(f"\n    {label}: no edges")
            continue
        print(f"\n    {label}:")
        print(f"      miRNAs: {edges['mirna'].nunique()}")
        print(f"      Genes:  {edges['gene'].nunique()}")
        print(f"      Edges:  {len(edges)}")
        gene_stats = edges.groupby("gene").agg(
            n_edges=("weight", "count"),
            max_w=("weight", "max"),
            mean_w=("weight", "mean"),
        ).sort_values("max_w", ascending=False)
        print(f"      Gene connectivity:")
        for gene, row in gene_stats.iterrows():
            pw = gene_to_pathway(gene)
            print(f"        {gene:<10s} edges={int(row['n_edges']):>4d}  "
                  f"max={row['max_w']:.2f}  mean={row['mean_w']:.2f}  [{pw}]")

    print(f"\n  {sep}")
    print(f"  Done.")
    print(f"  {sep}\n")


if __name__ == "__main__":
    main()

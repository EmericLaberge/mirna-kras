#!/usr/bin/env python3
"""
Top miRNAs per Gene Heatmap.

Selects the top-N miRNAs (by residual binding strength) for each of the 25
genes, takes the union of those miRNAs, and plots a clustered heatmap of
binding strength across all genes.

Usage:
    .venv/bin/python3 scripts/bipartite_top_heatmap.py
    .venv/bin/python3 scripts/bipartite_top_heatmap.py --top-n 5 --use-residual
    .venv/bin/python3 scripts/bipartite_top_heatmap.py --top-n 3 --output figures/custom.png
"""

import warnings
warnings.filterwarnings("ignore")

import os
import argparse
import numpy as np
import pandas as pd
import duckdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram

TRANSCRIPT_GENE_MAP = {
    "ENST00000215832": "RAF1", "ENST00000226574": "TP53", "ENST00000227507": "ERBB2",
    "ENST00000251849": "PTEN", "ENST00000256078": "BRAF", "ENST00000259523": "MAP2K1",
    "ENST00000262948": "KRAS", "ENST00000263025": "NRAS", "ENST00000263826": "HRAS",
    "ENST00000263967": "MAP2K2", "ENST00000267163": "MAPK3", "ENST00000269305": "MAPK1",
    "ENST00000275493": "RB1", "ENST00000307102": "DUSP6", "ENST00000308639": "MTOR",
    "ENST00000311278": "PIK3CB", "ENST00000315204": "APC", "ENST00000318493": "NF1",
    "ENST00000316629": "PIK3CA", "ENST00000331710": "AKT1", "ENST00000349310": "AKT2",
    "ENST00000361445": "EGFR", "ENST00000366790": "MYC", "ENST00000371953": "SPRY2",
    "ENST00000377687": "CDKN1A",
}

# Pathway groups for column ordering
PATHWAY_GROUPS = {
    "RAS": ["KRAS", "NRAS", "HRAS"],
    "RAF-MEK-ERK": ["BRAF", "RAF1", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "DUSP6", "SPRY2"],
    "PI3K-AKT-mTOR": ["PIK3CA", "PIK3CB", "AKT1", "AKT2", "MTOR"],
    "RTK": ["EGFR", "ERBB2"],
    "Tumor Suppressors": ["TP53", "PTEN", "RB1", "NF1", "APC", "CDKN1A"],
    "Oncogenes": ["MYC"],
}

PATHWAY_COLORS = {
    "RAS": "#e74c3c",
    "RAF-MEK-ERK": "#e67e22",
    "PI3K-AKT-mTOR": "#2ecc71",
    "RTK": "#3498db",
    "Tumor Suppressors": "#9b59b6",
    "Oncogenes": "#1abc9c",
}


def load_data(db_path: str = "rimap_results.duckdb"):
    """Load the miRNA-transcript binding matrix from DuckDB."""
    con = duckdb.connect(db_path, read_only=True)
    df = con.execute("SELECT * FROM results").fetchdf()
    con.close()

    n_total = len(df)
    print(f"  Total interactions: {n_total:,}")

    df = df[df["seed_type"] != "seedless"]
    n_filtered = len(df)
    print(f"  After removing seedless: {n_filtered:,}")

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

    return V_raw, mirna_ids, transcript_ids


def compute_residual(V_raw: np.ndarray) -> np.ndarray:
    """Compute residual matrix: clip(V_raw - col_means, 0, None), fill zeros with 1e-6."""
    col_means = V_raw.mean(axis=0)
    V_res = np.clip(V_raw - col_means, 0, None)
    V_res[V_res == 0] = 1e-6
    return V_res


def get_pathway_group(gene: str) -> str:
    """Find which pathway group a gene belongs to."""
    for group, genes in PATHWAY_GROUPS.items():
        if gene in genes:
            return group
    return "Other"


def order_columns_by_pathway(gene_names: list[str]) -> list[int]:
    """Return column indices sorted by pathway group, then by gene name within group."""
    group_order = list(PATHWAY_GROUPS.keys())
    gene_group = [(g, get_pathway_group(g), group_order.index(get_pathway_group(g)))
                  for g in gene_names]
    gene_group.sort(key=lambda x: (x[2], x[0]))
    # Return indices into original gene_names
    indexed = list(enumerate(gene_names))
    indexed.sort(key=lambda x: (get_pathway_group(x[1]), group_order.index(get_pathway_group(x[1])), x[1]))
    return [i for i, _ in indexed]


def parse_args() -> argparse.Namespace:
    """Parse CLI args: --top-n, --use-residual/--no-residual, --output, --dpi, --figsize.

    Returns:
        The parsed ``argparse.Namespace``.
    """
    parser = argparse.ArgumentParser(
        description="Heatmap of top miRNAs per gene by binding strength"
    )
    parser.add_argument("--top-n", type=int, default=5,
                        help="Top N miRNAs per gene to include")
    parser.add_argument("--use-residual", action="store_true", default=True,
                        help="Use residual matrix (default: True)")
    parser.add_argument("--no-residual", dest="use_residual", action="store_false",
                        help="Use raw matrix instead of residual")
    parser.add_argument("--output", type=str, default="figures/top_mirnas_heatmap.png",
                        help="Output file path")
    parser.add_argument("--dpi", type=int, default=200,
                        help="DPI for output figure")
    parser.add_argument("--figsize", type=str, default="16,14",
                        help="Figure size as width,height (e.g. 16,14)")
    return parser.parse_args()


def main():
    """CLI entry: load matrix -> top-N miRNAs per gene -> clustered heatmap."""
    args = parse_args()
    figsize = tuple(int(x) for x in args.figsize.split(","))

    plt.style.use("seaborn-v0_8-whitegrid")

    # ── 1. Load data ───────────────────────────────────────────
    print("Loading data...")
    V_raw, mirna_ids, transcript_ids = load_data()
    gene_names = [TRANSCRIPT_GENE_MAP.get(tid, tid) for tid in transcript_ids]

    # ── 2. Compute matrix ──────────────────────────────────────
    if args.use_residual:
        print("Computing residual matrix...")
        V = compute_residual(V_raw)
    else:
        V = V_raw.copy()

    # ── 3. Find top-N miRNAs per gene ──────────────────────────
    print(f"Selecting top-{args.top_n} miRNAs per gene...")
    top_mirna_indices = set()
    for j in range(V.shape[1]):
        col = V[:, j]
        top_idx = np.argsort(col)[::-1][:args.top_n]
        top_mirna_indices.update(top_idx.tolist())

    top_mirna_indices = sorted(top_mirna_indices)
    n_unique = len(top_mirna_indices)
    print(f"  Unique top miRNAs across all genes: {n_unique}")

    # ── 4. Build sub-matrix ────────────────────────────────────
    sub_matrix = V[np.ix_(top_mirna_indices, range(V.shape[1]))]
    sub_mirna_ids = [mirna_ids[i] for i in top_mirna_indices]

    # ── 5. Sort columns by pathway group ───────────────────────
    col_order = order_columns_by_pathway(gene_names)
    sub_matrix = sub_matrix[:, col_order]
    sorted_gene_names = [gene_names[i] for i in col_order]
    sorted_pathways = [get_pathway_group(g) for g in sorted_gene_names]

    # ── 6. Cluster rows (miRNAs) ───────────────────────────────
    print("Clustering miRNAs...")
    # Use only non-constant rows for clustering
    row_vars = sub_matrix.var(axis=1)
    if row_vars.max() > 0:
        Z = linkage(sub_matrix, method="average", metric="euclidean")
        row_order = dendrogram(Z, no_plot=True)["leaves"]
    else:
        row_order = list(range(sub_matrix.shape[0]))

    sub_matrix = sub_matrix[row_order, :]
    sub_mirna_ids = [sub_mirna_ids[i] for i in row_order]

    # ── 7. Build DataFrame for seaborn ─────────────────────────
    df_heat = pd.DataFrame(
        sub_matrix,
        index=sub_mirna_ids,
        columns=sorted_gene_names,
    )

    # ── 8. Plot ────────────────────────────────────────────────
    print("Plotting heatmap...")

    n_rows = len(sub_mirna_ids)
    n_cols = len(sorted_gene_names)

    # Dynamic row label size based on number of miRNAs
    label_fontsize = max(5, min(9, 300 // n_rows))

    fig, ax = plt.subplots(figsize=figsize)

    cmap = "YlOrRd"
    vmax = float(np.percentile(sub_matrix[sub_matrix > 1e-6], 95))

    sns.heatmap(
        df_heat,
        cmap=cmap,
        linewidths=0.3,
        linecolor="white",
        vmin=0,
        vmax=vmax,
        ax=ax,
        cbar_kws={
            "label": "Residual binding strength" if args.use_residual else "Mean |binding_dG|",
            "shrink": 0.5,
        },
        xticklabels=True,
        yticklabels=True,
    )

    ax.set_yticklabels(ax.get_yticklabels(), fontsize=label_fontsize, fontfamily="monospace")
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=9, rotation=45, ha="right")

    title = f"Top {args.top_n} miRNAs per Gene — Binding Strength Heatmap"
    title += f" ({n_unique} unique miRNAs)"
    ax.set_title(title, fontsize=13, pad=40)
    ax.set_xlabel("Gene", fontsize=11)
    ax.set_ylabel("miRNA", fontsize=11)

    # ── Column color bar for pathway groups ────────────────────
    # Draw colored bars at the top of the heatmap
    group_boundaries = []
    current_group = sorted_pathways[0]
    start = 0
    for i, grp in enumerate(sorted_pathways + ["END"]):
        if grp != current_group:
            group_boundaries.append((current_group, start, i - 1))
            current_group = grp
            start = i

    # Draw colored patches above the heatmap
    y_top = ax.get_ylim()[1]
    patch_height = (y_top / n_rows) * 1.2
    for grp, s, e in group_boundaries:
        color = PATHWAY_COLORS.get(grp, "#95a5a6")
        rect = mpatches.FancyBboxPatch(
            (s, y_top + patch_height * 0.3),
            e - s + 1,
            patch_height,
            boxstyle="square,pad=0",
            facecolor=color,
            edgecolor="white",
            linewidth=0.5,
            clip_on=False,
            transform=ax.transData,
        )
        ax.add_patch(rect)
        # Label
        mid = (s + e) / 2.0
        ax.text(
            mid, y_top + patch_height * 1.5, grp,
            ha="center", va="bottom", fontsize=7,
            fontweight="bold", color=color,
            clip_on=False,
        )

    # ── Legend for pathway colors ──────────────────────────────
    legend_patches = [
        mpatches.Patch(facecolor=PATHWAY_COLORS[g], edgecolor="gray", label=g)
        for g in PATHWAY_GROUPS if g in set(sorted_pathways)
    ]
    ax.legend(
        handles=legend_patches,
        loc="upper left",
        bbox_to_anchor=(1.15, 1.0),
        fontsize=8,
        frameon=True,
        title="Pathway",
        title_fontsize=9,
    )

    fig.tight_layout()

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    print(f"\n  ✓ Saved: {args.output}")
    print(f"    miRNAs: {n_unique}, Genes: {n_cols}")
    print(f"    Value range: [{sub_matrix[sub_matrix > 1e-6].min():.4f}, {sub_matrix.max():.4f}]")

    # ── Print top miRNAs per gene ──────────────────────────────
    print(f"\n  Top {args.top_n} miRNAs per gene:")
    for j, gene in enumerate(sorted_gene_names):
        col = sub_matrix[:, j]
        top_idx = np.argsort(col)[::-1][:args.top_n]
        names = [sub_mirna_ids[i] for i in top_idx]
        vals = [col[i] for i in top_idx]
        print(f"    {gene:>10s}: {', '.join(f'{n} ({v:.3f})' for n, v in zip(names, vals))}")


if __name__ == "__main__":
    main()

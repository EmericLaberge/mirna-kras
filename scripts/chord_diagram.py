#!/usr/bin/env python3
"""
Chord Diagram: Gene-Gene Co-targeting by miRNAs.

Shows which gene pairs are co-targeted by the same miRNAs,
based on shared above-threshold residual binding.

Usage:
    .venv/bin/python3 scripts/chord_diagram.py
    .venv/bin/python3 scripts/chord_diagram.py --threshold 150 --dpi 300
"""

import warnings
warnings.filterwarnings("ignore")

import argparse
import os
import numpy as np
import duckdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
import matplotlib.patheffects as pe

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

PATHWAY_COLORS = {
    # RAS family — red
    "KRAS": "#e41a1c", "NRAS": "#e41a1c", "HRAS": "#e41a1c",
    # Tumor suppressors — blue
    "TP53": "#377eb8", "PTEN": "#377eb8", "RB1": "#377eb8",
    "APC": "#377eb8", "NF1": "#377eb8", "CDKN1A": "#377eb8",
    # MAPK cascade — green
    "RAF1": "#4daf4a", "BRAF": "#4daf4a", "MAP2K1": "#4daf4a",
    "MAP2K2": "#4daf4a", "MAPK3": "#4daf4a", "MAPK1": "#4daf4a",
    "DUSP6": "#4daf4a", "SPRY2": "#4daf4a",
    # PI3K/AKT — orange
    "PIK3CA": "#ff7f00", "PIK3CB": "#ff7f00", "AKT1": "#ff7f00",
    "AKT2": "#ff7f00", "MTOR": "#ff7f00",
    # Receptors — purple
    "EGFR": "#984ea3", "ERBB2": "#984ea3",
    # Oncogenes / other — gray
    "MYC": "#999999",
}


# ────────────────────────────────────────────────────────────────
# Data Loading
# ────────────────────────────────────────────────────────────────
def load_matrix(db_path: str = "rimap_results.duckdb"):
    """Load binding matrix following nmf_analysis.py pattern."""
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

    col_means = V_raw.mean(axis=0)
    V_residual = np.clip(V_raw - col_means, 0, None)
    V_residual[V_residual == 0] = 1e-6

    return V_raw, V_residual, mirna_ids, transcript_ids


# ────────────────────────────────────────────────────────────────
# Co-targeting Matrix
# ────────────────────────────────────────────────────────────────
def compute_cotargeting(V, threshold):
    """Build 25×25 co-targeting matrix.

    For each gene pair (i,j), count how many miRNAs have residual
    binding > threshold to BOTH genes.
    """
    n_genes = V.shape[1]
    # Binary mask: which miRNAs target each gene above threshold
    binary = (V > threshold).astype(np.float64)
    print(f"  Residual threshold: {threshold}")
    print(f"  miRNAs above threshold per gene: "
          f"min={binary.sum(axis=0).min():.0f}, "
          f"max={binary.sum(axis=0).max():.0f}, "
          f"mean={binary.sum(axis=0).mean():.0f}")

    co_target = np.zeros((n_genes, n_genes))
    for i in range(n_genes):
        for j in range(i + 1, n_genes):
            shared = (binary[:, i] * binary[:, j]).sum()
            co_target[i, j] = shared
            co_target[j, i] = shared

    return co_target


# ────────────────────────────────────────────────────────────────
# Chord Diagram Drawing
# ────────────────────────────────────────────────────────────────
def draw_chord_diagram(co_target, gene_names, arc_threshold, output_path,
                       figsize=12, dpi=200):
    """Draw chord diagram with matplotlib.

    Places genes on a circle and draws bezier arcs between co-targeted pairs.
    """
    n = len(gene_names)

    # --- filter arcs ---
    mask = co_target > arc_threshold
    scores = co_target[mask]
    if len(scores) == 0:
        print(f"  WARNING: No arcs above threshold {arc_threshold}. "
              f"Max score: {co_target.max():.0f}")
        return
    print(f"  Arcs drawn: {len(scores)} (threshold={arc_threshold})")
    print(f"  Score range: [{scores.min():.0f}, {scores.max():.0f}]")

    # normalize scores for visual mapping
    s_min, s_max = scores.min(), scores.max()
    s_range = s_max - s_min if s_max > s_min else 1.0
    scores_norm = (scores - s_min) / s_range

    # --- node positions on circle ---
    radius = 1.0
    gap = 2 * np.pi * 0.01  # small gap between nodes
    angles = np.linspace(gap, 2 * np.pi - gap, n, endpoint=False)
    # start from top
    angles = angles - np.pi / 2

    xs = radius * np.cos(angles)
    ys = radius * np.sin(angles)

    # --- figure setup ---
    fig, ax = plt.subplots(figsize=(figsize, figsize), subplot_kw={"aspect": "equal"})
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.45, 1.45)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    # --- outer ring segments ---
    segment_width = 2 * np.pi / n * 0.8
    for i, gene in enumerate(gene_names):
        color = PATHWAY_COLORS.get(gene, "#999999")
        arc_angles = np.linspace(
            angles[i] - segment_width / 2,
            angles[i] + segment_width / 2,
            30,
        )
        inner_r = 0.92
        outer_r = 1.02
        x_inner = inner_r * np.cos(arc_angles)
        y_inner = inner_r * np.sin(arc_angles)
        x_outer = outer_r * np.cos(arc_angles[::-1])
        y_outer = outer_r * np.sin(arc_angles[::-1])
        ax.fill(
            np.concatenate([x_inner, x_outer]),
            np.concatenate([y_inner, y_outer]),
            color=color,
            edgecolor="white",
            linewidth=0.5,
            zorder=3,
        )

    # --- gene labels ---
    label_r = 1.18
    for i, gene in enumerate(gene_names):
        lx = label_r * np.cos(angles[i])
        ly = label_r * np.sin(angles[i])
        color = PATHWAY_COLORS.get(gene, "#999999")

        angle = angles[i]
        is_upper = ly >= 0
        if is_upper:
            rot = np.degrees(angle) - 90
        else:
            rot = np.degrees(angle) + 90

        ha = "center"
        va = "center"
        ax.text(
            lx, ly, gene, fontsize=7, fontweight="bold",
            color=color, ha=ha, va=va, rotation=rot,
            rotation_mode="anchor",
            path_effects=[pe.withStroke(linewidth=1.5, foreground="white")],
            zorder=5,
        )

    # --- draw arcs ---
    # collect arcs sorted by score (draw smallest first)
    arc_pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            if co_target[i, j] > arc_threshold:
                arc_pairs.append((i, j, co_target[i, j]))
    arc_pairs.sort(key=lambda x: x[2])

    for idx, (i, j, score) in enumerate(arc_pairs):
        sn = (score - s_min) / s_range
        alpha = 0.12 + 0.58 * sn
        lw = 0.4 + 4.0 * sn

        # bezier control point: midpoint pulled toward center
        mx = (xs[i] + xs[j]) / 2
        my = (ys[i] + ys[j]) / 2
        pull = 0.4 + 0.3 * sn
        cx = mx * pull
        cy = my * pull

        # cubic bezier through Path
        verts = [
            (xs[i], ys[i]),
            (cx, cy),
            (cx, cy),
            (xs[j], ys[j]),
        ]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        path = Path(verts, codes)

        # gradient: strong arcs → deep red, weak → pale gray (but keep pathway tint)
        strong_rgb = (0.80, 0.15, 0.15)
        weak_rgb = (0.75, 0.75, 0.75)
        r = weak_rgb[0] + (strong_rgb[0] - weak_rgb[0]) * sn
        g = weak_rgb[1] + (strong_rgb[1] - weak_rgb[1]) * sn
        b = weak_rgb[2] + (strong_rgb[2] - weak_rgb[2]) * sn

        c_i = PATHWAY_COLORS.get(gene_names[i], "#999999")
        c_j = PATHWAY_COLORS.get(gene_names[j], "#999999")
        import matplotlib.colors as mcolors
        r_i, g_i, b_i, _ = mcolors.to_rgba(c_i)
        r_j, g_j, b_j, _ = mcolors.to_rgba(c_j)
        blend = 0.35 * sn
        r = r + (r_i + r_j) / 2 * blend
        g = g + (g_i + g_j) / 2 * blend
        b = b + (b_i + b_j) / 2 * blend
        r, g, b = min(r, 1), min(g, 1), min(b, 1)

        arc_color = (r, g, b, alpha)

        patch = mpatches.PathPatch(
            path, facecolor="none", edgecolor=arc_color,
            linewidth=lw, zorder=2,
        )
        ax.add_patch(patch)

        # subtle shadow for top arcs
        if sn > 0.7:
            shadow_verts = [
                (xs[i] + 0.005, ys[i] - 0.005),
                (cx + 0.005, cy - 0.005),
                (cx + 0.005, cy - 0.005),
                (xs[j] + 0.005, ys[j] - 0.005),
            ]
            shadow_codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
            shadow_path = Path(shadow_verts, shadow_codes)
            shadow_patch = mpatches.PathPatch(
                shadow_path, facecolor="none",
                edgecolor=(0.0, 0.0, 0.0, 0.08),
                linewidth=lw + 1.5, zorder=1,
            )
            ax.add_patch(shadow_patch)

    # --- legend ---
    legend_items = [
        mpatches.Patch(color="#e41a1c", label="RAS family"),
        mpatches.Patch(color="#377eb8", label="Tumor suppressors"),
        mpatches.Patch(color="#4daf4a", label="MAPK cascade"),
        mpatches.Patch(color="#ff7f00", label="PI3K/AKT"),
        mpatches.Patch(color="#984ea3", label="Receptors"),
        mpatches.Patch(color="#999999", label="Other"),
    ]
    leg = ax.legend(
        handles=legend_items, loc="lower center",
        ncol=3, fontsize=8, frameon=True,
        fancybox=True, shadow=False,
        bbox_to_anchor=(0.5, -0.02),
    )
    leg.get_frame().set_alpha(0.9)

    ax.set_title(
        "Gene–Gene Co-targeting by miRNAs\n",
        fontsize=14, fontweight="bold", pad=25,
    )
    ax.text(
        0.5, 1.01,
        "Arc width ∝ number of miRNAs that bind both genes above their transcript-level mean affinity.\n"
        "Hub genes: KRAS, HRAS, PIK3CA.",
        fontsize=8.5, ha="center", va="top",
        color="#555555", transform=ax.transAxes,
        linespacing=1.6,
    )

    # --- save ---
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def _blend_colors(hex1, hex2, alpha=0.5):
    """Blend two hex colors and return rgba tuple."""
    import matplotlib.colors as mcolors
    r1, g1, b1, _ = mcolors.to_rgba(hex1)
    r2, g2, b2, _ = mcolors.to_rgba(hex2)
    return ((r1 + r2) / 2, (g1 + g2) / 2, (b1 + b2) / 2, alpha)


# ────────────────────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────────────────────
def main():
    """CLI entry: load matrix -> co-targeting scores -> draw chord diagram."""
    parser = argparse.ArgumentParser(
        description="Chord diagram of gene-gene co-targeting by miRNAs"
    )
    parser.add_argument("--threshold", type=float, default=100,
                        help="Minimum co-targeting score to draw an arc (default: 100)")
    parser.add_argument("--use-residual", action="store_true", default=True,
                        help="Use residual matrix (default: True)")
    parser.add_argument("--no-residual", dest="use_residual", action="store_false",
                        help="Use raw matrix instead of residual")
    parser.add_argument("--residual-threshold", type=float, default=0.5,
                        help="Min residual binding for miRNA to count as targeting (default: 0.5)")
    parser.add_argument("--output", type=str, default="figures/chord_diagram.png",
                        help="Output file path (default: figures/chord_diagram.png)")
    parser.add_argument("--dpi", type=int, default=200,
                        help="DPI for output (default: 200)")
    parser.add_argument("--figsize", type=float, default=12,
                        help="Figure size in inches, square (default: 12)")
    args = parser.parse_args()

    print("Loading data...")
    V_raw, V_residual, mirna_ids, transcript_ids = load_matrix()

    V = V_residual if args.use_residual else V_raw
    matrix_type = "residual" if args.use_residual else "raw"
    print(f"\nUsing {matrix_type} matrix")

    # Map transcript IDs to gene names (preserve column order)
    gene_names = [TRANSCRIPT_GENE_MAP.get(tid, tid) for tid in transcript_ids]
    print(f"Genes: {len(gene_names)}")

    print("\nComputing co-targeting matrix...")
    co_target = compute_cotargeting(V, args.residual_threshold)

    # Print top pairs
    print("\nTop 10 co-targeted gene pairs:")
    pairs = []
    for i in range(len(gene_names)):
        for j in range(i + 1, len(gene_names)):
            pairs.append((gene_names[i], gene_names[j], co_target[i, j]))
    pairs.sort(key=lambda x: x[2], reverse=True)
    for g1, g2, s in pairs[:10]:
        print(f"  {g1:>8} — {g2:<8} : {s:.0f}")

    print("\nDrawing chord diagram...")
    draw_chord_diagram(
        co_target, gene_names,
        arc_threshold=args.threshold,
        output_path=args.output,
        figsize=args.figsize,
        dpi=args.dpi,
    )
    print("Done.")


if __name__ == "__main__":
    main()

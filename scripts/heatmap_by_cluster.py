#!/usr/bin/env python3
"""
Binding Matrix Heatmap Sorted by NMF Signature.

Creates a full 2516×25 heatmap where rows (miRNAs) are sorted by
dominant NMF signature and membership weight, and columns (genes)
are sorted by H matrix loading. Produces a clear block-diagonal
structure showing how the two signatures separate binding profiles.

Usage:
    .venv/bin/python3 scripts/heatmap_by_cluster.py
    .venv/bin/python3 scripts/heatmap_by_cluster.py --output figures/heatmap.png --dpi 300
    .venv/bin/python3 scripts/heatmap_by_cluster.py --no-use-residual
    .venv/bin/python3 scripts/heatmap_by_cluster.py --show-mixed
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
from matplotlib.colors import ListedColormap
import seaborn as sns
from scipy.stats import entropy


# ────────────────────────────────────────────────────────────────
# Constants
# ────────────────────────────────────────────────────────────────
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

PATHWAY_GROUPS = {
    "RAS/MAPK":  {"RAF1", "BRAF", "MAP2K1", "KRAS", "NRAS", "HRAS",
                  "MAP2K2", "MAPK3", "MAPK1"},
    "RTK":       {"EGFR", "ERBB2"},
    "PI3K/AKT":  {"PIK3CA", "PIK3CB", "AKT1", "AKT2", "MTOR", "PTEN"},
    "Regulator": {"NF1", "DUSP6", "SPRY2"},
    "Supp/Cell": {"TP53", "RB1", "APC", "CDKN1A"},
    "Oncogene":  {"MYC"},
}

PATHWAY_COLORS = {
    "RAS/MAPK":  "#e74c3c",
    "RTK":       "#e67e22",
    "PI3K/AKT":  "#2980b9",
    "Regulator": "#27ae60",
    "Supp/Cell": "#8e44ad",
    "Oncogene":  "#f39c12",
}

SIGNATURE_PALETTE = ["#3498db", "#e74c3c"]  # Sig 1 blue, Sig 2 red


def gene_pathway(gene: str) -> str:
    """Classify a gene into a pathway group (kras_pathway, tumor_suppressor, mapk_pathway, Other).

    Args:
        gene: HGNC gene symbol.

    Returns:
        The pathway group name.
    """
    for group, genes in PATHWAY_GROUPS.items():
        if gene in genes:
            return group
    return "Other"


# ────────────────────────────────────────────────────────────────
# Data Loader (mirrors nmf_analysis.py)
# ────────────────────────────────────────────────────────────────
def load_binding_matrix(
    db_path: str = "rimap_results.duckdb", use_residual: bool = True
) -> tuple[np.ndarray, list[str], list[str]]:
    """Load the miRNA×transcript binding matrix from DuckDB.

    Returns:
        V: matrix (n_mirna, n_transcript)
        mirna_ids: row labels
        transcript_ids: column labels
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

    if use_residual:
        col_means = V_raw.mean(axis=0)
        V = np.clip(V_raw - col_means, 0, None)
        V[V == 0] = 1e-6
    else:
        V = V_raw

    return V, mirna_ids, transcript_ids


# ────────────────────────────────────────────────────────────────
# Ordering logic
# ────────────────────────────────────────────────────────────────
def order_rows_by_signature(
    W: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return row indices sorted by dominant signature then weight.

    Returns:
        row_order: sorted row indices
        dominant_sig: dominant signature per row (original order)
        sig_weight: membership weight for dominant signature (original order)
    """
    dominant_sig = W.argmax(axis=1)
    sig_weight = W[np.arange(len(W)), dominant_sig]

    row_order = np.lexsort((-sig_weight, dominant_sig))
    return row_order, dominant_sig, sig_weight


def order_cols_by_signature(
    H: np.ndarray,
) -> np.ndarray:
    """Return column indices sorted: Sig-1-high genes first, then Sig-2-high.

    Within each block, sorted by the loading to that signature (descending).
    """
    dominant_col_sig = H.argmax(axis=0)
    col_loading = H[dominant_col_sig, np.arange(H.shape[1])]

    col_order = np.lexsort((-col_loading, dominant_col_sig))
    return col_order


def detect_mixed(W: np.ndarray, threshold: float = 0.8) -> np.ndarray:
    """Boolean mask for miRNAs with high normalized entropy."""
    row_sums = W.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    W_prob = W / row_sums
    row_ent = entropy(W_prob.T)
    max_ent = np.log(W.shape[1])
    if max_ent == 0:
        return np.zeros(W.shape[0], dtype=bool)
    normalized = np.nan_to_num(row_ent / max_ent, nan=0.0)
    return normalized > threshold


# ────────────────────────────────────────────────────────────────
# Plotting
# ────────────────────────────────────────────────────────────────
def draw_heatmap(
    V: np.ndarray,
    W: np.ndarray,
    H: np.ndarray,
    gene_names: list[str],
    output_path: str,
    dpi: int,
    show_mixed: bool,
):
    """Draw the sorted binding matrix heatmap with signature color bar."""
    n_mirna, n_genes = V.shape
    k = W.shape[1]

    # ── Ordering ────────────────────────────────────────────────
    row_order, dominant_sig, sig_weight = order_rows_by_signature(W)
    col_order = order_cols_by_signature(H)

    V_sorted = V[row_order][:, col_order]
    dominant_sorted = dominant_sig[row_order]

    gene_names_ordered = [gene_names[i] for i in col_order]

    # Build signature color bar (1-D, one color per miRNA row)
    sig_colors = np.array([SIGNATURE_PALETTE[s] for s in dominant_sorted])

    # ── Mixed markers ───────────────────────────────────────────
    is_mixed = detect_mixed(W)
    is_mixed_sorted = is_mixed[row_order]

    # ── Figure layout ───────────────────────────────────────────
    fig_width = max(14, n_genes * 0.6 + 3)
    fig_height = max(10, n_mirna * 0.005 + 2)

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = fig.add_gridspec(
        1, 3,
        width_ratios=[0.25, 0.15, n_genes * 0.5],
        wspace=0.02,
    )

    ax_cbar = fig.add_subplot(gs[0, 0])   # signature color bar
    ax_gap = fig.add_subplot(gs[0, 1])     # gap / mixed markers
    ax_heat = fig.add_subplot(gs[0, 2])    # heatmap

    # ── Signature color bar ─────────────────────────────────────
    sig_bar = np.array(dominant_sorted).reshape(-1, 1).astype(float)
    sig_cmap = ListedColormap(SIGNATURE_PALETTE[:k])
    ax_cbar.imshow(sig_bar, aspect="auto", cmap=sig_cmap,
                   vmin=0, vmax=k - 1, interpolation="nearest")
    ax_cbar.set_xticks([])
    ax_cbar.set_yticks([])
    ax_cbar.set_ylabel("")
    ax_cbar.set_title("Sig.", fontsize=9, pad=4)
    ax_cbar.tick_params(left=False, bottom=False)

    # Draw boundary line between signatures
    sig1_count = (dominant_sorted == 0).sum()
    if 0 < sig1_count < n_mirna:
        boundary = sig1_count - 0.5
        ax_cbar.axhline(boundary, color="white", linewidth=1.5)
        ax_heat.axhline(boundary, color="white", linewidth=1.0, alpha=0.6)

    # ── Mixed markers column ────────────────────────────────────
    ax_gap.set_xlim(0, 1)
    ax_gap.set_ylim(0, n_mirna)
    ax_gap.axis("off")

    if show_mixed and is_mixed_sorted.any():
        mixed_y = np.where(is_mixed_sorted)[0]
        ax_gap.scatter(
            np.full(len(mixed_y), 0.5),
            n_mirna - mixed_y - 0.5,  # imshow origin=top
            marker="|",
            s=15,
            c="black",
            linewidths=0.5,
            zorder=5,
        )
        ax_gap.set_title("Mix", fontsize=8, pad=4)

    # ── Heatmap ─────────────────────────────────────────────────
    vmax = np.percentile(V_sorted[V_sorted > 0], 97)
    im = ax_heat.imshow(
        V_sorted,
        aspect="auto",
        cmap="inferno",
        vmin=0,
        vmax=vmax,
        interpolation="nearest",
    )

    ax_heat.set_yticks([])
    ax_heat.set_ylabel(f"{n_mirna} miRNAs", fontsize=11)

    # Column labels with pathway coloring
    ax_heat.set_xticks(range(n_genes))
    xlabels = []
    for gn in gene_names_ordered:
        pw = gene_pathway(gn)
        color = PATHWAY_COLORS.get(pw, "#333333")
        xlabels.append(ax_heat.text(
            0, 0, gn, rotation=45, ha="right", va="top",
            fontsize=8, color=color,
        ))

    ax_heat.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=8)

    # Color each tick label by pathway
    for i, label in enumerate(ax_heat.get_xticklabels()):
        pw = gene_pathway(gene_names_ordered[i])
        label.set_color(PATHWAY_COLORS.get(pw, "#333333"))

    # ── Colorbar for values ─────────────────────────────────────
    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.02, pad=0.02)
    cbar.set_label("Binding weight (|dG| residual)", fontsize=10)

    # ── Legend ───────────────────────────────────────────────────
    legend_handles = []
    for sig_idx in range(k):
        n_sig = (dominant_sorted == sig_idx).sum()
        patch = mpatches.Patch(
            color=SIGNATURE_PALETTE[sig_idx],
            label=f"Signature {sig_idx + 1} ({n_sig})",
        )
        legend_handles.append(patch)

    if show_mixed:
        n_mix = is_mixed_sorted.sum()
        legend_handles.append(
            mpatches.Patch(color="black", label=f"Mixed ({n_mix})")
        )

    # Pathway legend
    for pw, color in PATHWAY_COLORS.items():
        legend_handles.append(mpatches.Patch(color=color, label=pw))

    ax_heat.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(1.18, 1.0),
        fontsize=8,
        frameon=True,
        framealpha=0.9,
    )

    # ── Title ───────────────────────────────────────────────────
    ax_heat.set_title(
        "miRNA–Transcript Binding Matrix (sorted by NMF signature)",
        fontsize=13, pad=10,
    )

    fig.tight_layout()
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ────────────────────────────────────────────────────────────────
# CLI
# ────────────────────────────────────────────────────────────────
def parse_args() -> argparse.Namespace:
    """Parse CLI args: --output, --dpi, --use-residual/--no-use-residual, --show-mixed.

    Returns:
        The parsed ``argparse.Namespace``.
    """
    parser = argparse.ArgumentParser(
        description="Binding matrix heatmap sorted by NMF signature"
    )
    parser.add_argument(
        "--output", type=str, default="figures/heatmap_sorted_nmf.png",
        help="Output path for the heatmap figure",
    )
    parser.add_argument(
        "--dpi", type=int, default=200,
        help="DPI for the output figure",
    )
    parser.add_argument(
        "--use-residual", action=argparse.BooleanOptionalAction,
        default=True,
        help="Use residual matrix (subtract column means, clip to ≥0). "
             "Use --no-use-residual for raw matrix.",
    )
    parser.add_argument(
        "--show-mixed", action="store_true",
        help="Highlight mixed miRNAs (entropy > 0.8) with markers",
    )
    return parser.parse_args()


def main():
    """CLI entry: load V (residual), order rows/cols by NMF signature, draw heatmap."""
    args = parse_args()

    plt.style.use("seaborn-v0_8-whitegrid")

    # ── Load binding matrix ─────────────────────────────────────
    print("Loading binding matrix from DuckDB ...")
    V, mirna_ids, transcript_ids = load_binding_matrix(
        use_residual=args.use_residual
    )
    gene_names = [TRANSCRIPT_GENE_MAP.get(tid, tid) for tid in transcript_ids]
    print(f"  Matrix: {V.shape[0]} miRNAs × {V.shape[1]} transcripts")

    # ── Load NMF factors ────────────────────────────────────────
    print("Loading NMF factors ...")
    W = np.load("figures_nmf/W_matrix.npy")
    H = np.load("figures_nmf/H_matrix.npy")
    print(f"  W: {W.shape},  H: {H.shape}")

    assert W.shape[0] == V.shape[0], (
        f"Row mismatch: W has {W.shape[0]}, matrix has {V.shape[0]}"
    )
    assert H.shape[1] == V.shape[1], (
        f"Column mismatch: H has {H.shape[1]}, matrix has {V.shape[1]}"
    )

    # ── Draw ────────────────────────────────────────────────────
    print("Drawing heatmap ...")
    draw_heatmap(
        V=V,
        W=W,
        H=H,
        gene_names=gene_names,
        output_path=args.output,
        dpi=args.dpi,
        show_mixed=args.show_mixed,
    )

    print("Done.")


if __name__ == "__main__":
    main()

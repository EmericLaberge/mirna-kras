#!/usr/bin/env python3
"""
NMF Analysis Pipeline for miRNA-Transcript Binding Data.

Non-Negative Matrix Factorization on the 2516×25 binding matrix
to discover miRNA binding signatures via soft clustering.

Usage:
    .venv/bin/python3 scripts/nmf_analysis.py --output-dir figures_nmf --n-seeds 20
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
import seaborn as sns
from sklearn.decomposition import NMF
from scipy.cluster.hierarchy import linkage, cophenet
from scipy.spatial.distance import squareform
from scipy.stats import entropy
from time import perf_counter

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


# ────────────────────────────────────────────────────────────────
# Data Loader
# ────────────────────────────────────────────────────────────────
class NMFDataLoader:
    """Load miRNA-transcript binding matrix from DuckDB."""

    def __init__(self, db_path: str = "rimap_results.duckdb"):
        """Path to the DuckDB file. Read-only, never modified by this loader."""
        self.db_path = db_path

    def load(self) -> tuple[np.ndarray, np.ndarray, list[str], list[str]]:
        """Read DuckDB, filter seedless, aggregate mean |dG|, return (V_residual, V_raw, mirna_ids, transcript_ids)."""
        con = duckdb.connect(self.db_path, read_only=True)
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
        print(f"  Value range: [{V_raw.min():.3f}, {V_raw.max():.3f}]")
        print(f"  Density: {(V_raw > 0).sum() / V_raw.size * 100:.1f}%")

        col_means = V_raw.mean(axis=0)
        V_residual = np.clip(V_raw - col_means, 0, None)
        V_residual[V_residual == 0] = 1e-6

        nonzero_frac = (V_residual > 1e-6).sum() / V_residual.size
        print(f"  Residual matrix non-zero fraction: {nonzero_frac:.1%}")
        print(f"  Residual value range: [{V_residual[V_residual > 1e-6].min():.4f}, {V_residual.max():.4f}]")

        return V_residual, V_raw, mirna_ids, transcript_ids


# ────────────────────────────────────────────────────────────────
# Rank Selector
# ────────────────────────────────────────────────────────────────
class NMFRankSelector:
    """Evaluate NMF stability across ranks via consensus clustering."""

    def __init__(self, k_range: range = range(2, 9), n_seeds: int = 20, max_iter: int = 500):
        """``n_seeds`` random factorizations per k for the consensus matrix."""
        self.k_range = k_range
        self.n_seeds = n_seeds
        self.max_iter = max_iter

    def evaluate(self, V: np.ndarray) -> dict:
        """Run multi-seed NMF for each k and compute consensus metrics.

        Returns:
            dict mapping k → {cophenetic, error, consensus, W, H}
        """
        n = V.shape[0]
        results = {}

        total_k = len(self.k_range)
        for idx, k in enumerate(self.k_range, 1):
            t0 = perf_counter()
            connectivity = np.zeros((n, n), dtype=np.float64)
            errors = []

            for seed in range(self.n_seeds):
                model = NMF(
                    n_components=k,
                    init="random",
                    solver="mu",
                    beta_loss="frobenius",
                    max_iter=self.max_iter,
                    random_state=seed,
                    tol=1e-4,
                )
                W = model.fit_transform(V)
                H = model.components_
                errors.append(model.reconstruction_err_)

                # Vectorized connectivity
                labels = W.argmax(axis=1)
                same_cluster = labels[:, None] == labels[None, :]
                connectivity += same_cluster.astype(np.float64)

            consensus = connectivity / self.n_seeds

            # Cophenetic correlation from consensus matrix
            dist_matrix = 1.0 - consensus
            np.fill_diagonal(dist_matrix, 0.0)
            condensed = squareform(dist_matrix)
            Z = linkage(condensed, method="average")
            coph, _ = cophenet(Z, condensed)

            avg_error = np.mean(errors)

            # Refit one model for the stored W, H
            best_model = NMF(
                n_components=k,
                init="random",
                solver="mu",
                beta_loss="frobenius",
                max_iter=self.max_iter,
                random_state=0,
                tol=1e-4,
            )
            W_best = best_model.fit_transform(V)
            H_best = best_model.components_

            results[k] = {
                "cophenetic": coph,
                "error": avg_error,
                "consensus": consensus,
                "W": W_best,
                "H": H_best,
            }

            elapsed = perf_counter() - t0
            print(f"    k={k}: cophenetic={coph:.4f}, error={avg_error:.2f}  ({elapsed:.1f}s)")

        return results


# ────────────────────────────────────────────────────────────────
# Analyzer
# ────────────────────────────────────────────────────────────────
class NMFAnalyzer:
    """Final NMF fit and biological interpretation."""

    def __init__(self, n_components: int, max_iter: int = 1000):
        """``n_components=k`` is the chosen NMF rank from the rank selection phase."""
        self.n_components = n_components
        self.max_iter = max_iter

    def fit(self, V: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
        """Fit final NMF with NNDSVD initialization.

        Returns:
            W: membership matrix (n_mirna, k)
            H: signature matrix (k, n_transcripts)
            reconstruction_error: Frobenius norm of V - WH
        """
        model = NMF(
            n_components=self.n_components,
            init="nndsvd",
            solver="mu",
            beta_loss="frobenius",
            max_iter=self.max_iter,
            random_state=42,
        )
        W = model.fit_transform(V)
        H = model.components_
        return W, H, model.reconstruction_err_

    @staticmethod
    def get_top_transcripts(
        H: np.ndarray, transcript_ids: list[str], n_top: int = 5
    ) -> list[list[tuple[str, float]]]:
        """For each signature (row of H), return top-n transcripts by weight."""
        results = []
        for i in range(H.shape[0]):
            idx = np.argsort(H[i])[::-1][:n_top]
            top = [(transcript_ids[j], H[i, j]) for j in idx]
            results.append(top)
        return results

    def detect_mixed(W: np.ndarray, threshold: float = 0.8) -> tuple[np.ndarray, np.ndarray]:
        """Return (normalized_entropy_per_row, is_mixed_mask). A row with H_norm > threshold is "mixed"."""
        row_sums = W.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        W_prob = W / row_sums
        row_ent = entropy(W_prob.T)
        max_ent = np.log(W.shape[1])
        if max_ent == 0:
            return np.zeros(W.shape[0]), np.zeros(W.shape[0], dtype=bool)
        normalized = np.nan_to_num(row_ent / max_ent, nan=0.0)
        return normalized, normalized > threshold


# ────────────────────────────────────────────────────────────────
# Plotter
# ────────────────────────────────────────────────────────────────
class NMFPlotter:
    """Publication-quality NMF visualizations."""

    @staticmethod
    def rank_selection_curves(results: dict, optimal_k: int, output_dir: str):
        """Reconstruction error and cophenetic correlation vs k."""
        ks = sorted(results.keys())
        errors = [results[k]["error"] for k in ks]
        cophenetics = [results[k]["cophenetic"] for k in ks]

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        ax1.plot(ks, errors, "o-", color="#2c3e50", linewidth=2, markersize=8)
        ax1.set_xlabel("Number of signatures (k)", fontsize=12)
        ax1.set_ylabel("Reconstruction error (Frobenius)", fontsize=12)
        ax1.set_title("Reconstruction Error vs Rank", fontsize=13)
        ax1.set_xticks(ks)

        ax2.plot(ks, cophenetics, "s-", color="#8e44ad", linewidth=2, markersize=8)
        ax2.set_xlabel("Number of signatures (k)", fontsize=12)
        ax2.set_ylabel("Cophenetic correlation", fontsize=12)
        ax2.set_title("Cophenetic Correlation vs Rank", fontsize=13)
        ax2.set_xticks(ks)

        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, "rank_selection.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)

    @staticmethod
    def signature_heatmap(H: np.ndarray, transcript_ids: list[str], output_dir: str):
        """Heatmap of H matrix (k signatures × 25 transcripts)."""
        fig, ax = plt.subplots(figsize=(14, 6))

        # Sort transcripts globally by total weight across signatures
        order = np.argsort(H.sum(axis=0))[::-1]
        H_sorted = H[:, order]
        labels_sorted = [transcript_ids[i] for i in order]

        sns.heatmap(
            H_sorted,
            xticklabels=labels_sorted,
            yticklabels=[f"Signature {i+1}" for i in range(H.shape[0])],
            annot=True,
            fmt=".3f",
            cmap="YlOrRd",
            linewidths=0.5,
            ax=ax,
        )
        ax.set_title("NMF Signature × Transcript Weights", fontsize=14)
        ax.set_xlabel("Transcript", fontsize=12)
        ax.set_ylabel("Signature", fontsize=12)
        plt.xticks(rotation=45, ha="right")

        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, "signature_heatmap.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)

    @staticmethod
    def membership_scatter(
        W: np.ndarray,
        mirna_ids: list[str],
        is_mixed: np.ndarray,
        output_dir: str,
    ):
        """2-D scatter of miRNAs colored by dominant signature, with mixed in a separate panel.

        Args:
            W: Membership matrix (n_mirna, k) from the final NMF fit.
            mirna_ids: List of miRNA IDs in the same order as W rows.
            is_mixed: Boolean mask of length n_mirna marking mixed miRNAs.
            output_dir: Directory where the PNG will be written.
        """
        k = W.shape[1]
        palette = sns.color_palette("tab10", min(k, 10))
        if k > 10:
            palette = palette + sns.color_palette("tab10", k - 10)
        dominant = W.argmax(axis=1)

        if k > 2:
            from umap import UMAP
            from sklearn.preprocessing import StandardScaler
            W_scaled = StandardScaler().fit_transform(W)
            W_2d = UMAP(n_components=2, n_neighbors=30, min_dist=0.1, random_state=42).fit_transform(W_scaled)
            xlabel, ylabel = "UMAP 1", "UMAP 2"
        else:
            W_2d = W[:, :2]
            xlabel, ylabel = "Signature 1 weight", "Signature 2 weight"

        fig, ax = plt.subplots(figsize=(12, 8))

        for sig in range(k):
            mask = (dominant == sig) & (~is_mixed)
            marker = "^" if sig >= 10 else "o"
            ax.scatter(
                W_2d[mask, 0],
                W_2d[mask, 1],
                c=[palette[sig]],
                marker=marker,
                label=f"Signature {sig+1} ({mask.sum()})",
                alpha=0.6,
                s=20,
            )

        if is_mixed.any():
            ax.scatter(
                W_2d[is_mixed, 0],
                W_2d[is_mixed, 1],
                c="gray",
                marker="x",
                label=f"Mixed ({is_mixed.sum()})",
                alpha=0.4,
                s=30,
            )

        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_title("miRNA Membership in NMF Signatures", fontsize=14)
        ax.legend(fontsize=10, markerscale=2)

        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, "membership_scatter.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)

    @staticmethod
    def signature_bars(H: np.ndarray, transcript_ids: list[str], output_dir: str):
        """Horizontal bar charts: one per signature showing top-10 transcripts."""
        palette = sns.color_palette("tab10", H.shape[0])
        n_top = min(10, H.shape[1])

        for i in range(H.shape[0]):
            idx = np.argsort(H[i])[::-1][:n_top]
            names = [transcript_ids[j] for j in idx]
            values = H[i, idx]

            fig, ax = plt.subplots(figsize=(6, 4))
            ax.barh(range(n_top), values, color=palette[i], edgecolor="black", linewidth=0.5)
            ax.set_yticks(range(n_top))
            ax.set_yticklabels(names)
            ax.invert_yaxis()
            ax.set_xlabel("Weight", fontsize=11)
            ax.set_title(f"Signature {i+1} — Top {n_top} Transcripts", fontsize=12)

            fig.tight_layout()
            fig.savefig(
                os.path.join(output_dir, f"signature_{i+1}_bars.png"),
                dpi=200,
                bbox_inches="tight",
            )
            plt.close(fig)


# ────────────────────────────────────────────────────────────────
# CLI
# ────────────────────────────────────────────────────────────────
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the NMF analysis pipeline.

    Returns:
        The parsed arguments: --output-dir, --k-min, --k-max, --n-seeds,
        --max-iter, --mixed-threshold.
    """
    parser = argparse.ArgumentParser(
        description="NMF analysis pipeline for miRNA-transcript binding data"
    )
    parser.add_argument("--output-dir", type=str, default="figures_nmf",
                        help="Output directory for figures and data")
    parser.add_argument("--k-min", type=int, default=2,
                        help="Minimum rank to evaluate")
    parser.add_argument("--k-max", type=int, default=8,
                        help="Maximum rank to evaluate")
    parser.add_argument("--n-seeds", type=int, default=20,
                        help="Number of random seeds per rank")
    parser.add_argument("--max-iter", type=int, default=500,
                        help="Max iterations per NMF fit")
    parser.add_argument("--mixed-threshold", type=float, default=0.8,
                        help="Normalized entropy threshold for 'mixed' miRNAs")
    return parser.parse_args()

def main():
    """CLI entry: load DuckDB -> residual -> rank selection -> final NMF -> figures + W/H matrices."""

# ────────────────────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────────────────────
def main():
    """CLI entry: load DuckDB -> residual -> rank selection -> final NMF -> figures + W/H matrices."""
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    plt.style.use("seaborn-v0_8-whitegrid")

    separator = "=" * 60

    # ── 1. Load data ───────────────────────────────────────────
    print(f"\n  {separator}")
    print(f"  1. Loading data from DuckDB")
    print(f"  {separator}")
    loader = NMFDataLoader()
    X, X_raw, mirna_ids, transcript_ids = loader.load()
    gene_names = [TRANSCRIPT_GENE_MAP.get(tid, tid) for tid in transcript_ids]

    # ── 2. Rank selection ──────────────────────────────────────
    k_range = range(args.k_min, args.k_max + 1)
    print(f"\n  {separator}")
    print(f"  2. Rank selection (k={args.k_min}..{args.k_max}, {args.n_seeds} seeds each)")
    print(f"  {separator}")
    t0 = perf_counter()
    selector = NMFRankSelector(k_range=k_range, n_seeds=args.n_seeds, max_iter=args.max_iter)
    rank_results = selector.evaluate(X)
    elapsed_rank = perf_counter() - t0
    print(f"  Total rank selection time: {elapsed_rank:.1f}s")

    # Pick optimal k by cophenetic correlation (highest)
    optimal_k = max(rank_results, key=lambda k: rank_results[k]["cophenetic"])

    # Print rank selection table
    print(f"\n  {'k':>3}  {'Cophenetic':>12}  {'Error':>12}")
    print(f"  {'─'*3}  {'─'*12}  {'─'*12}")
    for k in sorted(rank_results.keys()):
        marker = " ◄ optimal" if k == optimal_k else ""
        print(f"  {k:>3}  {rank_results[k]['cophenetic']:>12.4f}  {rank_results[k]['error']:>12.2f}{marker}")

    # ── 3. Final NMF fit ──────────────────────────────────────
    print(f"\n  {separator}")
    print(f"  3. Fitting final NMF at k={optimal_k} (NNDSVD init)")
    print(f"  {separator}")
    analyzer = NMFAnalyzer(n_components=optimal_k)
    W, H, recon_err = analyzer.fit(X)
    print(f"  Reconstruction error: {recon_err:.2f}")

    # ── 4. Generate plots ─────────────────────────────────────
    print(f"\n  {separator}")
    print(f"  4. Generating figures")
    print(f"  {separator}")

    NMFPlotter.rank_selection_curves(rank_results, optimal_k, args.output_dir)
    print(f"    ✓ rank_selection.png")

    NMFPlotter.signature_heatmap(H, gene_names, args.output_dir)
    print(f"    ✓ signature_heatmap.png")

    _, is_mixed = NMFAnalyzer.detect_mixed(W, threshold=args.mixed_threshold)
    NMFPlotter.membership_scatter(W, mirna_ids, is_mixed, args.output_dir)
    print(f"    ✓ membership_scatter.png")

    NMFPlotter.signature_bars(H, gene_names, args.output_dir)
    for i in range(H.shape[0]):
        print(f"    ✓ signature_{i+1}_bars.png")

    # ── 5. Signature interpretation ────────────────────────────
    print(f"\n  {separator}")
    print(f"  5. Signature interpretation (k={optimal_k})")
    print(f"  {separator}")

    top_transcripts = NMFAnalyzer.get_top_transcripts(H, gene_names, n_top=5)
    for i, sig in enumerate(top_transcripts):
        print(f"\n    Signature {i+1}:")
        for name, weight in sig:
            print(f"      {name:>10s}  {weight:.4f}")

    # ── 5b. Biological interpretation with raw values ─────────
    print(f"\n  {separator}")
    print(f"  5b. Raw binding strength per signature (mean |dG|)")
    print(f"  {separator}")
    dominant = W.argmax(axis=1)
    for sig in range(optimal_k):
        mask = dominant == sig
        if mask.sum() == 0:
            continue
        mean_binding = X_raw[mask].mean(axis=0)
        idx = np.argsort(mean_binding)[::-1]
        print(f"\n    Signature {sig+1} ({mask.sum()} miRNAs) — mean |dG| per transcript:")
        for rank, j in enumerate(idx[:8], 1):
            print(f"      {rank}. {gene_names[j]:>10s}  {mean_binding[j]:.2f}")

    # ── 6. Cluster sizes ───────────────────────────────────────
    print(f"\n  {separator}")
    print(f"  6. Cluster sizes")
    print(f"  {separator}")
    dominant = W.argmax(axis=1)
    for i in range(optimal_k):
        count = (dominant == i).sum()
        pct = count / len(dominant) * 100
        print(f"    Signature {i+1}: {count:,} miRNAs ({pct:.1f}%)")

    # ── 7. Mixed miRNAs ────────────────────────────────────────
    print(f"\n  {separator}")
    print(f"  7. Mixed miRNAs (entropy > {args.mixed_threshold})")
    print(f"  {separator}")
    norm_ent, is_mixed_arr = NMFAnalyzer.detect_mixed(W, threshold=args.mixed_threshold)
    n_mixed = is_mixed_arr.sum()
    print(f"    Mixed miRNAs: {n_mixed:,} / {len(mirna_ids):,} ({n_mixed/len(mirna_ids)*100:.1f}%)")
    print(f"    Mean entropy: {norm_ent.mean():.4f}")

    # ── 8. Save matrices ───────────────────────────────────────
    print(f"\n  {separator}")
    print(f"  8. Saving matrices")
    print(f"  {separator}")
    np.save(os.path.join(args.output_dir, "W_matrix.npy"), W)
    np.save(os.path.join(args.output_dir, "H_matrix.npy"), H)
    print(f"    ✓ W_matrix.npy  shape={W.shape}")
    print(f"    ✓ H_matrix.npy  shape={H.shape}")

    print(f"\n  {separator}")
    print(f"  Done. All outputs in {args.output_dir}/")
    print(f"  {separator}\n")


if __name__ == "__main__":
    main()

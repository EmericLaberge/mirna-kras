#!/usr/bin/env python3
# venv: .venv/bin/python3
"""
Transcript-profile miRNA clustering pipeline (OOP):
  DuckDB -> miRNA x transcript binding matrix (25D) -> residual matrix
  -> row-normalized profiles -> k-NN graph (cosine)
  -> Leiden / Spectral (GPU) / Girvan-Newman
  -> parameter sweep with Silhouette + Modularity selection
  -> UMAP visualization + comparison report

Only generates:
  - 3 param_selection plots (one per algorithm)
  - 3 best scatter plots (one per algorithm at optimal param)
  - 1 comparison figure
"""

import warnings
warnings.filterwarnings("ignore")

import os
import argparse
import multiprocessing
from dataclasses import dataclass, field
from time import perf_counter

import duckdb
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

import igraph as ig
import leidenalg
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.neighbors import kneighbors_graph
from umap import UMAP

HAS_GPU = False
try:
    import torch
    if torch.cuda.is_available():
        HAS_GPU = True
except ImportError:
    pass


# ============================================================
# Data structures
# ============================================================

@dataclass
class ClusteringResult:
    """Holds the output of any clustering algorithm."""
    membership: np.ndarray
    n_clusters: int
    silhouette: float
    modularity: float
    runtime: float
    params: dict = field(default_factory=dict)
    algorithm_name: str = ""
    optimal_param_value: float = 0.0
    optimal_param_name: str = ""


@dataclass
class SweepRecord:
    """One row from a parameter sweep."""
    param_value: float
    n_clusters: int
    silhouette: float
    modularity: float
    membership: np.ndarray = field(default_factory=lambda: np.array([]))


# ============================================================
# DataLoader
# ============================================================

class DataLoader:
    """Loads and prepares miRNA-transcript data from DuckDB."""

    def __init__(self, db_path: str, exclude_seedless: bool = True):
        """``exclude_seedless=True`` drops interactions without a canonical seed."""
        self.db_path = db_path
        self.exclude_seedless = exclude_seedless
        self.df: pd.DataFrame | None = None

    def load(self) -> pd.DataFrame:
        """Open DuckDB read-only, fetch the ``results`` table, optionally filter seedless."""
        con = duckdb.connect(self.db_path, read_only=True)
        self.df = con.execute("SELECT * FROM results").fetchdf()
        con.close()

        print(f"  Total interactions : {len(self.df):,}")
        print(f"  Unique miRNAs      : {self.df['mirna_id'].nunique()}")
        print(f"  Unique transcripts : {self.df['transcript_id'].nunique()}")

        if self.exclude_seedless:
            self.df = self.df[self.df["seed_type"] != "seedless"].copy()
            print(f"  After seedless removal : {len(self.df):,} interactions")

        return self.df

    def build_binding_matrix(self, power: float = 1.0) -> tuple[np.ndarray, list[str], list[str]]:
        """Aggregate mean |dG| per (mirna, transcript), pivot to wide, apply power weighting."""
        agg = (
            self.df.groupby(["mirna_id", "transcript_id"])
            .agg(weight=("binding_dG", lambda x: abs(x).mean()))
            .reset_index()
        )

        print(f"  Unique miRNA-transcript pairs: {len(agg):,}")
        print(f"  |dG| range: [{agg['weight'].min():.3f}, {agg['weight'].max():.3f}]")
        print(f"  Power weighting: {power}")

        matrix_df = agg.pivot(index="mirna_id", columns="transcript_id", values="weight").fillna(0)
        matrix = matrix_df.values.copy()

        # Apply power weighting after computing mean |dG|
        matrix = matrix ** power

        mirna_ids = matrix_df.index.tolist()
        transcript_ids = matrix_df.columns.tolist()

        sparsity = (matrix == 0).sum() / matrix.size
        print(f"  Matrix shape: {matrix.shape[0]} x {matrix.shape[1]}")
        print(f"  Sparsity: {sparsity:.1%}")

        return matrix, mirna_ids, transcript_ids

    @staticmethod
    def build_residual_matrix(V_raw: np.ndarray) -> np.ndarray:
        """Same preprocessing as NMF: subtract column means, clip to 0, fill zeros."""
        col_means = V_raw.mean(axis=0)
        V_residual = np.clip(V_raw - col_means, 0, None)
        V_residual[V_residual == 0] = 1e-6
        nonzero_frac = (V_residual > 1e-6).sum() / V_residual.size
        print(f"  Residual matrix non-zero fraction: {nonzero_frac:.1%}")
        return V_residual

    @staticmethod
    def row_normalize(matrix: np.ndarray) -> np.ndarray:
        """Divide each row by its sum (L1 normalize). Rows summing to 0 are left as 0."""
        row_sums = matrix.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        return matrix / row_sums


# ============================================================
# GraphBuilder
# ============================================================

class GraphBuilder:
    """Builds k-NN graphs from binding profiles."""
    def __init__(self, n_neighbors: int = 15, metric: str = "cosine"):
        """``n_neighbors=15`` and cosine similarity are the project-wide defaults."""
        self.n_neighbors = n_neighbors
        self.metric = metric

    def build_knn_graph(self, X: np.ndarray) -> tuple[ig.Graph, np.ndarray]:
        """Build weighted igraph from k-NN on the profiles. Returns (graph, cosine_sim_matrix)."""
        n = X.shape[0]
        knn_sparse = kneighbors_graph(X, n_neighbors=self.n_neighbors, metric=self.metric, include_self=False)
        cos_sim = cosine_similarity(X)

        edges_set: set[tuple[int, int]] = set()
        edge_weights: list[float] = []

        knn_coo = knn_sparse.tocoo()
        for i, j in zip(knn_coo.row, knn_coo.col):
            edge_key = (min(i, j), max(i, j))
            if edge_key not in edges_set:
                edges_set.add(edge_key)
                edge_weights.append(float(cos_sim[i, j]))

        graph = ig.Graph(n=n, edges=list(edges_set), directed=False)
        graph.es["weight"] = edge_weights

        n_edges = graph.ecount()
        density = 2 * n_edges / (n * (n - 1)) if n > 1 else 0
        print(f"  k-NN k={self.n_neighbors}, metric={self.metric}")
        print(f"  Nodes: {n}, Edges: {n_edges:,}, Density: {density:.6f}")

        return graph, cos_sim

    def sparsify(self, graph: ig.Graph, top_k: int = 15) -> ig.Graph:
        """Keep only top-k edges per node by weight."""
        sparsify_edges: set[tuple[int, int]] = set()
        for v in tqdm(range(graph.vcount()), desc="  Sparsifying", ncols=80):
            neighbors = graph.neighbors(v)
            if not neighbors:
                continue
            weights_v = [(n, graph.es[graph.get_eid(v, n)]["weight"]) for n in neighbors]
            weights_v.sort(key=lambda x: x[1], reverse=True)
            for n, _ in weights_v[:top_k]:
                sparsify_edges.add((min(v, n), max(v, n)))

        weights_list = []
        for u, v in sparsify_edges:
            eid = graph.get_eid(u, v)
            weights_list.append(graph.es[eid]["weight"])

        g_sparse = ig.Graph(n=graph.vcount(), edges=list(sparsify_edges), directed=False)
        g_sparse.es["weight"] = weights_list
        print(f"  Sparsified (top-{top_k}/node): {g_sparse.vcount()} nodes, {g_sparse.ecount():,} edges")
        return g_sparse


# ============================================================
# UMAPProjector
# ============================================================

class UMAPProjector:
    """Project 25-D miRNA binding profiles to 2-D for visualization.

    Wraps UMAP with a StandardScaler pre-fit so all features contribute
    equally regardless of their original |binding_dG| magnitude.
    """

    def __init__(self, n_neighbors: int = 30, min_dist: float = 0.1, n_components: int = 2, random_state: int = 42):
        """Store UMAP hyperparameters. Call ``fit_transform(X)`` to project."""
        self.n_neighbors = n_neighbors
        self.min_dist = min_dist
        self.n_components = n_components
        self.random_state = random_state

    def fit_transform(self, X: np.ndarray) -> np.ndarray:
        """Standardize ``X`` then project via UMAP. Returns shape (n, n_components)."""
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        model = UMAP(
            n_components=self.n_components,
            n_neighbors=self.n_neighbors,
            min_dist=self.min_dist,
            metric="euclidean",
            random_state=self.random_state,
        )
        return model.fit_transform(X_scaled)


# ============================================================
# Plotter
# ============================================================

class Plotter:
    """Static plotting utilities."""

    @staticmethod
    def param_selection(sweep: list[SweepRecord], param_name: str, title: str,
                        filepath: str, optimal_value: float) -> None:
        """Dual-panel plot: Silhouette + Modularity vs parameter. Like NMF rank_selection.png."""
        params = [r.param_value for r in sweep]
        sils = [r.silhouette for r in sweep]
        mods = [r.modularity for r in sweep]

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        ax1.plot(params, sils, "o-", color="#2c3e50", linewidth=2, markersize=8)
        ax1.axvline(optimal_value, color="red", linestyle="--", alpha=0.7,
                     label=f"optimal {param_name}={optimal_value}")
        ax1.set_xlabel(param_name, fontsize=12)
        ax1.set_ylabel("Silhouette coefficient", fontsize=12)
        ax1.set_title("Silhouette vs " + param_name, fontsize=13)
        ax1.legend(fontsize=11)

        ax2.plot(params, mods, "s-", color="#8e44ad", linewidth=2, markersize=8)
        ax2.axvline(optimal_value, color="red", linestyle="--", alpha=0.7,
                     label=f"optimal {param_name}={optimal_value}")
        ax2.set_xlabel(param_name, fontsize=12)
        ax2.set_ylabel("Modularity", fontsize=12)
        ax2.set_title("Modularity vs " + param_name, fontsize=13)
        ax2.legend(fontsize=11)

        fig.tight_layout()
        fig.savefig(filepath, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {filepath}")

    @staticmethod
    def scatter(X_umap: np.ndarray, membership: np.ndarray, mirna_ids: list[str],
                title: str, filepath: str) -> None:
        """Scatter UMAP projection colored by cluster membership. One figure per algorithm."""
        membership = np.array(membership)
        n_clusters = len(set(membership))
        cmap = plt.cm.get_cmap("tab20", max(n_clusters, 1))

        fig, ax = plt.subplots(figsize=(14, 10))
        for c in sorted(set(membership)):
            mask = membership == c
            color = cmap(c % 20)
            ax.scatter(
                X_umap[mask, 0], X_umap[mask, 1],
                c=[color], s=25, alpha=0.7,
                label=f"Cluster {c} ({mask.sum()})",
            )

        ax.set_title(title, fontsize=16, fontweight="bold")
        ax.set_xlabel("UMAP 1", fontsize=12)
        ax.set_ylabel("UMAP 2", fontsize=12)
        ax.legend(markerscale=2, fontsize=9, loc="best", framealpha=0.9)
        ax.set_facecolor("#f8f8f8")
        fig.patch.set_facecolor("white")
        plt.tight_layout()
        fig.savefig(filepath, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {filepath}")

    def comparison(X_umap: np.ndarray, results: list[ClusteringResult],
                   mirna_ids: list[str], filepath: str) -> None:
        """Side-by-side UMAP scatters, one panel per algorithm, for cross-comparison."""
        fig, axes = plt.subplots(1, len(results), figsize=(7 * len(results), 7))
        if len(results) == 1:
            axes = [axes]

        for ax, res in zip(axes, results):
            membership = np.array(res.membership)
            n_c = len(set(membership))
            cmap = plt.cm.get_cmap("tab20", max(n_c, 1))
            for c in sorted(set(membership)):
                mask = membership == c
                ax.scatter(
                    X_umap[mask, 0], X_umap[mask, 1],
                    c=[cmap(c % 20)], s=15, alpha=0.7,
                    label=f"C{c} ({mask.sum()})",
                )
            ax.set_title(res.algorithm_name, fontsize=13, fontweight="bold")
            ax.set_xlabel("UMAP 1", fontsize=10)
            ax.set_ylabel("UMAP 2", fontsize=10)
            ax.legend(fontsize=7, loc="best", framealpha=0.9, ncol=2)
            ax.set_facecolor("#f8f8f8")

        fig.patch.set_facecolor("white")
        plt.tight_layout()
        fig.savefig(filepath, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {filepath}")


# ============================================================
# LeidenClustering
# ============================================================

class LeidenClustering:
    """Leiden community detection with a sweep over the resolution parameter.

    For each resolution gamma in ``resolutions``, run Leiden on the k-NN
    graph and record silhouette + modularity. The best partition is the
    one with the highest silhouette.
    """

    def __init__(self, resolutions: list[float] | None = None):
        """Defaults to sweeping gamma from 0.1 to 2.0 in 0.1 steps."""
        if resolutions is None:
            self.resolutions = [round(0.1 * i, 1) for i in range(1, 21)]  # 0.1 to 2.0
        else:
            self.resolutions = resolutions

    def run(self, graph: ig.Graph, X: np.ndarray, X_umap: np.ndarray,
            mirna_ids: list[str], output_dir: str) -> tuple[ClusteringResult, list[SweepRecord]]:
        """Sweep resolutions, select by highest Silhouette. Returns (best_result, sweep_records)."""
        print(f"  Leiden clustering (CPU) — sweeping {len(self.resolutions)} resolutions")
        t0 = perf_counter()
        sweep: list[SweepRecord] = []

        for res in self.resolutions:
            try:
                p = leidenalg.find_partition(
                    graph,
                    leidenalg.RBConfigurationVertexPartition,
                    resolution_parameter=res,
                    weights="weight",
                    seed=42,
                )
                mem = np.array(p.membership)
                n = len(set(mem))
                if n < 2:
                    continue
                mod = graph.modularity(mem, weights="weight")
                sil = silhouette_score(X, mem)
                print(f"    res={res:.1f}: {n} clusters, sil={sil:.4f}, mod={mod:.4f}")
                sweep.append(SweepRecord(
                    param_value=res, n_clusters=n,
                    silhouette=sil, modularity=mod,
                    membership=mem.copy(),
                ))
            except Exception:
                continue

        if not sweep:
            raise RuntimeError("Leiden: no valid partition found")

        # Select by highest Silhouette
        best = max(sweep, key=lambda r: r.silhouette)
        runtime = perf_counter() - t0

        print(f"  Best Leiden: res={best.param_value:.1f}, {best.n_clusters} clusters, "
              f"sil={best.silhouette:.4f}, mod={best.modularity:.4f}, time={runtime:.2f}s")

        # Save best figure
        Plotter.scatter(X_umap, best.membership, mirna_ids,
                        f"miRNA Clustering — Leiden (γ={best.param_value:.1f}, "
                        f"{best.n_clusters} clusters)\nSil={best.silhouette:.3f}, Mod={best.modularity:.3f}",
                        f"{output_dir}/best_leiden.png")

        # Save param selection plot
        Plotter.param_selection(sweep, "Resolution (γ)",
                                "Leiden Parameter Selection",
                                f"{output_dir}/leiden_param_selection.png",
                                best.param_value)

        result = ClusteringResult(
            membership=best.membership,
            n_clusters=best.n_clusters,
            silhouette=best.silhouette,
            modularity=best.modularity,
            runtime=runtime,
            params={"resolutions": self.resolutions, "optimal_resolution": best.param_value},
            algorithm_name=f"Leiden (γ={best.param_value:.1f}, {best.n_clusters} clusters)",
            optimal_param_value=best.param_value,
            optimal_param_name="γ",
        )
        return result, sweep


# ============================================================
# SpectralClustering
# ============================================================

class SpectralClustering:
    """Spectral clustering with a sweep over the number of eigenvectors k.

    Uses the normalized graph Laplacian. Selects k by highest silhouette.
    Uses GPU (PyTorch) when available and ``use_gpu=True``, else falls
    back to scikit-learn.
    """

    def __init__(self, k_range: list[int] | None = None, use_gpu: bool = True):
        """``k_range`` defaults to 2..10. GPU auto-disabled if CUDA is absent."""
        self.k_range = k_range or list(range(2, 11))  # 2 to 10
        self.use_gpu = use_gpu and HAS_GPU

    def run(self, adj_matrix: np.ndarray, X: np.ndarray, graph: ig.Graph,
            X_umap: np.ndarray, mirna_ids: list[str], output_dir: str
            ) -> tuple[ClusteringResult, list[SweepRecord]]:
        """Sweep k in ``self.k_range``, select by highest silhouette, render figures.

        Args:
            adj_matrix: Cosine similarity matrix (n, n).
            X: Row-normalized feature matrix used for silhouette scoring.
            graph: igraph object used for modularity scoring.
            X_umap: 2-D UMAP embedding for the scatter plot.
            mirna_ids: miRNA labels for the scatter plot legend.
            output_dir: Where to write ``best_spectral.png`` and ``spectral_param_selection.png``.

        Returns:
            Tuple of (best ClusteringResult, full sweep records).
        """
        backend = "GPU (CUDA)" if self.use_gpu else "CPU (sklearn)"
        print(f"  Spectral clustering ({backend}) — sweeping k={self.k_range}")
        t0 = perf_counter()
        sweep: list[SweepRecord] = []

        n = adj_matrix.shape[0]
        W = (adj_matrix + adj_matrix.T) / 2.0
        np.fill_diagonal(W, 0.0)

        if self.use_gpu:
            W_tensor = torch.tensor(W, dtype=torch.float32).to("cuda")
            print(f"  GPU memory: {torch.cuda.memory_allocated() / 1e6:.1f} MB")

            for k in self.k_range:
                t_k = perf_counter()
                epsilon = 1e-6
                W_gpu = W_tensor + epsilon * torch.eye(n, device="cuda")

                d = W_gpu.sum(dim=1)
                D_inv_sqrt = torch.diag(1.0 / torch.sqrt(d))
                I = torch.eye(n, device="cuda")
                L = I - D_inv_sqrt @ W_gpu @ D_inv_sqrt
                L = (L + L.T) / 2.0

                _, eigenvectors = torch.linalg.eigh(L)
                U = eigenvectors[:, 1:k + 1].cpu().numpy()

                row_norms = np.linalg.norm(U, axis=1, keepdims=True)
                row_norms[row_norms == 0] = 1.0
                U_norm = U / row_norms

                km = KMeans(n_clusters=k, random_state=42, n_init=20)
                labels_k = km.fit_predict(U_norm)

                sil = silhouette_score(X, labels_k)
                mod = graph.modularity(labels_k, weights="weight")
                dt = perf_counter() - t_k
                print(f"    k={k}: sil={sil:.4f}, mod={mod:.4f}, time={dt:.2f}s")

                sweep.append(SweepRecord(
                    param_value=float(k), n_clusters=k,
                    silhouette=sil, modularity=mod,
                    membership=labels_k.copy(),
                ))

            torch.cuda.empty_cache()
        else:
            from sklearn.cluster import SpectralClustering as SkSpectral

            for k in self.k_range:
                t_k = perf_counter()
                sc = SkSpectral(n_clusters=k, affinity="precomputed", random_state=42, n_init=20)
                labels_k = sc.fit_predict(W)

                sil = silhouette_score(X, labels_k)
                mod = graph.modularity(labels_k, weights="weight")
                dt = perf_counter() - t_k
                print(f"    k={k}: sil={sil:.4f}, mod={mod:.4f}, time={dt:.2f}s")

                sweep.append(SweepRecord(
                    param_value=float(k), n_clusters=k,
                    silhouette=sil, modularity=mod,
                    membership=labels_k.copy(),
                ))

        if not sweep:
            raise RuntimeError("Spectral: no valid k found")

        # Select by highest Silhouette
        best = max(sweep, key=lambda r: r.silhouette)
        runtime = perf_counter() - t0

        print(f"  Best Spectral: k={int(best.param_value)}, sil={best.silhouette:.4f}, "
              f"mod={best.modularity:.4f}, time={runtime:.2f}s")

        Plotter.scatter(X_umap, best.membership, mirna_ids,
                        f"miRNA Clustering — Spectral (k={int(best.param_value)})\n"
                        f"Sil={best.silhouette:.3f}, Mod={best.modularity:.3f}",
                        f"{output_dir}/best_spectral.png")

        Plotter.param_selection(sweep, "k (clusters)",
                                "Spectral Parameter Selection",
                                f"{output_dir}/spectral_param_selection.png",
                                best.param_value)

        result = ClusteringResult(
            membership=best.membership,
            n_clusters=best.n_clusters,
            silhouette=best.silhouette,
            modularity=best.modularity,
            runtime=runtime,
            params={"best_k": int(best.param_value), "k_range": self.k_range},
            algorithm_name=f"Spectral (k={int(best.param_value)})",
            optimal_param_value=best.param_value,
            optimal_param_name="k",
        )
        return result, sweep


# ============================================================
# GirvanNewmanClustering
# ============================================================

class GirvanNewmanClustering:
    """Girvan–Newman edge-betweenness community detection.

    Sparsifies the k-NN graph (top-k edges per node) before running, since GN
    is O(m²) and the dense graph would never terminate. Even with sparsification
    this method usually collapses to a single cluster on our data.
    """

    def __init__(self, target_range: list[int] | None = None,
                 sparsify_top_k: int = 15, max_time_per_target: int = 120):
        """``target_range`` defaults to 2..10. Times out per target to avoid hangs."""
        self.target_range = target_range or list(range(2, 11))  # 2 to 10
        self.sparsify_top_k = sparsify_top_k
        self.max_time_per_target = max_time_per_target

    def _run_single(self, graph: ig.Graph, target: int, max_time: int) -> np.ndarray | None:
        """Run GN for a single target_communities value. Returns membership or None."""
        t0 = perf_counter()

        builder = GraphBuilder()
        g_sparse = builder.sparsify(graph, top_k=self.sparsify_top_k)

        # Single-pass edge betweenness + weight-normalized ranking
        eb_scores = g_sparse.edge_betweenness(weights="weight")
        edge_weights = g_sparse.es["weight"]
        max_weight = max(edge_weights) if edge_weights else 1.0
        eb_normalized = [eb / max(w / max_weight, 0.01) for eb, w in zip(eb_scores, edge_weights)]

        edge_order = sorted(range(len(eb_normalized)), key=lambda i: eb_normalized[i], reverse=True)

        # Batch removal until target communities reached
        g_gn = g_sparse.copy()
        batch_size = max(1, len(edge_order) // 200)
        max_remove = len(edge_order) // 2

        pbar = tqdm(total=max_remove, desc=f"  GN target={target}", ncols=80, unit="edges")
        removed_total = 0
        n_components = 1

        for batch_start in range(0, min(len(edge_order), max_remove), batch_size):
            if perf_counter() - t0 > max_time:
                pbar.close()
                break

            batch_end = min(batch_start + batch_size, len(edge_order))
            batch_indices = edge_order[batch_start:batch_end]

            removed_batch = 0
            for idx in batch_indices:
                try:
                    e = g_gn.es[idx]
                    src, tgt = e.tuple
                    if g_gn.are_connected(src, tgt):
                        g_gn.delete_edges(idx)
                        removed_batch += 1
                        removed_total += 1
                except (IndexError, ValueError):
                    continue

            pbar.update(removed_batch)

            components = g_gn.connected_components()
            n_components = len(set(components.membership))
            elapsed = perf_counter() - t0
            pbar.set_postfix({"comms": n_components, "elapsed": f"{elapsed:.0f}s"})

            if n_components >= target:
                pbar.close()
                break
        else:
            pbar.close()

        if n_components < 2:
            return None

        gn_membership = np.array(components.membership)

        # Merge tiny components (<5 miRNAs) into largest cluster
        min_size = 5
        large_clusters = {c for c in set(gn_membership) if (gn_membership == c).sum() >= min_size}
        if len(large_clusters) < n_components:
            tiny_mask = np.array([(gn_membership == c).sum() < min_size for c in gn_membership])
            if tiny_mask.any() and len(large_clusters) > 0:
                largest = max(set(gn_membership), key=lambda c: (gn_membership == c).sum())
                gn_membership[tiny_mask] = largest
                unique_labels = sorted(set(gn_membership))
                label_map = {old: new for new, old in enumerate(unique_labels)}
                gn_membership = np.array([label_map[l] for l in gn_membership])

        return gn_membership

    def run(self, graph: ig.Graph, X: np.ndarray, X_umap: np.ndarray,
            mirna_ids: list[str], output_dir: str) -> tuple[ClusteringResult, list[SweepRecord]]:
        """Sweep target_communities, select by highest Silhouette."""
        print(f"  Girvan-Newman — sweeping target_communities={self.target_range}")
        t0 = perf_counter()
        sweep: list[SweepRecord] = []

        for target in self.target_range:
            print(f"\n  --- GN target={target} ---")
            mem = self._run_single(graph, target, self.max_time_per_target)
            if mem is None:
                print(f"    target={target}: failed to partition")
                continue

            n_clusters = len(set(mem))
            sil = silhouette_score(X, mem) if n_clusters > 1 else 0.0
            mod = graph.modularity(mem, weights="weight")
            runtime_target = perf_counter() - t0
            print(f"    target={target}: {n_clusters} clusters, sil={sil:.4f}, mod={mod:.4f}")

            sweep.append(SweepRecord(
                param_value=float(target), n_clusters=n_clusters,
                silhouette=sil, modularity=mod,
                membership=mem.copy(),
            ))

        runtime = perf_counter() - t0

        if not sweep:
            # Fallback: single cluster
            print("  Girvan-Newman: all targets failed, returning single cluster")
            mem = np.zeros(graph.vcount(), dtype=int)
            result = ClusteringResult(
                membership=mem, n_clusters=1,
                silhouette=0.0, modularity=0.0,
                runtime=runtime, params={"targets": self.target_range},
                algorithm_name="Girvan-Newman (failed)",
            )
            return result, []

        # Select by highest Silhouette
        best = max(sweep, key=lambda r: r.silhouette)

        print(f"\n  Best GN: target={int(best.param_value)}, {best.n_clusters} clusters, "
              f"sil={best.silhouette:.4f}, mod={best.modularity:.4f}, total time={runtime:.1f}s")

        Plotter.scatter(X_umap, best.membership, mirna_ids,
                        f"miRNA Clustering — Girvan-Newman (target={int(best.param_value)}, "
                        f"{best.n_clusters} clusters)\nSil={best.silhouette:.3f}, Mod={best.modularity:.3f}",
                        f"{output_dir}/best_girvan_newman.png")

        if len(sweep) > 1:
            Plotter.param_selection(sweep, "Target communities",
                                    "Girvan-Newman Parameter Selection",
                                    f"{output_dir}/gn_param_selection.png",
                                    best.param_value)

        result = ClusteringResult(
            membership=best.membership,
            n_clusters=best.n_clusters,
            silhouette=best.silhouette,
            modularity=best.modularity,
            runtime=runtime,
            params={"optimal_target": int(best.param_value), "targets": self.target_range},
            algorithm_name=f"Girvan-Newman (target={int(best.param_value)}, {best.n_clusters} clusters)",
            optimal_param_value=best.param_value,
            optimal_param_name="target",
        )
        return result, sweep


# ============================================================
# PipelineReport
# ============================================================

class PipelineReport:
    """Prints comparison table and ARI between all results."""

    @staticmethod
    def compare(results: list[ClusteringResult], mirna_ids: list[str], df: pd.DataFrame) -> None:
        """Print a side-by-side metrics table, pairwise ARI, and top-5 transcripts per cluster."""
        header = f"  {'Algorithm':<40} {'Clusters':>10} {'Silhouette':>12} {'Modularity':>12} {'Time (s)':>10}"
        print(header)
        print("  " + "-" * (len(header) - 2))
        for r in results:
            print(f"  {r.algorithm_name:<40} {r.n_clusters:>10} {r.silhouette:>12.4f} "
                  f"{r.modularity:>12.4f} {r.runtime:>10.2f}")
        print()

        print("  Adjusted Rand Index (pairwise agreement):")
        for i in range(len(results)):
            for j in range(i + 1, len(results)):
                ari = adjusted_rand_score(results[i].membership, results[j].membership)
                print(f"    {results[i].algorithm_name} vs {results[j].algorithm_name}: {ari:.4f}")
        print()

        for r in results:
            if r.n_clusters < 2:
                continue
            print(f"  {r.algorithm_name} — top 5 transcripts per cluster:")
            mirna_cluster_df = pd.DataFrame({"mirna_id": mirna_ids, "cluster": r.membership})
            df_clustered = df.merge(mirna_cluster_df, on="mirna_id", how="inner")

            for c in sorted(set(r.membership)):
                cluster_df = df_clustered[df_clustered["cluster"] == c]
                top_trans = cluster_df["transcript_id"].value_counts().head(5)
                n_mirnas_c = cluster_df["mirna_id"].nunique()
                print(f"\n    Cluster {c} ({n_mirnas_c} miRNAs, {len(cluster_df):,} interactions):")
                print(f"    {'Rank':<6} {'Transcript':<20} {'Count':>10} {'%':>8}")
                print(f"    {'-' * 6} {'-' * 20} {'-' * 10} {'-' * 8}")
                for rank, (tid, cnt) in enumerate(top_trans.items(), 1):
                    pct = 100 * cnt / len(cluster_df)
                    print(f"    {rank:<6} {tid:<20} {cnt:>10,} {pct:>7.1f}%")
            print()


# ============================================================
# main
# ============================================================

def main():
    """CLI entry: load DuckDB -> residual -> k-NN -> Leiden/Spectral/GN sweep -> figures."""
    parser = argparse.ArgumentParser(description="miRNA transcript-profile clustering pipeline")
    parser.add_argument("--power", type=float, default=1.0,
                        help="Power weighting for |dG| (1=linear, 2=quadratic)")
    parser.add_argument("--k-neighbors", type=int, default=15,
                        help="k for k-NN graph construction")
    parser.add_argument("--output-dir", type=str, default="figures",
                        help="Output directory for figures")
    parser.add_argument("--n-components", type=int, default=2,
                        help="UMAP dimensions (2 or 3)")
    args = parser.parse_args()

    n_cores = multiprocessing.cpu_count()
    print(f"  CPU cores: {n_cores}, GPU: {HAS_GPU}")
    print(f"  Power: {args.power}, k-neighbors: {args.k_neighbors}, dims: {args.n_components}")
    print()

    plt.style.use("seaborn-v0_8-whitegrid")
    os.makedirs(args.output_dir, exist_ok=True)

    # 1. Load data
    print("=" * 60)
    print("1. Loading data from DuckDB")
    print("=" * 60)
    loader = DataLoader("rimap_results.duckdb", exclude_seedless=True)
    df = loader.load()
    X_raw, mirna_ids, transcript_ids = loader.build_binding_matrix(power=args.power)
    print()

    # 2. Residual matrix (same as NMF)
    print("=" * 60)
    print("2. Computing residual matrix (same as NMF)")
    print("=" * 60)
    V_residual = DataLoader.build_residual_matrix(X_raw)
    X_profiles = DataLoader.row_normalize(V_residual)
    print()

    # 3. UMAP on residual matrix
    print("=" * 60)
    print("3. Computing UMAP 2D embedding (on residual matrix)")
    print("=" * 60)
    t0 = perf_counter()
    projector = UMAPProjector(n_neighbors=30, min_dist=0.1, n_components=args.n_components, random_state=42)
    X_umap = projector.fit_transform(V_residual)
    print(f"  UMAP shape: {X_umap.shape}, time: {perf_counter() - t0:.2f}s")
    print()

    # 4. Graph on residual profiles
    print("=" * 60)
    print("4. Building k-NN graph (on residual profiles)")
    print("=" * 60)
    t0 = perf_counter()
    builder = GraphBuilder(n_neighbors=args.k_neighbors, metric="cosine")
    graph, cos_sim = builder.build_knn_graph(X_profiles)
    print(f"  Graph build time: {perf_counter() - t0:.2f}s")
    print()

    # 5. Leiden sweep
    print("=" * 60)
    print("5. Leiden Clustering — parameter sweep")
    print("=" * 60)
    leiden = LeidenClustering(resolutions=[round(0.1 * i, 1) for i in range(1, 21)])
    r_leiden, sweep_leiden = leiden.run(graph, V_residual, X_umap, mirna_ids, args.output_dir)
    print()

    # 6. Spectral sweep
    print("=" * 60)
    print("6. Spectral Clustering — parameter sweep")
    print("=" * 60)
    spectral = SpectralClustering(k_range=list(range(2, 11)), use_gpu=True)
    r_spectral, sweep_spectral = spectral.run(cos_sim, V_residual, graph, X_umap, mirna_ids, args.output_dir)
    print()

    # 7. Girvan-Newman sweep
    print("=" * 60)
    print("7. Girvan-Newman — parameter sweep")
    print("=" * 60)
    gn = GirvanNewmanClustering(target_range=list(range(2, 11)), sparsify_top_k=15, max_time_per_target=120)
    r_gn, sweep_gn = gn.run(graph, V_residual, X_umap, mirna_ids, args.output_dir)
    print()

    # 8. Compare
    print("=" * 60)
    print("8. Algorithm comparison")
    print("=" * 60)
    results = [r_leiden, r_spectral, r_gn]
    PipelineReport.compare(results, mirna_ids, df)
    Plotter.comparison(X_umap, results, mirna_ids, f"{args.output_dir}/comparison.png")
    print()

    print("=" * 60)
    print("DONE")
    print("=" * 60)


if __name__ == "__main__":
    main()

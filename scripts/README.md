# `scripts/` — Python pipeline

All scripts run with `.venv/bin/python3` (Python 3.14). Most read `rimap_results.duckdb` from the repo root.

## Main pipeline (execution order)

### Stage 0 — Data

| Script | Role |
|---|---|
| `gene_ids.py` | List of 25 transcript Ensembl IDs (constants, no execution). |
| `fetch_gencode.py` | Fetches transcript sequences via the Ensembl REST API. |
| `extract_transcripts.py` | Alternative: extracts transcripts from a local GENCODE archive. |

### Stage 1 — Binding generation (long, hours)

| Script | Role |
|---|---|
| `rimap_pipeline.py` | Computes affinities via the RIMap-RISC web API → `rimap_results.duckdb`. |
| `rimap_pipeline_local.py` | Local alternative — Java MiRScanS + MC-Flashfold (faster, but requires install). |

### Stage 2 — Analyses (on the `.duckdb`)

| Script | Output | Role |
|---|---|---|
| `transcript_clustering.py` | `figures/best_*.png`, `param_selection_*.png`, `comparison.png` | Graph-based clustering (Leiden, Spectral, Girvan–Newman) with hyperparameter sweep. |
| `nmf_analysis.py` | `figures_nmf/{rank_selection,k2..k15}/`, `W_matrix.npy`, `H_matrix.npy` | NMF with rank selection by cophenetic correlation (20 seeds per k). |
| `chord_diagram.py` | `figures/chord_diagram.png` | Gene–gene co-targeting chord diagram. |
| `heatmap_by_cluster.py` | `figures/heatmap_sorted_nmf.png` | Heatmap sorted by dominant NMF signature. |
| `bipartite_graph.py` | `figures/bipartite_graph*.png` | Bipartite miRNA–transcript graph (raw and residual versions). |
| `bipartite_by_cluster.py` | `figures/bipartite_by_cluster.png` | Bipartite graph filtered by dominant NMF signature. |
| `bipartite_top_heatmap.py` | `figures/top_mirnas_heatmap.png` | Top-N miRNAs per gene, clustered heatmap. |
| `mirna_interactive_plot.py` | `viz/index.html` (if the folder exists) | Interactive Plotly UMAP + HDBSCAN plot. |

### Utilities (not part of the pipeline)

| Script | Role |
|---|---|
| `test_curl.py` | One-shot smoke test of the RIMap-RISC API (single miRNA–transcript pair). |

## Typical execution

```bash
# Regenerate all figures + heatmaps from the existing .duckdb
.venv/bin/python3 scripts/transcript_clustering.py
.venv/bin/python3 scripts/nmf_analysis.py --output-dir figures_nmf --n-seeds 20
.venv/bin/python3 scripts/chord_diagram.py
.venv/bin/python3 scripts/heatmap_by_cluster.py
.venv/bin/python3 scripts/bipartite_by_cluster.py
.venv/bin/python3 scripts/bipartite_top_heatmap.py
```

## Dependencies

See [`requirements.txt`](../requirements.txt) at the root (create one if it doesn't exist yet, see the *Setup* section of the [README](../README.md)). Main dependencies: `duckdb`, `numpy`, `pandas`, `scikit-learn`, `scipy`, `igraph`, `leidenalg`, `matplotlib`, `seaborn`, `plotly`, `umap-learn`, `hdbscan`, `tqdm`, `requests`.
# Getting started

## Prerequisites

- Python 3.12+
- Java 17+ (only for `rimap_pipeline_local.py`)
- ~5 GB disk space for the DuckDB database

## Installation

```bash
# Virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Pipeline dependencies
pip install duckdb numpy pandas scikit-learn scipy \
            igraph leidenalg matplotlib seaborn \
            plotly umap-learn hdbscan tqdm requests

# Documentation site dependencies (optional)
pip install mkdocs mkdocs-material mkdocstrings mkdocstrings-python
```

## Pipeline

The pipeline has **three stages**. All commands run from the repo root.

### Stage 0 — Data

```bash
# Ensembl IDs of the 25 transcripts (constants, no execution)
python scripts/gene_ids.py

# One of the two variants:
python scripts/fetch_gencode.py        # via the Ensembl REST API
# or
python scripts/extract_transcripts.py  # from a local GENCODE archive (gencode.v46.transcripts.fa.gz)
```

→ Produces `data/fasta/*.fasta` (25 files, a few KB each).

### Stage 1 — Interaction computation (~hours)

```bash
# Web API variant (recommended) — 4 parallel workers, checkpoint
python scripts/rimap_pipeline.py

# Local variant — Java MiRScanS + MC-Flashfold (faster if installed)
python scripts/rimap_pipeline_local.py
```

→ Produces `rimap_results.duckdb` (~214 MB, gitignored).

### Stage 2 — Analyses (seconds to minutes)

```bash
# Graph clustering
python scripts/transcript_clustering.py
# NMF + rank selection (k=2..15, 20 seeds per k)
python scripts/nmf_analysis.py --output-dir figures_nmf --n-seeds 20
# Visualizations
python scripts/chord_diagram.py
python scripts/heatmap_by_cluster.py
python scripts/bipartite_by_cluster.py
python scripts/bipartite_graph.py
python scripts/bipartite_top_heatmap.py
python scripts/mirna_interactive_plot.py
```

→ Produces `figures/*.png` and `figures_nmf/{k2..k15}/*.png`.

## Documentation site

To build this site locally:

```bash
mkdocs serve              # dev server, hot reload
mkdocs build              # static site to ./site
mkdocs gh-deploy          # GitHub Pages deploy (if configured)
```

## Conventions

See [Architecture — Conventions](architecture.md#conventions-and-invariants).

## Troubleshooting

- **DB not found**: the repo doesn't include `rimap_results.duckdb` (gitignored). Regenerate via Stage 1.
- **ModuleNotFoundError**: did you activate the venv? `source .venv/bin/activate`
- **CUDA not available**: `transcript_clustering.py` detects this automatically and falls back to CPU.
- **RIMap-RISC API down**: the pipeline has a text checkpoint (`pipeline_checkpoint.txt`) — re-running resumes where it left off.
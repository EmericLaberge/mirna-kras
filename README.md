# miRNA–KRas: binding-profile classification



**License:** [GPL-3.0](LICENSE)

## What this project does

We predict miRNA–transcript binding affinities for **2,516 mature miRNAs × 25 KRas-network transcripts** using [RIMap-RISC](https://rimap-risc.api.major.iric.ca), then cluster the miRNAs into **functional signatures** via non-negative matrix factorization (NMF) and graph-based clustering.

## Main result

NMF at **k = 15** extracts 15 latent signatures; the top 5 (69.5% of miRNAs) separate the PI3K regulators (PIK3CA, PIK3CB) from the direct RAS effectors (KRAS, NRAS, HRAS). Only **3.5% of miRNAs are "mixed"** (membership entropy > 0.8) — most miRNAs slot cleanly into one signature.

Hard clustering (Leiden, Spectral, Girvan–Newman) on the same residual matrix converges to a coarse bipartition and fails to isolate sub-structures — confirming the signal is better captured with soft membership.

## Navigating the repo

| You're looking for… | Go to |
|---|---|
| **Quick-start** (install + run commands) | [`docs/getting-started.md`](docs/getting-started.md) — rendered at the [mkdocs site](docs/index.md) |
| **Architecture** (data flow, DuckDB schema, modules) | [`docs/architecture.md`](docs/architecture.md) |
| **Data model** (DuckDB schema, 4 tables) | [`docs/data-model.md`](docs/data-model.md) |
| **API reference** (every script, auto-generated from docstrings) | [`docs/reference/index.md`](docs/reference/index.md) |
| **Extension points** (add a transcript, a clustering method, …) | [`docs/extension.md`](docs/extension.md) |
| **NMF method** (algorithm, rank selection, interpretation) | [`doc/nmf_analysis.md`](doc/nmf_analysis.md) |
| **Indices** | [`doc/README.md`](doc/README.md) · [`scripts/README.md`](scripts/README.md) · [`data/README.md`](data/README.md) |
| **Figures** | `figures/` (graph clustering) · `figures_nmf/` (NMF) |
| **Database** | `rimap_results.duckdb` — gitignored, ~214 MB; regenerate via the pipeline |

## Pipeline — overview

```
data/mirnas.fa + data/fasta/*.fasta
            │
            ▼
  [1] gene_ids.py                  Ensembl IDs for the 25 transcripts
  [2] extract_transcripts.py       FASTA sequences (from local archive)
       or fetch_gencode.py         …or via Ensembl REST API
  [3] rimap_pipeline.py            Compute bindings via RIMap-RISC API
       or rimap_pipeline_local.py  …or locally (MiRScanS + MC-Flashfold)
            │
            ▼
       rimap_results.duckdb  (2,516 × 25, mean |binding_dG|)
            │
            ├─► [4] transcript_clustering.py   Leiden / Spectral / GN
            ├─► [5] nmf_analysis.py            NMF + rank selection
            ├─► [6] chord_diagram.py           Gene–gene co-targeting
            ├─► [7] heatmap_by_cluster.py      Heatmap sorted by signature
            └─► [8] mirna_interactive_plot.py  UMAP + HDBSCAN interactive
```

See [`scripts/README.md`](scripts/README.md) for the full list and execution order.

## Clustering methods tested

| Method | Type | Best result | Metric |
|---|---|---|---|
| **NMF** (Multiplicative Updates + NNDSVD) | Soft membership | k = 15, RAS/PI3K signatures | Cophenetic ≈ 0.9 on the consensus matrix |
| Leiden (k-NN cosine) | Hard membership | γ = 0.2, 2 clusters | Silhouette = 0.06, Modularity = 0.34 |
| Spectral (k-NN cosine) | Hard membership | k = 2 | Silhouette = 0.12, Modularity = 0.37 |
| Girvan–Newman | Hard membership | 1 cluster — fails | Graph too dense |

## Key preprocessing — residual matrix

Raw binding profiles are nearly identical (CV = 0.06) — a typical miRNA binds every transcript with roughly the same strength. To extract the **differential signal**:

```
V_res[i,j] = max(0, V[i,j] − mean_col(V)[j])
```

Only positive deviations from the per-transcript mean are kept. → 58.7% non-zero entries, compatible with NMF.

## Quick reproduction

Results under `figures/` and `figures_nmf/` are already generated. To regenerate from `rimap_results.duckdb`:

```bash
.venv/bin/python3 scripts/transcript_clustering.py
.venv/bin/python3 scripts/nmf_analysis.py --output-dir figures_nmf --n-seeds 20
.venv/bin/python3 scripts/chord_diagram.py
.venv/bin/python3 scripts/heatmap_by_cluster.py
.venv/bin/python3 scripts/bipartite_by_cluster.py
.venv/bin/python3 scripts/bipartite_top_heatmap.py
.venv/bin/python3 scripts/mirna_interactive_plot.py
```

To regenerate the DB from scratch (long — API calls):

```bash
.venv/bin/python3 scripts/gene_ids.py > /dev/null
.venv/bin/python3 scripts/fetch_gencode.py         # or extract_transcripts.py
.venv/bin/python3 scripts/rimap_pipeline.py        # hours
```

## Repo layout

```
.
├── README.md                  This file
├── LICENSE                    GPL-3.0
├── requirements.txt           Python deps (see Setup in docs)
├── mkdocs.yml                 mkdocs config (Google-style docstrings)
├── rimap_results.duckdb       Interactions DB (gitignored, 214 MB)
├── data/
│   ├── mirnas.fa              2,588 mature miRNAs (MIMAT IDs, miRBase)
│   └── fasta/                 25 transcripts (ENST IDs, GENCODE)
├── scripts/                   Python pipeline (see scripts/README.md)
├── figures/                   Graph-clustering results
├── figures_nmf/               NMF results (rank_selection, k2..k15)
├── doc/
│   ├── README.md              Index of doc/
│   ├── nmf_analysis.md        NMF method details
│   └── presentation/          Beamer slides
├── kami-out/                  Stylized dev-doc (HTML + PDF)
└── docs/                      mkdocs source (built site in ./site)
```

## Documentation

| Format | How to read |
|---|---|
| **mkdocs site** (recommended) | `mkdocs serve` → open `http://127.0.0.1:8000`. Build static via `mkdocs build` (output in `./site`). |
| **Plain Markdown** | Each `*.md` is self-contained. The `docs/` folder mirrors the mkdocs nav. |
| **Kami-styled PDF** | `kami-out/architecture.pdf` — parchment design system, dev-focused |

## Dependencies (Python)

`duckdb`, `numpy`, `pandas`, `scikit-learn`, `scipy`, `igraph`, `leidenalg`, `matplotlib`, `seaborn`, `plotly`, `umap-learn`, `hdbscan`, `tqdm`, `requests`. Optional (for the docs site): `mkdocs`, `mkdocs-material`, `mkdocstrings`, `mkdocstrings-python`.

## Keywords

miRNA, KRas, NMF, Leiden clustering, binding affinity, RIMap-RISC, RAS-MAPK, PI3K-AKT, non-negative factorization, cancer.

## License

This project is licensed under the **GNU General Public License v3.0** — see [`LICENSE`](LICENSE) for the full text. In short: you can use, modify, and redistribute, but any derivative work must also be GPL-3.0 and the authors provide no warranty.
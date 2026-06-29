# mirna-kras — miRNA binding profile classification

Unsupervised classification of 2,516 mature miRNAs by binding-affinity profile against 25 transcripts of the KRas oncogenic signaling network. The pipeline produces **functional signatures** via non-negative matrix factorization (NMF) and graph-based clustering.

## Main result

NMF at **k = 15** extracts 15 latent signatures; the top 5 (69.5% of miRNAs) separate the PI3K regulators (PIK3CA, PIK3CB) from the direct RAS effectors (KRAS, NRAS, HRAS). Only **3.5% of miRNAs are "mixed"** (membership entropy > 0.8).

Hard clustering methods (Leiden, Spectral, Girvan–Newman) all converge to a coarse bipartition on the same residual matrix and fail to isolate sub-structures — confirming the signal is better captured with soft membership.

## Navigating

<div class="grid cards" markdown>

- :material-rocket-launch: **[Getting started](getting-started.md)**

    Installation, dependencies, reproduction commands

- :material-sitemap: **[Architecture](architecture.md)**

    Data flow, DuckDB schema, conventions

- :material-database: **[Data model](data-model.md)**

    DuckDB schema for the 4 tables (transcripts, mirnas, results, leiden_clusters)

- :material-book-open-variant: **[API reference](reference/index.md)**

    Auto-generated docs from Google-style docstrings on the Python scripts

- :material-puzzle: **[Extension points](extension.md)**

    Add a transcript, a clustering method, try sparse-NMF

</div>

## Repo layout

```
.
├── scripts/              Python pipeline (15 scripts)
├── data/                 FASTA sources (miRNAs + transcripts)
├── figures/              Graph-clustering results
├── figures_nmf/          NMF results
├── doc/                  Human-facing documentation (markdown, Beamer)
├── kami-out/             Kami-styled developer doc (HTML + PDF)
├── mkdocs.yml            This documentation site config
└── LICENSE               GPL-3.0
```

## Building the docs

```bash
mkdocs serve              # dev server, hot reload at http://127.0.0.1:8000
mkdocs build              # static site to ./site
```

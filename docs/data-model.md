# Data model — DuckDB

Four tables in `rimap_results.duckdb`. Schema is fixed at creation time by `rimap_pipeline.py`.

| Table | Rows | Role |
|---|---|---|
| `transcripts` | 25 | Target transcript metadata |
| `mirnas` | 2,588 | Index of mature miRNAs |
| `results` | 1,208,803 | **Central table** — one binding site per row |
| `leiden_clusters` | 2,516 | Pre-computed Leiden clustering cache |

## `transcripts`

| Column | Type | Description |
|---|---|---|
| `transcript_id` | VARCHAR PK | ENST id (e.g. `ENST00000256078`) |
| `gene_name` | VARCHAR | Gene symbol (KRAS, EGFR, …) |
| `category` | VARCHAR | `kras_pathway` / `tumor_suppressor` / `mapk_pathway` |
| `enst_id` | VARCHAR | Versioned Ensembl ID (e.g. `KRAS-201`) |

## `mirnas`

| Column | Type | Description |
|---|---|---|
| `mirna_id` | VARCHAR PK | MIMAT id (e.g. `MIMAT0000062`) |
| `name` | VARCHAR | Symbolic name (e.g. `let-7a-5p`) |

## `results` — central table

Composite PK: `(transcript_id, mirna_id, position_start)`.

| Column | Type | Description |
|---|---|---|
| `transcript_id` | VARCHAR PK | FK to `transcripts` |
| `mirna_id` | VARCHAR PK | FK to `mirnas` |
| `binding_dG` | FLOAT | Binding free energy (kcal/mol) — negative, lower = more stable |
| `Kd` | FLOAT | Dissociation constant |
| `position_start` | BIGINT PK | Site start position on the transcript |
| `position_end` | BIGINT | Site end position |
| `region` | VARCHAR | `5UTR` / `CDS` / `3UTR` |
| `seed_type` | VARCHAR | Seed type (6mer, 7mer-m8, 7mer-a1, 8mer) |
| `phastcons` | FLOAT | Conservation score |
| `bridge` | INT | 1 if compensatory 3' bridging interaction |
| `validated` | BOOLEAN | Site experimentally validated |
| `interaction_string` | VARCHAR | Formatted site detail |

## `leiden_clusters` — cache

Pre-computed by `transcript_clustering.py`. Regenerated on every run. Caches the clustering output so a viz-only change doesn't recompute it.
# Docker — running the full pipeline in a container

The image builds everything from source at build time, so a single `docker build` is enough — no host setup required beyond Docker itself.

## What's bundled

| Dependency | Source in repo | Built inside the image by |
|---|---|---|
| Python 3.14 + Python deps | `requirements.txt` | `pip install -r` |
| Java 21 (OpenJDK) | upstream Debian packages | `apt-get` |
| MC-Flashfold (C) | `docker/build-context/MC-Flashfold-v37.0/` (33 MB) | `gcc src/flashfold.c` |
| RINexus Java classes | `docker/build-context/rimap_java_code/` (1.2 MB, pom.xml + src/) | `mvn package` |

The image is fully reproducible offline once the build context is in place. No internet is required at build time beyond the initial apt/pip downloads (which are cached in Docker layers).

## Quick start

```bash
cd mirna-kras

# Build the image (~5-10 min the first time)
docker build -t mirna-kras .

# Run the docs site at http://localhost:8000
docker run --rm -p 8000:8000 mirna-kras

# Drop into a bash shell inside the container, with a volume for outputs
docker run --rm -it -v $(pwd)/work:/app/work mirna-kras bash

# Build the static site
docker run --rm -v $(pwd)/site:/app/site mirna-kras mkdocs build --strict

# Run the full pipeline (Stages 0 → 1 → 2)
docker run --rm -v $(pwd)/work:/app/work mirna-kras bash -c '
  set -euo pipefail
  cd /app/work
  python /app/scripts/fetch_gencode.py
  python /app/scripts/rimap_pipeline.py            # or rimap_pipeline_local.py
  for s in nmf_analysis transcript_clustering chord_diagram heatmap_by_cluster \
           bipartite_by_cluster bipartite_graph bipartite_top_heatmap mirna_interactive_plot; do
    python /app/scripts/$s.py || true
  done
'
```

## What works inside the image

| Stage | Script | Status |
|---|---|---|
| 0 | `scripts/fetch_gencode.py` | needs internet to hit Ensembl REST |
| 0 | `scripts/extract_transcripts.py` | needs the GENCODE archive as a local file (drop in `/app/work`) |
| 1 | `scripts/rimap_pipeline.py` (API) | needs internet to hit RIMap-RISC |
| 1 | `scripts/rimap_pipeline_local.py` | needs `gencode.v46.transcripts.fa.gz` mounted; uses MC-Flashfold + RINexus from the image |
| 2 | `transcript_clustering.py`, `nmf_analysis.py`, … | reads the DuckDB at `/app/work/rimap_results.duckdb` if mounted |
| — | `mkdocs serve` / `mkdocs build` | no external deps needed |

## Mount points

| Host path (suggested) | Container path | Purpose |
|---|---|---|
| `$(pwd)/work` | `/app/work` | Drop generated `.duckdb`, `figures/`, `figures_nmf/` here between container runs |
| `(none)` | `/opt/MC-Flashfold/tables/` | Energy parameter tables — already inside the image |
| `(none)` | `/app/archive/rimap_java_code/target/classes/` | RINexus Java classes — already inside the image |

## Image size

Approximate, after a clean build:

| Layer | Size |
|---|---|
| `python:3.14-slim-bookworm` | ~150 MB |
| `openjdk-21-jdk-headless` | ~250 MB |
| `maven` (Stage 2 build only — discarded) | 0 (final image has no Maven) |
| `gcc` (Stage 1 build only — discarded) | 0 (final image has no gcc) |
| Python deps from `requirements.txt` | ~500 MB |
| Repo + MC-Flashfold tables | ~70 MB |

**Total: ~970 MB** — fits comfortably in a 2 GB CI cache.

## Ports and environment

| | |
|---|---|
| Port `8000` | mkdocs dev server (only when you run `mkdocs serve`) |
| `MCFLASHFOLD_PATH` | `/opt/MC-Flashfold` (points inside the image) |
| `JAVA_CLASSPATH` | `/app/archive/rimap_java_code/target/classes` |
| `PATH` | includes `/usr/local/bin` so `mcff` and `db2cm` are on it |

## Updating the bundled sources

| Source | When to refresh | How |
|---|---|---|
| MC-Flashfold | when the upstream releases a new version | extract the new ZIP into `docker/build-context/MC-Flashfold-vX.Y.Z/`, update the `COPY` line in the Dockerfile, bump the path |
| RINexus | when you push a fix upstream | pull the new source into `docker/build-context/rimap_java_code/`, commit, rebuild |

## Troubleshooting

- **`mcff not found` inside the container** — the MC-Flashfold C build failed silently. Check `docker build -t mirna-kras .` output for gcc errors.
- **`ClassNotFoundException` for `ca.iric.major.rinexus.rimaprisc.MiRScanS`** — the Maven build failed. Run `docker build --progress=plain` to see Maven logs.
- **`mkdocs serve` exits immediately** — you forgot to publish the port: `docker run --rm -p 8000:8000 mirna-kras`.
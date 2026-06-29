# Docker — running the full pipeline in a container

The image builds everything from source at build time. A fresh third-party user only needs three commands:

```bash
git clone <repo-url> mirna-kras
cd mirna-kras
./bin/mirna-kras build      # ~5–10 min the first time
./bin/mirna-kras test-local # 3-pair smoke test (~15 s)
```

The `bin/mirna-kras` CLI wraps every common Docker interaction (build, serve, shell, run, test, pipeline, clean) so the user does not have to remember the long `docker run -v ... -e ...` invocations.

## The `bin/mirna-kras` CLI

```
mirna-kras [-h] <subcommand> ...

Subcommands:
  build      Build the Docker image
  serve      Start the mkdocs dev server at http://localhost:8000
  shell      Open an interactive shell inside the container
  test-local Run the local pipeline smoke test (3 pairs, new DB)
  test-api   Run the API smoke test (test_curl.py)
  run        Run any script from scripts/ inside the container
  pipeline   Run the full local pipeline (WARNING: takes hours)
  python     Open a Python REPL with the pipeline scripts on sys.path
  clean      Remove the Docker image and dangling build cache
```

Every subcommand accepts `--work DIR` (default `./work`) for the work directory mounted at `/tmp/work` in the container. Files written there are owned by the host user (the CLI runs Docker as `-u $(id -u):$(id -g)`).

### Typical third-party workflow

```bash
# 1. Build the image (one-time, 5-10 min)
./bin/mirna-kras build

# 2. Smoke-test the local pipeline (3 KRAS-201 pairs, new DB, no API)
./bin/mirna-kras test-local
# → 48 sites in a new DB at ./work/rimap_results.duckdb
# → Same counts as a manual run, no production DB touched

# 3. Read the docs locally
./bin/mirna-kras serve    # → http://localhost:8000

# 4. Drop into a shell to run analysis scripts
./bin/mirna-kras shell    # → bash inside the image, with scripts/ on sys.path
# In the shell:
#   $ python /app/repo/scripts/nmf_analysis.py --output-dir /tmp/work/nmf

# 5. Clean up when done
./bin/mirna-kras clean
```

### Running a one-off script

```bash
./bin/mirna-kras run nmf_analysis.py --output-dir /tmp/work/nmf --n-seeds 20
```

This is equivalent to:
```bash
docker run --rm -u $(id -u):$(id -g) \
  -v $(pwd):/app/repo -v $(pwd)/work:/tmp/work \
  -e DB_PATH=/tmp/work/rimap_results.duckdb \
  -e GENCODE_FASTA=/tmp/work/gencode.v46.transcripts.fa.gz \
  -e JAVA_CLASSPATH=/app/archive/rimap_java_code/target/classes:/opt/rinexus-deps/* \
  mirna-kras python3 /app/repo/scripts/nmf_analysis.py --output-dir /tmp/work/nmf --n-seeds 20
```

## What's bundled

| Dependency | Source in repo | Built inside the image by |
|---|---|---|
| Python 3.14 + Python deps | `requirements.txt` | `pip install -r` |
| Java 21 (Temurin) | `eclipse-temurin:21-jdk` (upstream Adoptium) | base image |
| MC-Flashfold (C) | `docker/build-context/MC-Flashfold-v37.0/` (33 MB) | `gcc src/flashfold.c` |
| RINexus Java classes | `docker/build-context/rimap_java_code/` (1.2 MB) | `mvn package` + `dependency:copy-dependencies` |
| miRBase miRNA file | `data/mirnas.fa` (205 KB) | `COPY` |

The image is fully reproducible offline. The only network access is during `apt-get` / `pip install` for system deps; after the first build, subsequent runs use cached layers.

## What's required at runtime

The local pipeline (`rimap_pipeline_local.py`) needs the GENCODE v46 transcript archive. It's not bundled in the image (83 MB is too big). Provide it via the `--gencode` flag:

```bash
./bin/mirna-kras pipeline --gencode /path/to/gencode.v46.transcripts.fa.gz --yes
```

Download it from <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz>.

If you only want to use the API variant (`rimap_pipeline.py`, which queries RIMap-RISC over HTTPS), no GENCODE archive is needed.

## What works inside the image

| Stage | Script | Status |
|---|---|---|
| 0 | `scripts/fetch_gencode.py` | needs internet to hit Ensembl REST |
| 0 | `scripts/extract_transcripts.py` | needs the GENCODE archive as a local file (drop in `/tmp/work`) |
| 1 | `scripts/rimap_pipeline.py` (API) | needs internet to hit RIMap-RISC |
| 1 | `scripts/rimap_pipeline_local.py` | needs `gencode.v46.transcripts.fa.gz` mounted; uses MC-Flashfold + RINexus from the image |
| 2 | `transcript_clustering.py`, `nmf_analysis.py`, … | reads the DuckDB at `/tmp/work/rimap_results.duckdb` if mounted |
| — | `mkdocs serve` / `mkdocs build` | no external deps needed |

## Mount points

| Host path (suggested) | Container path | Purpose |
|---|---|---|
| `$(pwd)/work` (or `--work DIR`) | `/tmp/work` | Drop generated `.duckdb`, `figures/`, `figures_nmf/` here between container runs |
| `(none)` | `/opt/MC-Flashfold/tables/` | Energy parameter tables — already inside the image |
| `(none)` | `/app/archive/rimap_java_code/target/classes/` | RINexus Java classes + Jackson + friends — already inside the image |
| `--gencode PATH` | `/tmp/work/gencode.v46.transcripts.fa.gz` | GENCODE v46 transcript archive (mount the host file) |

## Image size

| Layer | Size |
|---|---|
| `python:3.14-slim-bookworm` | ~150 MB |
| `openjdk-21-jre-headless` (from Temurin) | ~50 MB |
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
| `JAVA_CLASSPATH` | `/app/archive/rimap_java_code/target/classes:/opt/rinexus-deps/*` |
| `PATH` | includes `/usr/local/bin` so `mcff` and `db2cm` are on it |

## Updating the bundled sources

| Source | When to refresh | How |
|---|---|---|
| MC-Flashfold | when the upstream releases a new version | extract the new ZIP into `docker/build-context/MC-Flashfold-vX.Y.Z/`, update the `COPY` line in the Dockerfile, bump the path |
| RINexus | when you push a fix upstream | pull the new source into `docker/build-context/rimap_java_code/`, commit, rebuild |
| `data/mirnas.fa` | when miRBase updates | replace the file in the repo, commit, rebuild |

## Troubleshooting

- **`mcff not found` inside the container** — the MC-Flashfold C build failed silently. Check `docker build -t mirna-kras .` output for gcc errors.
- **`ClassNotFoundException: com.fasterxml.jackson.databind.ObjectMapper`** — the Maven build didn't extract Jackson deps. Re-run `./bin/mirna-kras build` and look for "RINexus deps: N jars" in the build log.
- **`mkdocs serve` exits immediately** — you forgot to publish the port: `./bin/mirna-kras serve` already does `-p 8000:8000`.
- **`Permission denied` writing to `./work/`** — a previous run created the dir as root. `sudo chown -R $(id -u):$(id -g) ./work` to fix.
- **GENCODE archive not found** — the CLI prints a clear error with the download URL. Pass `--gencode PATH` to override the default.
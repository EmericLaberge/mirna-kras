# `data/` — source sequences

Two types of FASTA sequences, consumed by `scripts/fetch_gencode.py` (transcripts) and by `scripts/rimap_pipeline*.py` (transcripts + miRNAs).

## `mirnas.fa`

- **2,588 mature human miRNAs** (miRBase format, `MIMAT*` IDs).
- One sequence per miRNA (~22 nt).
- Used as-is by RIMap-RISC to predict binding sites on transcripts.

## `fasta/`

- **25 transcripts** of the KRas signaling network, one file per transcript (Ensembl ENST ID).
- Sequences from GENCODE v46 (via `fetch_gencode.py` or `extract_transcripts.py`).
- Three functional sub-groups:

| Pathway | # | Transcripts |
|---|---|---|
| KRas (direct effectors + MAPK) | 11 | KRAS, NRAS, HRAS, PIK3CA, PIK3CB, MTOR, BRAF, RAF1, MAP2K1, MAP2K2, MAPK1 |
| Tumor suppressors | 3 | TP53, PTEN, RB1 |
| MAPK / receptor tyrosine kinases | 10 | SOS1, GRB2, SHC1, EGFR, ERBB2, ERBB3, MET, IGF1R, FGFR1, FGFR2 |

> **Note:** the LaTeX report referenced a different set of 24 transcripts (with NRAS, HRAS, BRAF, PIK3CB, SOS1, …). The list above is the one actually present in `data/fasta/` and used by the pipeline. If you regenerate the analysis, update any external documentation accordingly.

## Provenance

- Transcripts: GENCODE v46 → Ensembl REST API (or local archive `gencode.v46.transcripts.fa.gz` for `extract_transcripts.py`).
- miRNAs: miRBase release 22, file `mirnas.fa` provided in the repo.

See [`scripts/README.md`](../scripts/README.md) for the download pipeline.
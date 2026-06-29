# Extension points

The code is deliberately modular: adding an analysis or a method means adding a script to `scripts/`, without touching the rest.

## You want… → Touch…

| You want… | Action |
|---|---|
| Add a transcript to the panel | `gene_ids.py` + `extract_transcripts.py` / `fetch_gencode.py` + **6× `TRANSCRIPT_GENE_MAP`** + re-run stage 1 |
| Add a clustering method | Create a script in `scripts/` that reads DuckDB, writes to `figures/` |
| Change the residualization threshold | `RESIDUAL_THRESHOLD` in `nmf_analysis.py` (search `V_res`) — default: 0 |
| Try another decomposition (ICA, sparse-NMF) | Reuse the `W` / `H` matrices, compare with cophenetic correlation |
| Experimentally validate a signature | Read `W_matrix.npy` (k=15) + `results` to identify the miRNAs, follow standard protocol |
| Switch underlying DB (SQLite, Postgres) | Replace `duckdb` in the imports — the SQL schema is portable |

## Adding a transcript — detailed checklist

1. **Constant**: add the ENST to `KRAS_PATHWAY` (or another category) at the top of `scripts/extract_transcripts.py` and `scripts/fetch_gencode.py`.
2. **Constant**: add `(symbol, ENSG_id)` to `scripts/gene_ids.py`.
3. **Mapping**: add `"ENSTxxxxxxxx": "GENE_NAME"` to the 6 `TRANSCRIPT_GENE_MAP` dictionaries (or equivalent `gene_pathway`):
   - `scripts/chord_diagram.py`
   - `scripts/heatmap_by_cluster.py`
   - `scripts/bipartite_top_heatmap.py`
   - `scripts/nmf_analysis.py`
   - `scripts/bipartite_by_cluster.py`
   - `scripts/bipartite_graph.py`
4. **Re-run stage 0 + 1**: `extract_transcripts.py` then `rimap_pipeline.py`.
5. **Verify**: `SELECT COUNT(*) FROM results WHERE transcript_id = 'ENSTxxxxxxxx'` — must be > 0.

## Adding a clustering method

Create a script that:

1. Opens `rimap_results.duckdb` with `read_only=True`.
2. Reads the `results` table, aggregates with `mean(|binding_dG|)` per pair.
3. Applies the residualization `max(0, V - V.mean(axis=0))`.
4. Runs your algorithm.
5. Writes its figures to `figures/`.

Minimal pattern (pseudo-code):

```python
import duckdb
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

con = duckdb.connect("rimap_results.duckdb", read_only=True)
df = con.execute("SELECT * FROM results").fetchdf()
con.close()

# aggregate mean(|binding_dG|)
V = (df.assign(w=df["binding_dG"].abs())
       .groupby(["mirna_id", "transcript_id"])["w"].mean()
       .unstack(fill_value=0).values)

# residualize
V_res = np.clip(V - V.mean(axis=0), 0, None)

# your algorithm...
# plt.savefig("figures/my_new_method.png")
```

## Known limitations

- **No unit tests.** Validation is visual (figures) + sanity checks on the metrics.
- **No shared cache between methods.** Each analysis re-aggregates `mean(|binding_dG|)` from `results`. Cost is negligible (~5 s per script) but not optimal.
- **Transcript discrepancy.** Some external documentation references a slightly different set of 25 transcripts. The pipeline and figures use the contents of `data/fasta/`. If you regenerate the analysis, align any external docs with what's actually on disk.
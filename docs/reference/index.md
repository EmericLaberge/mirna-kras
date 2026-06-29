# API Reference

Auto-generated documentation from the Google-style docstrings on the Python scripts. Each page covers a thematic subset of the project's scripts.

## Pages

| Page | Scripts | Role |
|---|---|---|
| [Data](data.md) | `gene_ids`, `fetch_gencode`, `extract_transcripts` | Stage 0 — FASTA sources |
| [Rimap computation](rimap.md) | `rimap_pipeline`, `rimap_pipeline_local` | Stage 1 — interaction computation |
| [Graph clustering](clustering.md) | `transcript_clustering` | Stage 2 — Leiden/Spectral/GN |
| [NMF](nmf.md) | `nmf_analysis` | Stage 2 — non-negative factorization |
| [Visualizations](viz.md) | `chord_diagram`, `heatmap_by_cluster`, `bipartite_*`, `mirna_interactive_plot` | Stage 2 — figures |
| [Utilities](utilities.md) | `test_curl` | Off-pipeline |

## Docstring conventions

Every module, class, and function follows the **Google style** format:

```python
def example(arg1: str, arg2: int = 5) -> bool:
    """Short summary.

    Optional longer description.

    Args:
        arg1: Description of the first argument.
        arg2: Description of the second. Defaults to 5.

    Returns:
        Description of the return value.

    Raises:
        ValueError: When something is wrong.
    """
```

This format is parsed by mkdocstrings + griffe and rendered as HTML on these pages.
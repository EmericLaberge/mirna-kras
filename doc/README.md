# `doc/` — project documentation

| File / folder | Contents | When to read |
|---|---|---|
| [`nmf_analysis.md`](nmf_analysis.md) | Technical NMF doc: algorithm, rank selection (cophenetic correlation), signatures 1 & 2, comparison with hard clustering. | To understand *why* k = 15 (or k = 2) and how normalized entropy measures "mixed" miRNAs. |
| [`presentation/`](presentation/) | Beamer presentation (12 sections, PDF + LaTeX sources). | For preparing the oral defense. |

## How it fits together

The mkdocs site (`docs/`) is the canonical, browsable technical doc. This `doc/` folder carries:

1. **Method-level detail** for NMF — `nmf_analysis.md` is the deep dive.
2. **Defense material** — `presentation/` for the slide deck.

If you don't know where to start: read [`README.md`](../README.md) for the project overview, then jump to the [mkdocs site](../docs/index.md) for the technical walk-through.
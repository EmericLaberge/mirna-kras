# NMF — Stage 2

Non-negative matrix factorization with rank selection by cophenetic correlation (Brunet et al., 2004).

## `scripts.nmf_analysis`

Full pipeline: load → residualize → rank selection (k=2..15, 20 seeds) → final fit → figures + W/H matrices.

::: scripts.nmf_analysis
# NMF — Non-Negative Matrix Factorization for miRNA Signatures

## Goal

Identify **latent functional signatures** in miRNA–transcript binding profiles via NMF ($V \approx WH$). Unlike hard clustering (Leiden, Spectral) which assigns each miRNA to a single cluster, NMF produces a **soft membership**: each miRNA contributes to several signatures with continuous weights.

## Algorithm

### Input data

- **Matrix $V$**: 2,516 miRNAs × 25 transcripts, values = mean(|binding_dG|).
- Source: `rimap_results.duckdb`, RIMap-RISC interactions (excluding seedless ones).
- The |dG| values lie in [10, 26] — all miRNAs target almost all transcripts (98% density).

### Preprocessing — residual matrix

The raw data has a problem: binding profiles are nearly identical (coefficient of variation of row sums = 0.06). A naive NMF places 100% of miRNAs in a single signature.

**Solution**: compute the residual matrix $V_{res}$ where each element is the positive deviation from the column mean:

$$V_{res}[i,j] = \max(0,\ V[i,j] - \bar{V}[:,j])$$

This extracts the **differential signal** — miRNAs that bind more strongly than average for a given transcript. MiRNAs below the mean get zeroed out (weak binding = no signal). The residual matrix is 58.7% non-negative, which is compatible with NMF.

### NMF decomposition

$$V_{res} \approx W \cdot H$$

- **$W$** (2,516 × $k$): membership matrix — each row gives a miRNA's weight in each signature.
- **$H$** ($k$ × 25): signature matrix — each row is the transcript profile defining one signature.
- Solver: Multiplicative Updates (MU), Frobenius divergence.
- Initialization: NNDSVD for the final fit (deterministic, structured).

### Rank selection ($k$)

For each $k \in [2, 25]$:

1. Run 20 NMFs with different random seeds.
2. Build a **consensus matrix** (co-clustering frequency for each miRNA pair).
3. Compute the **cophenetic correlation** of this matrix via hierarchical clustering.
4. The optimal $k$ is the one with the highest cophenetic correlation (maximal stability).

> **Note**: $k = 25$ is the theoretical maximum because the matrix has only 25 features (transcripts). The NNDSVD initialization requires $k \leq \min(n_{\text{samples}}, n_{\text{features}})$.

### Mixed-miRNA detection

The normalized entropy of each row of $W$ measures whether a miRNA has a dominant signature or is shared:

$$H_{norm}(i) = \frac{-\sum_j w_{ij} \log w_{ij}}{\log k}$$

If $H_{norm} > 0.8$, the miRNA is considered "mixed" — it does not have a clear preferential affinity.

## Results

### Optimal rank: $k = 2$

| $k$ | Cophenetic | Error | $k$ | Cophenetic | Error |
|-----|-----------|-------|-----|-----------|-------|
| **2** | **0.9795** | 152.2 | 14 | 0.9329 | 93.8 |
| 3 | 0.9274 | 143.6 | 15 | 0.9416 | 91.1 |
| 4 | 0.9023 | 134.9 | 16 | 0.9297 | 88.6 |
| 5 | 0.9415 | 127.1 | 17 | 0.9325 | 85.6 |
| 6 | 0.9394 | 121.8 | 18 | 0.9183 | 82.8 |
| 7 | 0.9238 | 117.1 | 19 | 0.9268 | 80.5 |
| 8 | 0.9304 | 112.6 | 20 | 0.9309 | 78.0 |
| 9 | 0.9187 | 108.7 | 21 | 0.9295 | 75.8 |
| 10 | 0.9226 | 105.3 | 22 | 0.9182 | 73.6 |
| 11 | 0.9173 | 102.4 | 23 | 0.9148 | 71.7 |
| 12 | 0.9190 | 99.5 | 24 | 0.9006 | 70.2 |
| 13 | 0.9277 | 96.5 | 25 | 0.8905 | 68.9 |

The highest cophenetic correlation is at $k = 2$ (0.9795). There is a secondary plateau around $k = 5$–$6$ (0.939–0.942) and a local peak at $k = 15$ (0.942). The reconstruction error decreases monotonically, which is expected — more components = better fit. The cophenetic decreases globally after $k = 15$, suggesting that structures beyond are unstable.

### Signature 1 — RAS-MAPK pathway (1,524 miRNAs, 60.6%)

| Transcript | NMF weight | mean |dG| |
|-----------|-----------|-----------|
| KRAS | 6.21 | 18.96 |
| PIK3CB | 5.11 | 18.92 |
| NRAS | 3.98 | 18.97 |
| HRAS | 3.84 | 18.83 |
| MTOR | 3.29 | 18.83 |

miRNAs with differential affinity for **KRAS, NRAS, PIK3CB** — the main effectors of the RAS-MAPK pathway. These miRNAs target the oncogenic RAS signal more strongly than average.

### Signature 2 — PI3K-AKT pathway (992 miRNAs, 39.4%)

| Transcript | NMF weight | mean |dG| |
|-----------|-----------|-----------|
| PIK3CA | 8.26 | 18.83 |
| HRAS | 5.55 | 18.39 |

miRNAs with differential affinity for **PIK3CA** — the central catalyst of the PI3K-AKT pathway. These miRNAs show a bias toward regulating cell survival via PI3K.

### Mixed miRNAs

- **1,272 miRNAs (50.6%)** are "mixed" (entropy > 0.8).
- Mean entropy: 0.567.
- These miRNAs contribute to both signatures without clear dominance — they reflect the functional continuum of the KRas network.

## Biological interpretation

The bipartition NMF separates miRNAs into two functional groups consistent with KRas biology:

1. **RAS-MAPK group**: miRNAs whose differential affinity targets the direct effectors of RAS signaling (KRAS, NRAS) and the PI3K beta kinase (PIK3CB). These miRNAs may be involved in regulating cellular proliferation.

2. **PI3K-AKT group**: miRNAs preferentially targeting PIK3CA (the alpha catalytic subunit of PI3K). These miRNAs may modulate the AKT cell-survival pathway.

The fact that 50% of miRNAs are "mixed" is expected — the KRas network forms a functional continuum where most miRNAs have similar binding profiles. NMF captures the subtle differences at the top of the affinity distribution.

## Comparison with classic clustering

| Method | Clusters | Best metric | Membership |
|--------|----------|-------------|------------|
| Leiden (res=0.5) | 7 | Sil=0.496, Mod=0.450 | Hard |
| Spectral (k=3) | 3 | Sil=0.413, Mod=0.399 | Hard |
| NMF (k=2) | 2 | Coph=0.980 | Soft |
| Girvan–Newman | 1 | — | — |

NMF offers a complementary perspective: instead of forcing strict boundaries, it quantifies the degree to which each miRNA belongs to each signature. This is more nuanced biologically, especially for a network this dense.

## Files

```
scripts/nmf_analysis.py      # NMF pipeline (542 lines)
figures_nmf/
├── rank_selection.png        # Cophenetic + error curves vs k (k=2..25)
├── k2/                       # Results k=2
│   ├── heatmap.png           # Heatmap H (signatures × transcripts)
│   ├── scatter.png           # Scatter miRNAs (W space)
│   ├── signature_1_bars.png  # Top transcripts signature 1
│   └── signature_2_bars.png  # Top transcripts signature 2
├── k3/                       # Results k=3 (3 signatures)
├── k4/                       # Results k=4 (4 signatures)
└── k5/                       # Results k=5 (5 signatures)
```

## Usage

```bash
.venv/bin/python3 scripts/nmf_analysis.py --output-dir figures_nmf --n-seeds 20
```

Options: `--k-min`, `--k-max`, `--max-iter`, `--mixed-threshold`
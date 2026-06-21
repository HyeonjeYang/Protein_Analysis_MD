# Heatmaps

The unified heatmap helpers in `protein_analysis_md.visualization.heatmaps`
accept NumPy arrays or pandas matrices and add colorbars, units, and optional
sequence annotations.

Supported heatmap data include raw contact maps, observed/expected-like contact
maps, delta maps, residue-class contact matrices, PTM-site profiles,
fragment-fragment matrices, inter-chain matrices, density projections,
PCA/contact loading maps, sequence-feature correlation matrices, and
condition-by-metric summaries.

Observed/expected-like maps use local expected-by-separation logic:

```text
expected(s) = mean C_ij for |i-j| = s
OE_ij = log((C_ij + eps) / (expected(|i-j|) + eps))
```

Matrix smoothing is disabled by default. If a smoothed display matrix is
generated, raw matrix data must still be saved separately and the figure should
be labeled `display-smoothed`.

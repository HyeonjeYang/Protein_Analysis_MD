# PCA And Contact-Environment Decomposition

`protein_analysis_md` includes optional exploratory decomposition analyses for
protein/IDR trajectories. These methods are useful for summarizing ensemble
shifts, but they are not validated as standalone biological classifiers.

## Coordinate PCA / Essential Dynamics

Coordinate PCA flattens each frame's bead coordinates and decomposes the main
directions of structural motion. The default preprocessing removes each frame's
center of mass and does not align frames. This is intentional: IDRs may not have
a stable reference structure, so Kabsch alignment can introduce artifacts.

Outputs include PC scores, coordinate loadings, explained variance, and
representative frames.

## Contact-Map PCA

Contact PCA converts each frame's contact map into residue-pair features and
identifies dominant contact-pattern modes. This can help detect PTM- or
cleavage-induced shifts in contact ensembles.

## Distance-Map PCA

Distance PCA uses per-frame pairwise distance maps. A log transform is available
to reduce the influence of very large distances. Results should be interpreted
as geometric ensemble modes, not mechanistic pathways.

## Feature PCA

Feature PCA operates on per-frame scalar observables such as Rg, Ree, total
contacts, long-range contacts, and optional energy-like features. Missing
optional columns are ignored; PCA is computed from available numeric features.

## Contact-Environment Eigendecomposition

The contact-environment eigendecomposition is an analogy to Hi-C
observed/expected compartment decomposition:

1. Average the contact map.
2. Estimate P(s), where `s = |i - j|`.
3. Compute an observed/expected-like matrix by log-ratio or difference.
4. Mask short sequence separations.
5. Correlate rows/columns by contact environment.
6. Eigen-decompose the correlation matrix.
7. Orient eigenvector signs using sequence features when available.

The output is labeled as contact-environment eigenvectors. Do not call these
chromosome compartments, and do not interpret EV1/EV2 as genomic A/B
compartments.

## NMF Contact Modules

NMF contact modules are optional and experimental. They decompose a non-negative
contact map into residue module weights and module contact maps. Use this as a
hypothesis-generation tool only.

## Comparing States

For WT-vs-PTM or intact-vs-cleaved comparisons, decomposition outputs can
summarize:

- EV1 correlation against WT or intact state.
- Delta EV1 along the sequence.
- PC centroid shifts for feature/contact PCA.
- EV1-correlation and PC-centroid trends across cleavage cut number or event
  time when metadata exist.

## Limitations

- These analyses are exploratory unless validated for a specific system.
- Results depend on sampling, contact cutoff, sequence-separation mask, and
  preprocessing choices.
- IDR coordinate PCA can be sensitive to alignment and center-of-mass removal.
- Contact-environment eigenvectors are an analogy to Hi-C decomposition, not a
  direct biological equivalent.

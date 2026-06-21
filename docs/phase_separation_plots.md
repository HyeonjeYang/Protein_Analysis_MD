# Phase-Separation Plots

Phase visualizations are report aids for slab, droplet-like, and multi-chain
simulations. They do not by themselves prove phase separation.

Implemented helper panels include:

- protein/component density profile vs z
- dense/dilute concentration time series when available
- dense/dilute distributions and component partition coefficients
- 2D projected density heatmaps
- cluster-size distributions
- inter-chain and homotypic/heterotypic contact heatmaps

The report code emits a reliability warning when a trajectory has too few
frames or chains for robust dense/dilute inference. Concentration summaries
should be based on raw or explicitly defined binned profiles, not arbitrary
display smoothing.

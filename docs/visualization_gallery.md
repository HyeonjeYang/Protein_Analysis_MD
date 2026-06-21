# Visualization Gallery

`protein_analysis_md` report figures are built from raw analysis outputs and
optional visualization overlays. Every new visualization helper can save PNG,
PDF, and the plotting data through `save_visualization`.

Single-chain reports can include Rg/Ree distributions, Rg/Ree time series,
joint Rg/Ree density, P(s), R(s), local scaling exponent, contact maps,
observed/expected-like contact maps, sequence tracks, and energy plots when the
corresponding files exist.

PTM reports can include WT-vs-PTM distributions, raw delta contact maps,
PTM-site contact profiles, residue-class contact changes, and exploratory
decomposition panels.

Cleavage reports can include cleavage maps, fragment architecture, raw event
schedules, fragment length distributions, cut-number trends, and contact-map
comparisons in original sequence coordinates.

Phase and multi-protein reports can include density profiles, cluster-size
distributions, inter-chain contact heatmaps, and partitioning summaries. These
plots are labeled exploratory when the system is short or small.

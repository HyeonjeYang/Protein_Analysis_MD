# Presets

Simulation presets include `smoke_single_chain`, `short_single_chain`,
`production_single_chain`, `cleavage_smoke`, `cleavage_production`,
`phase_smoke`, and `phase_slab_production`.

Analysis presets include `minimal`, `standard_idr`, `standard_cleavage`,
`standard_phase`, and `full`.

Report presets include `minimal`, `standard`, and `publication_draft`.
Preset values are defaults; user overrides win and are recorded in lock metadata.

`minimal` leaves smoothing disabled. `standard_idr` enables log-space smoothing
for P(s) and R(s) trend columns, while `standard_cleavage` also enables energy
time-series smoothing. Report presets control whether smoothed trend lines are
shown alongside raw points.

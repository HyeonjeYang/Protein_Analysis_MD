# Free-Energy Landscapes

`protein_analysis_md.analysis.free_energy` computes exploratory 2D free-energy
surfaces from raw binned counts:

```text
F(x, y) = -kBT ln P(x, y) + constant
```

Supported variable pairs include Rg/Ree, Rg/total contacts, PC1/PC2,
fragment-distance/fragment-contacts, and other per-frame scalar pairs when
available.

Outputs preserve raw counts:

- `*_counts.npy`
- `*_free_energy.npy`
- `*_x_edges.npy`
- `*_y_edges.npy`
- `*_grid.csv`
- `*_metadata.json`

Display smoothing can be used for visualization, but raw counts remain the
primary data. These surfaces are sampling-dependent and must not be used to
infer kinetic barriers from short trajectories.

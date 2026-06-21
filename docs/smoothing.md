# Smoothing And Coarse-Grained Curves

`protein_analysis_md` can append smoothed curves for noisy observables while
preserving the raw measurements. Smoothing is intended for visualization,
trend estimation, and report readability. It is not a replacement for raw
measurements, and scalar summaries use raw data unless `analysis.fit_to:
smoothed` is explicitly configured for the Flory-like exponent fit.

## Raw Data Policy

Analysis writers never overwrite raw columns. For example:

- `ps.parquet` keeps `p_contact` and may add `p_contact_smooth`.
- `scaling.parquet` keeps `mean_distance_nm` and may add
  `mean_distance_nm_smooth`.
- `timeseries_rg.parquet` keeps `rg` and may add `rg_nm_smooth`.
- `timeseries_ree.parquet` keeps `ree` and may add `ree_nm_smooth`.
- `energy.parquet` keeps raw energy and temperature columns and may add
  `*_smooth` columns.

Each output also receives a `*.units.json` sidecar. Analysis `summary.json`
records the smoothing methods, parameters, raw columns, and smoothed columns.
Compiled projects record smoothing settings in both `project.lock.yaml` and
`config_resolved.json`.

## Curve Smoothing

P(s) and R(s) use log-space smoothing by default in standard report-oriented
presets. For each sequence separation `s_i`, the smoother uses points whose
`log10(s)` values fall within the configured window. This follows the same
general idea as Hi-C expected-cis smoothing practice, such as cooltools-style
P(s) smoothing, but the implementation here is local and lightweight using
NumPy, SciPy, and pandas. The package does not depend on cooltools.

Example:

```yaml
analysis:
  smoothing:
    ps:
      enabled: true
      method: logspace
      window_log10: 0.2
      min_points: 5
      robust: true
    rs:
      enabled: true
      method: logspace
      window_log10: 0.2
      min_points: 5
      robust: true
```

## Time Series Smoothing

Rg, Ree, and energy-like time series support rolling and Savitzky-Golay
smoothing. These smoothed traces are visualization aids by default.

```yaml
analysis:
  smoothing:
    rg:
      enabled: true
      method: rolling
      window: 25
    energy:
      enabled: true
      method: rolling
      window: 25
```

## Contact Map Smoothing

Contact map smoothing is disabled by default and is visualization-only unless a
future analysis mode explicitly opts into quantitative use. When enabled,
analysis writes both `contact_map.npy` and `contact_map_smoothed.npy`. WT-vs-PTM
delta contact maps use raw contact maps by default.

```yaml
analysis:
  smoothing:
    contact_map:
      enabled: true
      method: gaussian
      sigma: 1.0
      visualization_only: true
```

## Reports

Report presets can show raw points with smoothed trend lines when smoothed
columns are available:

```yaml
report:
  smoothing:
    use_smoothed_ps: true
    use_smoothed_rs: true
    show_raw_points: true
    show_smoothed_line: true
    show_smoothing_metadata: true
```

Figure titles and report text indicate when a smoothed trend is displayed.

## Strict Policy

Raw data are the source of truth. Smoothed values are optional display layers
unless a user explicitly configures otherwise. Default summary metrics use raw
values: mean Rg, mean Ree, mean energy, mean contact probability,
dense/dilute concentration tables, and delta contact maps are not computed from
display-smoothed data.

Allowed or recommended for visualization:

| Target | Default | Method | Notes |
| --- | --- | --- | --- |
| P(s) | On in standard presets | log-space | Stores `p_contact` and `p_contact_smooth`. |
| R(s) | On in standard presets | log-space | Stores `mean_distance_nm` and `mean_distance_nm_smooth`. |
| Energy | On in standard presets | rolling | Visual drift/equilibration check. |
| Rg/Ree time series | Off by default | rolling | Visual trend only. |
| Density profiles | On for phase presets | Gaussian/rolling | Raw binned profile remains available. |
| Free-energy display | Optional | binned/KDE-like display smoothing | Raw counts are stored separately. |

Allowed only with strong warnings:

| Target | Default | Rule |
| --- | --- | --- |
| Contact maps | Off | Write `contact_map_smoothed.npy`; never use for default deltas. |
| Observed/expected heatmaps | Off | Display smoothing only. |
| Delta contact maps | Off | Smooth only after raw delta calculation. |
| EV1/EV2 tracks | Off | EV correlations use raw eigenvectors. |

Do not smooth by default:

| Target | Reason |
| --- | --- |
| Binary contact time series | Would alter state definitions. |
| Contact lifetime calculations | Quantitative correlation must use raw states. |
| Cleavage event schedules and Poisson times | Discrete events. |
| Cut-site positions and fragment boundaries | Discrete sequence coordinates. |
| Cluster membership/state labels | Categorical labels. |
| Sequence features, pI/NCPR/FCR/aromaticity | Deterministic sequence values. |

Smoothing metadata records method, window, minimum points, robust/median/mean
choice, visualization-only status, and whether a value was used for a
quantitative fit. If smoothing is ever used for a quantitative metric, reports
must include a warning.

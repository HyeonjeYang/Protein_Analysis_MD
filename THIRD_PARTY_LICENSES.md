# Third-Party Licenses

`idr-ptm-md` is licensed as GPL-3.0-only. It is designed as a wrapper workflow
around external scientific software and libraries; it does not vendor CALVADOS
or dependency source code.

## Simulation Backend

| Component | Role | License | Source |
| --- | --- | --- | --- |
| CALVADOS | External simulation backend | GPL-3.0 | <https://github.com/KULL-Centre/CALVADOS> |

Users must install CALVADOS separately. `idr-ptm-md` prepares run directories,
execution scaffolds, and analysis workflows without modifying upstream CALVADOS.

## Python Dependencies

| Component | Role | License |
| --- | --- | --- |
| NumPy | Numerical arrays | BSD-3-Clause |
| pandas | Tabular data handling | BSD-3-Clause |
| SciPy | Scientific algorithms | BSD-3-Clause |
| Matplotlib | Plotting | PSF-based / Matplotlib license |
| Pydantic | Configuration schemas | MIT |
| Typer | Command-line interface | MIT |
| PyYAML | YAML parsing | MIT |
| MDTraj | Trajectory I/O and analysis helper | LGPL-2.1-or-later |
| PyArrow | Parquet output support | Apache-2.0 |
| pytest | Test runner | MIT |
| Ruff | Linting and formatting checks | MIT |

Dependency license details should be rechecked during release packaging,
especially if pinned versions change.

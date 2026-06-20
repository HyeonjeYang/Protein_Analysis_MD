from __future__ import annotations

import importlib


def test_package_imports() -> None:
    import idrptm

    assert idrptm.__version__
    assert idrptm.SUPPORTED_MVP_PTMS == ("pSer", "pThr")


def test_stage_1_modules_import() -> None:
    modules = [
        "idrptm.schema",
        "idrptm.sequence",
        "idrptm.ptm",
        "idrptm.residue_params",
        "idrptm.design",
        "idrptm.calvados_adapter",
        "idrptm.runner",
        "idrptm.hpc",
        "idrptm.analysis.io",
        "idrptm.analysis.rg",
        "idrptm.analysis.ree",
        "idrptm.analysis.contacts",
        "idrptm.analysis.ps",
        "idrptm.analysis.scaling",
        "idrptm.analysis.msd",
        "idrptm.analysis.lifetime",
        "idrptm.analysis.sequence_features",
        "idrptm.analysis.compare",
        "idrptm.plotting.plots",
        "idrptm.plotting.report",
    ]

    for module in modules:
        assert importlib.import_module(module)

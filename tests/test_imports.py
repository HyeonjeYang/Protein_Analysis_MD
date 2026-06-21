from __future__ import annotations

import importlib


def test_package_imports() -> None:
    import idrptm
    import protein_analysis_md

    assert idrptm.__version__
    assert idrptm.SUPPORTED_MVP_PTMS == ("pSer", "pThr")
    assert protein_analysis_md.__version__ == idrptm.__version__


def test_stage_1_modules_import() -> None:
    modules = [
        "idrptm.schema",
        "idrptm.sequence",
        "idrptm.ptm",
        "idrptm.cleavage",
        "idrptm.residue_params",
        "idrptm.design",
        "idrptm.calvados_adapter",
        "idrptm.runner",
        "idrptm.hpc",
        "idrptm.units",
            "idrptm.storage",
            "idrptm.uniprot",
            "idrptm.presets",
            "idrptm.configuration",
            "idrptm.recipe",
            "idrptm.project",
            "idrptm.registry",
            "idrptm.enzymes",
        "idrptm.environment",
        "idrptm.sequence_features",
        "protein_analysis_md",
            "protein_analysis_md.schema",
            "protein_analysis_md.presets",
            "protein_analysis_md.recipe",
            "protein_analysis_md.analysis.rg",
        "idrptm.analysis.io",
        "idrptm.analysis.cleavage",
            "idrptm.analysis.energy",
            "idrptm.analysis.equilibration",
            "idrptm.analysis.phase",
        "idrptm.analysis.rg",
        "idrptm.analysis.ree",
        "idrptm.analysis.contacts",
        "idrptm.analysis.ps",
        "idrptm.analysis.scaling",
        "idrptm.analysis.smoothing",
        "idrptm.analysis.msd",
        "idrptm.analysis.multichain",
        "idrptm.analysis.lifetime",
        "idrptm.analysis.sequence_features",
        "idrptm.analysis.compare",
            "protein_analysis_md.analysis.smoothing",
        "idrptm.plotting.plots",
        "idrptm.plotting.report",
    ]

    for module in modules:
        assert importlib.import_module(module)

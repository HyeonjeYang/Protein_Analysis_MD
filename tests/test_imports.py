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
        "idrptm.environment_check",
        "idrptm.repo_check",
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
        "protein_analysis_md.units",
        "protein_analysis_md.sequence",
        "protein_analysis_md.uniprot",
        "protein_analysis_md.ptm",
        "protein_analysis_md.cleavage",
        "protein_analysis_md.enzymes",
        "protein_analysis_md.design",
        "protein_analysis_md.residue_params",
        "protein_analysis_md.calvados_adapter",
        "protein_analysis_md.runner",
        "protein_analysis_md.storage",
        "protein_analysis_md.environment_check",
        "protein_analysis_md.repo_check",
        "protein_analysis_md.presets",
        "protein_analysis_md.recipe",
        "protein_analysis_md.analysis.rg",
        "idrptm.analysis.io",
        "idrptm.analysis.cleavage",
        "idrptm.analysis.energy",
        "idrptm.analysis.equilibration",
        "idrptm.analysis.free_energy",
        "idrptm.analysis.phase",
        "idrptm.analysis.rg",
        "idrptm.analysis.ree",
        "idrptm.analysis.contacts",
        "idrptm.analysis.decomposition",
        "idrptm.analysis.ps",
        "idrptm.analysis.scaling",
        "idrptm.analysis.smoothing",
        "idrptm.analysis.msd",
        "idrptm.analysis.multichain",
        "idrptm.analysis.lifetime",
        "idrptm.analysis.sequence_features",
        "idrptm.analysis.compare",
        "protein_analysis_md.analysis.smoothing",
        "protein_analysis_md.analysis.decomposition",
        "protein_analysis_md.analysis.free_energy",
        "idrptm.plotting.plots",
        "idrptm.plotting.report",
        "idrptm.visualization",
        "idrptm.visualization.single_chain",
        "idrptm.visualization.ptm",
        "idrptm.visualization.cleavage",
        "idrptm.visualization.phase",
        "idrptm.visualization.heatmaps",
        "idrptm.visualization.sequence_tracks",
        "idrptm.visualization.decomposition",
        "idrptm.visualization.smoothing_policy",
        "idrptm.visualization.free_energy",
        "protein_analysis_md.visualization",
        "protein_analysis_md.visualization.single_chain",
        "protein_analysis_md.visualization.ptm",
        "protein_analysis_md.visualization.cleavage",
        "protein_analysis_md.visualization.phase",
        "protein_analysis_md.visualization.heatmaps",
        "protein_analysis_md.visualization.sequence_tracks",
        "protein_analysis_md.visualization.decomposition",
        "protein_analysis_md.visualization.smoothing_policy",
        "protein_analysis_md.visualization.free_energy",
    ]

    for module in modules:
        assert importlib.import_module(module)

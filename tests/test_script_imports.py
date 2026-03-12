"""Import-level coverage tests for all CLI scripts.

These tests import each script module in-process (not via subprocess) so that
top-level code is measured by pytest-cov.  Each test also verifies that the
module exposes a callable ``main`` entry point.
"""

import importlib

import pytest

SCRIPT_MODULES = [
    "scripts.bam_stat",
    "scripts.bam2fq",
    "scripts.bam2wig",
    "scripts.clipping_profile",
    "scripts.deletion_profile",
    "scripts.divide_bam",
    "scripts.FPKM_count",
    "scripts.FPKM_UQ",
    "scripts.geneBody_coverage",
    "scripts.geneBody_coverage2",
    "scripts.infer_experiment",
    "scripts.inner_distance",
    "scripts.insertion_profile",
    "scripts.junction_annotation",
    "scripts.junction_saturation",
    "scripts.mismatch_profile",
    "scripts.normalize_bigwig",
    "scripts.overlay_bigwig",
    "scripts.read_distribution",
    "scripts.read_duplication",
    "scripts.read_GC",
    "scripts.read_hexamer",
    "scripts.read_NVC",
    "scripts.read_quality",
    "scripts.RNA_fragment_size",
    "scripts.RPKM_saturation",
    "scripts.sc_bamStat",
    "scripts.sc_editMatrix",
    "scripts.sc_seqLogo",
    "scripts.sc_seqQual",
    "scripts.split_bam",
    "scripts.split_paired_bam",
    "scripts.tin",
]


@pytest.mark.parametrize("module_name", SCRIPT_MODULES)
def test_script_importable(module_name: str) -> None:
    """Each script module should be importable and expose a callable main()."""
    mod = importlib.import_module(module_name)
    assert hasattr(mod, "main"), f"{module_name} missing main()"
    assert callable(mod.main)

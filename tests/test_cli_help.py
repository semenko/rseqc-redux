"""Smoke tests: verify all CLI scripts exit cleanly with --help."""

import subprocess
import sys

import pytest

SCRIPTS = [
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


@pytest.mark.parametrize("module", SCRIPTS, ids=[s.split(".")[-1] for s in SCRIPTS])
def test_script_help(module):
    """Each CLI script should exit 0 and print usage info with --help."""
    result = subprocess.run(
        [sys.executable, "-m", module, "--help"],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, f"{module} --help failed: {result.stderr}"
    output = (result.stdout + result.stderr).lower()
    assert "usage" in output or "option" in output, f"{module} --help produced no usage text"

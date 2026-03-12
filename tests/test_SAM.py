"""Tests for rseqc.SAM."""

import importlib
import io
import sys
from pathlib import Path

from rseqc import SAM

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def test_import():
    """Verify that rseqc.SAM can be imported."""
    mod = importlib.import_module("rseqc.SAM")
    assert mod is not None


def test_parsebam_exists():
    """ParseBAM should still be importable."""
    assert hasattr(SAM, "ParseBAM")


def test_parsesam_removed():
    """ParseSAM dead code class should no longer exist."""
    assert not hasattr(SAM, "ParseSAM")


def test_qcsam_removed():
    """QCSAM dead code class should no longer exist."""
    assert not hasattr(SAM, "QCSAM")


def test_parsebam_init(mini_bam):
    """ParseBAM can be instantiated with a BAM file."""
    obj = SAM.ParseBAM(str(mini_bam))
    assert obj is not None
    assert obj.bam_format is True


def test_parsebam_stat(mini_bam, capsys):
    """ParseBAM.stat() should report mapping statistics."""
    obj = SAM.ParseBAM(str(mini_bam))
    obj.stat(q_cut=30)
    captured = capsys.readouterr()
    assert "Total records:" in captured.out


def test_parsebam_configure_experiment(mini_bam):
    """ParseBAM.configure_experiment() with mini BED should not crash."""
    obj = SAM.ParseBAM(str(mini_bam))
    bed_file = str(FIXTURES_DIR / "mini.bed")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.configure_experiment(refbed=bed_file, sample_size=100)
    finally:
        sys.stderr = old_stderr

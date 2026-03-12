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


def test_parsebam_readsNVC(mini_bam, tmp_path):
    """readsNVC() should produce a .NVC.xls file with per-position base counts."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.readsNVC(outfile=outprefix, q_cut=30)
    finally:
        sys.stderr = old_stderr
    xls_file = Path(outprefix + ".NVC.xls")
    assert xls_file.exists()
    content = xls_file.read_text()
    assert "Position" in content
    # Should have A, C, G, T columns
    assert "\tA\t" in content


def test_parsebam_readGC(mini_bam, tmp_path):
    """readGC() should produce a .GC.xls file with GC content distribution."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.readGC(outfile=outprefix, q_cut=30)
    finally:
        sys.stderr = old_stderr
    xls_file = Path(outprefix + ".GC.xls")
    assert xls_file.exists()
    content = xls_file.read_text()
    assert "GC%" in content
    # Should have at least one data row
    lines = [ln for ln in content.strip().split("\n") if not ln.startswith("GC%")]
    assert len(lines) >= 1


def test_parsebam_readDupRate(mini_bam, tmp_path):
    """readDupRate() should produce seq and pos DupRate files."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.readDupRate(q_cut=30, outfile=outprefix)
    finally:
        sys.stderr = old_stderr
    seq_file = Path(outprefix + ".seq.DupRate.xls")
    pos_file = Path(outprefix + ".pos.DupRate.xls")
    assert seq_file.exists()
    assert pos_file.exists()
    seq_content = seq_file.read_text()
    assert "Occurrence" in seq_content
    pos_content = pos_file.read_text()
    assert "Occurrence" in pos_content

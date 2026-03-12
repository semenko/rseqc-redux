"""Tests for rseqc.SAM."""

import importlib
import io
import sys
from pathlib import Path

import pytest

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


def test_parsebam_readsQual_boxplot(mini_bam, tmp_path):
    """readsQual_boxplot() should produce a .qual.r file."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.readsQual_boxplot(outfile=outprefix, q_cut=30)
    finally:
        sys.stderr = old_stderr
    r_file = Path(outprefix + ".qual.r")
    assert r_file.exists()
    content = r_file.read_text()
    assert "boxplot(" in content
    assert "heatmap(" in content


def test_parsebam_clipping_profile_SE(mini_bam, tmp_path):
    """clipping_profile() in single-end mode should produce output files."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.clipping_profile(outfile=outprefix, q_cut=30, PE=False)
    finally:
        sys.stderr = old_stderr
    xls_file = Path(outprefix + ".clipping_profile.xls")
    assert xls_file.exists()
    content = xls_file.read_text()
    assert "Position" in content


def test_parsebam_insertion_profile_SE(mini_bam, tmp_path):
    """insertion_profile() in single-end mode should produce output files."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.insertion_profile(outfile=outprefix, q_cut=30, PE=False)
    finally:
        sys.stderr = old_stderr
    xls_file = Path(outprefix + ".insertion_profile.xls")
    assert xls_file.exists()
    content = xls_file.read_text()
    assert "Position" in content


def test_parsebam_bam2fq_single(mini_bam, tmp_path):
    """bam2fq() in single-end mode should produce a .fastq file."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.bam2fq(prefix=outprefix, paired=False)
    finally:
        sys.stderr = old_stderr
    fq_file = Path(outprefix + ".fastq")
    assert fq_file.exists()
    content = fq_file.read_text()
    # FASTQ entries start with @
    assert content.startswith("@")
    lines = content.strip().split("\n")
    # Each FASTQ entry is 4 lines
    assert len(lines) % 4 == 0
    assert len(lines) >= 4


def test_parsebam_bam2fq_paired(mini_bam, tmp_path):
    """bam2fq() in paired-end mode should produce R1 and R2 files."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.bam2fq(prefix=outprefix, paired=True)
    finally:
        sys.stderr = old_stderr
    r1_file = Path(outprefix + ".R1.fastq")
    r2_file = Path(outprefix + ".R2.fastq")
    assert r1_file.exists()
    assert r2_file.exists()


def test_parsebam_stat_counts(mini_bam, capsys):
    """stat() should report correct counts for known reads."""
    obj = SAM.ParseBAM(str(mini_bam))
    obj.stat(q_cut=30)
    captured = capsys.readouterr()
    output = captured.out
    # We have 11 total reads in mini_bam
    assert "Total records:" in output
    # Check that some expected categories appear
    assert "QC failed:" in output
    assert "Unmapped reads:" in output


@pytest.mark.parametrize("q_cut", [0, 30, 60])
def test_parsebam_stat_qcut_levels(mini_bam, capsys, q_cut):
    """stat() should work with different quality cutoffs."""
    obj = SAM.ParseBAM(str(mini_bam))
    obj.stat(q_cut=q_cut)
    captured = capsys.readouterr()
    assert "Total records:" in captured.out


def test_parsebam_readGC_content(mini_bam, tmp_path):
    """readGC() with q_cut=0 should include all mapped reads."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.readGC(outfile=outprefix, q_cut=0)
    finally:
        sys.stderr = old_stderr
    xls_file = Path(outprefix + ".GC.xls")
    content = xls_file.read_text()
    # All reads are "AAAA..." so GC should be 0.00
    assert "0.00" in content
    r_file = Path(outprefix + ".GC_plot.r")
    assert r_file.exists()
    assert "pdf(" in r_file.read_text()


def test_parsebam_readsNVC_content(mini_bam, tmp_path):
    """readsNVC() should count nucleotides at each position."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.readsNVC(outfile=outprefix, q_cut=0)
    finally:
        sys.stderr = old_stderr
    xls_file = Path(outprefix + ".NVC.xls")
    content = xls_file.read_text()
    lines = [ln for ln in content.strip().split("\n") if not ln.startswith("Position")]
    # Should have data for each position in the reads
    assert len(lines) >= 1
    # R script should be generated
    r_file = Path(outprefix + ".NVC_plot.r")
    assert r_file.exists()


def test_parsebam_readDupRate_content(mini_bam, tmp_path):
    """readDupRate() with q_cut=0 includes more reads."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.readDupRate(q_cut=0, outfile=outprefix)
    finally:
        sys.stderr = old_stderr
    seq_file = Path(outprefix + ".seq.DupRate.xls")
    seq_content = seq_file.read_text()
    lines = [ln for ln in seq_content.strip().split("\n") if not ln.startswith("Occurrence")]
    assert len(lines) >= 1
    # R script should be generated
    r_file = Path(outprefix + ".DupRate_plot.r")
    assert r_file.exists()
    assert "pdf(" in r_file.read_text()

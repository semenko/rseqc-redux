"""Integration tests: run actual CLI commands with mini BAM + mini BED."""

import subprocess
import sys
from pathlib import Path

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def test_bam_stat(mini_bam):
    """Run bam_stat with mini BAM file."""
    result = subprocess.run(
        [sys.executable, "-m", "scripts.bam_stat", "-i", str(mini_bam)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, f"bam_stat failed: {result.stderr}"
    output = result.stdout + result.stderr
    assert "Total records:" in output or "mapq" in output.lower() or "Total" in output


def test_infer_experiment(mini_bam):
    """Run infer_experiment with mini BAM + mini BED."""
    bed_file = str(FIXTURES_DIR / "mini.bed")
    result = subprocess.run(
        [sys.executable, "-m", "scripts.infer_experiment", "-i", str(mini_bam), "-r", bed_file],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, f"infer_experiment failed: {result.stderr}"
    output = result.stdout + result.stderr
    # Should report sample info
    assert "usable reads were sampled" in output or "Fraction" in output


def test_read_distribution(mini_bam):
    """Run read_distribution with mini BAM + mini BED."""
    bed_file = str(FIXTURES_DIR / "mini.bed")
    result = subprocess.run(
        [sys.executable, "-m", "scripts.read_distribution", "-i", str(mini_bam), "-r", bed_file],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, f"read_distribution failed: {result.stderr}"
    output = result.stdout + result.stderr
    assert "Total Assigned Tags" in output or "Group" in output or "Tags" in output


def test_read_GC(mini_bam, tmp_path):
    """Run read_GC with mini BAM file."""
    outprefix = str(tmp_path / "gc_test")
    result = subprocess.run(
        [sys.executable, "-m", "scripts.read_GC", "-i", str(mini_bam), "-o", outprefix],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, f"read_GC failed: {result.stderr}"


def test_read_quality(mini_bam, tmp_path):
    """Run read_quality with mini BAM file."""
    outprefix = str(tmp_path / "qual_test")
    result = subprocess.run(
        [sys.executable, "-m", "scripts.read_quality", "-i", str(mini_bam), "-o", outprefix],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, f"read_quality failed: {result.stderr}"


def test_junction_annotation(mini_bam, tmp_path):
    """Run junction_annotation with mini BAM + mini BED."""
    bed_file = str(FIXTURES_DIR / "mini.bed")
    outprefix = str(tmp_path / "junc_anno")
    result = subprocess.run(
        [
            sys.executable, "-m", "scripts.junction_annotation",
            "-i", str(mini_bam),
            "-r", bed_file,
            "-o", outprefix,
            "-q", "0",
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, f"junction_annotation failed: {result.stderr}"


def test_junction_saturation(mini_bam, tmp_path):
    """Run junction_saturation with mini BAM + mini BED."""
    bed_file = str(FIXTURES_DIR / "mini.bed")
    outprefix = str(tmp_path / "junc_sat")
    result = subprocess.run(
        [
            sys.executable, "-m", "scripts.junction_saturation",
            "-i", str(mini_bam),
            "-r", bed_file,
            "-o", outprefix,
            "-q", "0",
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, f"junction_saturation failed: {result.stderr}"


def test_inner_distance(mini_bam, tmp_path):
    """Run inner_distance with mini BAM + mini BED."""
    bed_file = str(FIXTURES_DIR / "mini.bed")
    outprefix = str(tmp_path / "inner_dist")
    result = subprocess.run(
        [
            sys.executable, "-m", "scripts.inner_distance",
            "-i", str(mini_bam),
            "-r", bed_file,
            "-o", outprefix,
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, f"inner_distance failed: {result.stderr}"


def test_mismatch_profile(mini_bam, tmp_path):
    """Run mismatch_profile with mini BAM file.

    The mini BAM has no MD tags and uniform bases, so the script will find
    no mismatches and exit with code 1.  We accept that as valid behavior.
    """
    outprefix = str(tmp_path / "mismatch")
    result = subprocess.run(
        [
            sys.executable, "-m", "scripts.mismatch_profile",
            "-i", str(mini_bam),
            "-l", "50",
            "-o", outprefix,
            "-q", "0",
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )
    # returncode 0 = normal, returncode 1 with "No mismatches found" = expected for mini BAM
    if result.returncode != 0:
        assert "No mismatches found" in result.stderr, f"mismatch_profile failed unexpectedly: {result.stderr}"


def test_deletion_profile(mini_bam, tmp_path):
    """Run deletion_profile with mini BAM file."""
    outprefix = str(tmp_path / "deletion")
    result = subprocess.run(
        [
            sys.executable, "-m", "scripts.deletion_profile",
            "-i", str(mini_bam),
            "-l", "50",
            "-o", outprefix,
            "-q", "0",
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, f"deletion_profile failed: {result.stderr}"

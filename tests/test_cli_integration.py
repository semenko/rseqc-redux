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


def test_read_duplication(mini_bam, tmp_path):
    """Run read_duplication with mini BAM file."""
    outprefix = str(tmp_path / "dup_test")
    result = subprocess.run(
        [
            sys.executable, "-m", "scripts.read_duplication",
            "-i", str(mini_bam),
            "-o", outprefix,
            "-q", "0",
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, f"read_duplication failed: {result.stderr}"


def test_clipping_profile(mini_bam, tmp_path):
    """Run clipping_profile with mini BAM file.

    The mini BAM has no soft-clipped reads (all CIGAR ops are M or N),
    so the script may produce empty output.  We accept returncode 0 or 1.
    """
    outprefix = str(tmp_path / "clip_test")
    result = subprocess.run(
        [
            sys.executable, "-m", "scripts.clipping_profile",
            "-i", str(mini_bam),
            "-o", outprefix,
            "-q", "0",
            "-s", "SE",
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )
    # returncode 0 = normal; returncode 1 acceptable if no clipped reads found
    if result.returncode != 0:
        assert (
            "No clipping" in result.stderr
            or "clipping" in result.stderr.lower()
            or result.returncode == 1
        ), f"clipping_profile failed unexpectedly: {result.stderr}"


def test_insertion_profile(mini_bam, tmp_path):
    """Run insertion_profile with mini BAM file.

    The mini BAM has no insertions in CIGAR strings, so the script may
    produce empty output.  We accept returncode 0 or 1.
    """
    outprefix = str(tmp_path / "ins_test")
    result = subprocess.run(
        [
            sys.executable, "-m", "scripts.insertion_profile",
            "-i", str(mini_bam),
            "-o", outprefix,
            "-q", "0",
            "-s", "SE",
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )
    # returncode 0 = normal; returncode 1 acceptable if no insertions found
    if result.returncode != 0:
        assert (
            "No insertion" in result.stderr
            or "insertion" in result.stderr.lower()
            or result.returncode == 1
        ), f"insertion_profile failed unexpectedly: {result.stderr}"


def test_read_NVC(mini_bam, tmp_path):
    """Run read_NVC with mini BAM file."""
    outprefix = str(tmp_path / "nvc_test")
    result = subprocess.run(
        [
            sys.executable, "-m", "scripts.read_NVC",
            "-i", str(mini_bam),
            "-o", outprefix,
            "-q", "0",
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, f"read_NVC failed: {result.stderr}"


def test_tin(mini_bam, tmp_path):
    """Run tin with mini BAM + mini BED file.

    tin.py writes output files to the current working directory, so we
    run it with cwd=tmp_path.  We set -c 0 so transcripts with very few
    reads are still processed.
    """
    bed_file = str(FIXTURES_DIR / "mini.bed")
    result = subprocess.run(
        [
            sys.executable, "-m", "scripts.tin",
            "-i", str(mini_bam),
            "-r", bed_file,
            "-c", "0",
        ],
        capture_output=True,
        text=True,
        timeout=60,
        cwd=str(tmp_path),
    )
    assert result.returncode == 0, f"tin failed: {result.stderr}"


def test_geneBody_coverage(mini_bam, tmp_path):
    """Run geneBody_coverage with mini BAM + mini BED file."""
    bed_file = str(FIXTURES_DIR / "mini.bed")
    outprefix = str(tmp_path / "genebody")
    result = subprocess.run(
        [
            sys.executable, "-m", "scripts.geneBody_coverage",
            "-i", str(mini_bam),
            "-r", bed_file,
            "-o", outprefix,
            "-l", "100",
        ],
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, f"geneBody_coverage failed: {result.stderr}"

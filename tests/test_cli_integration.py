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

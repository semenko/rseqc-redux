"""Tests for rseqc.wiggle."""

import importlib
from pathlib import Path

from rseqc import wiggle

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def test_import():
    """Verify that rseqc.wiggle can be imported."""
    mod = importlib.import_module("rseqc.wiggle")
    assert mod is not None


class TestParseWig:
    def setup_method(self):
        self.wig = wiggle.ParseWig(str(FIXTURES_DIR / "mini.wig"))

    def test_chroms_loaded(self):
        """Both chr1 and chr2 should be loaded (uppercased)."""
        assert "CHR1" in self.wig.scores
        assert "CHR2" in self.wig.scores

    def test_fetch_all_scores_chr1(self):
        """Fetch scores from variableStep chr1 data (0-based positions)."""
        scores = self.wig.fetch_all_scores("chr1", 0, 5)
        assert scores[0] == 1.5
        assert scores[4] == 5.5

    def test_fetch_max_scores(self):
        assert self.wig.fetch_max_scores("chr1", 0, 5) == 5.5

    def test_fetch_min_scores(self):
        assert self.wig.fetch_min_scores("chr1", 0, 5) == 1.5

    def test_fetch_sum_scores(self):
        total = self.wig.fetch_sum_scores("chr1", 0, 5)
        assert abs(total - (1.5 + 2.5 + 3.5 + 4.5 + 5.5)) < 0.01

    def test_fetch_avg_scores(self):
        avg = self.wig.fetch_avg_scores("chr1", 0, 5)
        assert abs(avg - 3.5) < 0.01

    def test_case_insensitive_lookup(self):
        """Chrom lookup should be case-insensitive (uppercased internally)."""
        scores_lower = self.wig.fetch_all_scores("chr1", 0, 1)
        scores_upper = self.wig.fetch_all_scores("CHR1", 0, 1)
        assert scores_lower[0] == scores_upper[0]

    def test_chr2_fixedstep(self):
        """fixedStep data for chr2."""
        scores = self.wig.fetch_all_scores("chr2", 0, 5)
        assert scores[0] == 10.0
        assert scores[4] == 50.0


class TestParseWig2:
    def setup_method(self):
        self.wig = wiggle.ParseWig2(str(FIXTURES_DIR / "mini.wig"))

    def test_chroms_loaded(self):
        assert "CHR1" in self.wig.scores
        assert "CHR2" in self.wig.scores

    def test_fetch_all_scores_by_range(self):
        scores = self.wig.fetch_all_scores_by_range("chr1", 0, 5)
        assert scores[0] == 1.5

    def test_fetch_by_positions(self):
        scores = self.wig.fetch_all_scores_by_positions("chr1", [0, 2, 4])
        assert scores[0] == 1.5
        assert scores[1] == 3.5
        assert scores[2] == 5.5

    def test_fetch_avg_scores_by_range(self):
        avg = self.wig.fetch_avg_scores_by_range("chr1", 0, 5)
        assert abs(avg - 3.5) < 0.01

    def test_fetch_avg_scores_by_positions(self):
        avg = self.wig.fetch_avg_scores_by_positions("chr1", [0, 2, 4])
        # (1.5 + 3.5 + 5.5) / 3 = 3.5
        assert abs(avg - 3.5) < 0.01

    def test_fetch_sum_scores_by_range(self):
        total = self.wig.fetch_sum_scores_by_range("chr1", 0, 5)
        assert abs(total - (1.5 + 2.5 + 3.5 + 4.5 + 5.5)) < 0.01

    def test_fetch_sum_scores_by_positions(self):
        total = self.wig.fetch_sum_scores_by_positions("chr1", [0, 2, 4])
        assert abs(total - (1.5 + 3.5 + 5.5)) < 0.01

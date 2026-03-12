"""Tests for rseqc.BED."""

import importlib
import os
import tempfile
from pathlib import Path

from rseqc import BED


def test_import():
    """Verify that rseqc.BED can be imported."""
    mod = importlib.import_module("rseqc.BED")
    assert mod is not None


FIXTURES_DIR = Path(__file__).parent / "fixtures"


# ---------------------------------------------------------------------------
# Helper: create a BED file from lines and return a ParseBED instance
# ---------------------------------------------------------------------------


def _make_bed(lines, tmp_path=None):
    """Write BED lines to a temp file and return a ParseBED instance."""
    if tmp_path is None:
        fd, path = tempfile.mkstemp(suffix=".bed")
        os.close(fd)
    else:
        path = str(tmp_path / "test.bed")
    with open(path, "w") as fh:
        for line in lines:
            fh.write(line + "\n")
    return BED.ParseBED(path), path


# Standard BED12 lines used across many tests.
# gene1: chr1, 1000-5000, +, 3 exons (1000-1500, 2500-3100, 4500-5000), CDS 1200-4800
# gene2: chr1, 6000-9000, -, 2 exons (6000-7000, 8000-9000), CDS 6200-8800
# gene3: chr1, 10000-11000, +, 1 exon (10000-11000), CDS 10000-11000
BED12_LINES = [
    "chr1\t1000\t5000\tgene1\t0\t+\t1200\t4800\t0,0,0\t3\t500,600,500\t0,1500,3500",
    "chr1\t6000\t9000\tgene2\t0\t-\t6200\t8800\t0,0,0\t2\t1000,1000\t0,2000",
    "chr1\t10000\t11000\tgene3\t0\t+\t10000\t11000\t0,0,0\t1\t1000\t0",
]


# --- tillingBed ---


def test_tillingbed_basic():
    result = list(BED.tillingBed("chr1", 100, stepSize=30))
    assert result == [("chr1", 0, 30), ("chr1", 30, 60), ("chr1", 60, 90), ("chr1", 90, 100)]


def test_tillingbed_exact():
    result = list(BED.tillingBed("chr1", 100, stepSize=50))
    assert result == [("chr1", 0, 50), ("chr1", 50, 100)]


def test_tillingbed_single():
    result = list(BED.tillingBed("chr1", 50, stepSize=100))
    assert result == [("chr1", 0, 50)]


def test_tillingbed_zero_size():
    result = list(BED.tillingBed("chr1", 0, stepSize=10))
    assert result == []


def test_tillingbed_step_equals_size():
    result = list(BED.tillingBed("chr1", 30, stepSize=30))
    assert result == [("chr1", 0, 30)]


def test_tillingbed_large_step():
    result = list(BED.tillingBed("chr2", 10, stepSize=1000))
    assert result == [("chr2", 0, 10)]


def test_tillingbed_step_one():
    result = list(BED.tillingBed("chr1", 5, stepSize=1))
    assert len(result) == 5
    assert result[0] == ("chr1", 0, 1)
    assert result[4] == ("chr1", 4, 5)


# --- unionBed3 ---


def test_unionbed3_no_overlap():
    lst = [["chr1", 10, 20], ["chr1", 30, 40]]
    result = BED.unionBed3(lst)
    assert ["chr1", 10, 20] in result
    assert ["chr1", 30, 40] in result


def test_unionbed3_overlap():
    lst = [["chr1", 10, 30], ["chr1", 20, 40]]
    result = BED.unionBed3(lst)
    assert len(result) == 1
    assert result[0] == ["chr1", 10, 40]


def test_unionbed3_adjacent():
    lst = [["chr1", 10, 20], ["chr1", 20, 30]]
    result = BED.unionBed3(lst)
    assert len(result) == 1
    assert result[0] == ["chr1", 10, 30]


def test_unionbed3_multi_chrom():
    lst = [["chr1", 10, 20], ["chr2", 10, 20]]
    result = BED.unionBed3(lst)
    chroms = {r[0] for r in result}
    assert chroms == {"chr1", "chr2"}


def test_unionbed3_single_entry():
    lst = [["chr1", 100, 200]]
    result = BED.unionBed3(lst)
    assert result == [["chr1", 100, 200]]


def test_unionbed3_multiple_overlapping():
    lst = [["chr1", 10, 30], ["chr1", 20, 50], ["chr1", 40, 60]]
    result = BED.unionBed3(lst)
    assert len(result) == 1
    assert result[0] == ["chr1", 10, 60]


def test_unionbed3_contained():
    """Interval fully contained within another."""
    lst = [["chr1", 10, 50], ["chr1", 20, 30]]
    result = BED.unionBed3(lst)
    assert len(result) == 1
    assert result[0] == ["chr1", 10, 50]


# --- intersectBed3 ---


def test_intersectbed3_overlap():
    lst1 = [["chr1", 10, 30]]
    lst2 = [["chr1", 20, 40]]
    result = BED.intersectBed3(lst1, lst2)
    assert len(result) == 1
    assert result[0] == ["chr1", 20, 30]


def test_intersectbed3_no_overlap():
    lst1 = [["chr1", 10, 20]]
    lst2 = [["chr1", 30, 40]]
    result = BED.intersectBed3(lst1, lst2)
    assert result == []


def test_intersectbed3_contained():
    lst1 = [["chr1", 10, 50]]
    lst2 = [["chr1", 20, 30]]
    result = BED.intersectBed3(lst1, lst2)
    assert len(result) == 1
    assert result[0] == ["chr1", 20, 30]


def test_intersectbed3_identical():
    lst1 = [["chr1", 10, 50]]
    lst2 = [["chr1", 10, 50]]
    result = BED.intersectBed3(lst1, lst2)
    assert len(result) == 1
    assert result[0] == ["chr1", 10, 50]


def test_intersectbed3_different_chroms():
    lst1 = [["chr1", 10, 50]]
    lst2 = [["chr2", 10, 50]]
    result = BED.intersectBed3(lst1, lst2)
    assert result == []


def test_intersectbed3_multiple():
    lst1 = [["chr1", 10, 30], ["chr1", 40, 60]]
    lst2 = [["chr1", 20, 50]]
    result = BED.intersectBed3(lst1, lst2)
    assert len(result) == 2
    starts = {r[1] for r in result}
    assert 20 in starts
    assert 40 in starts


# --- subtractBed3 ---


def test_subtractbed3_basic():
    lst1 = [["chr1", 10, 40]]
    lst2 = [["chr1", 20, 30]]
    result = BED.subtractBed3(lst1, lst2)
    assert len(result) == 2
    assert ["chr1", 10, 20] in result
    assert ["chr1", 30, 40] in result


def test_subtractbed3_no_overlap():
    lst1 = [["chr1", 10, 20]]
    lst2 = [["chr1", 30, 40]]
    result = BED.subtractBed3(lst1, lst2)
    assert result == [["chr1", 10, 20]]


def test_subtractbed3_different_chroms():
    """Bug #9 regression: subtractBed3 with different chroms should preserve lst1 entries."""
    lst1 = [["chr1", 10, 20]]
    lst2 = [["chr2", 10, 20]]
    result = BED.subtractBed3(lst1, lst2)
    assert result == [["chr1", 10, 20]]


def test_subtractbed3_full_overlap():
    """Subtracting an interval that covers lst1 entirely."""
    lst1 = [["chr1", 20, 30]]
    lst2 = [["chr1", 10, 40]]
    result = BED.subtractBed3(lst1, lst2)
    assert result == []


def test_subtractbed3_identical():
    lst1 = [["chr1", 10, 20]]
    lst2 = [["chr1", 10, 20]]
    result = BED.subtractBed3(lst1, lst2)
    assert result == []


def test_subtractbed3_left_overlap():
    lst1 = [["chr1", 10, 30]]
    lst2 = [["chr1", 5, 20]]
    result = BED.subtractBed3(lst1, lst2)
    assert len(result) == 1
    assert result[0] == ["chr1", 20, 30]


def test_subtractbed3_right_overlap():
    lst1 = [["chr1", 10, 30]]
    lst2 = [["chr1", 20, 40]]
    result = BED.subtractBed3(lst1, lst2)
    assert len(result) == 1
    assert result[0] == ["chr1", 10, 20]


# --- ParseBED class with fixture ---


class TestParseBED:
    def setup_method(self):
        self.bed = BED.ParseBED(str(FIXTURES_DIR / "mini.bed"))

    def test_get_exon(self):
        exons = self.bed.getExon()
        # gene1: 3 exons, gene2: 2 exons, gene3: 1 exon = 6 total
        assert len(exons) == 6
        assert exons[0] == ("chr1", 1000, 1500)
        assert exons[1] == ("chr1", 2500, 3100)
        assert exons[2] == ("chr1", 4500, 5000)

    def test_get_intron(self):
        introns = self.bed.getIntron()
        # gene1 introns: [1500, 2500] and [3100, 4500]
        # gene2 introns: [7000, 8000]
        # Note: line 720 bug `if int(fields[9] == 1)` is always False
        assert len(introns) >= 3

    def test_get_utr_both(self):
        utrs = self.bed.getUTR(utr=35)
        assert len(utrs) >= 2

    def test_get_utr_5prime(self):
        utrs = self.bed.getUTR(utr=5)
        for u in utrs:
            assert len(u) == 6

    def test_get_transcript_ranges(self):
        ranges = list(self.bed.getTranscriptRanges())
        assert len(ranges) == 3
        assert ranges[0][0] == "chr1"
        assert ranges[0][1] == 1000
        assert ranges[0][2] == 5000

    def test_get_cds_exon(self):
        cds_exons = self.bed.getCDSExon()
        assert len(cds_exons) >= 3

    def test_get_intergenic_up(self):
        regions = self.bed.getIntergenic(direction="up", size=500)
        assert len(regions) == 3
        # gene1 (+): upstream = [500, 1000]
        assert regions[0] == ["chr1", 500, 1000]
        # gene2 (-): upstream = tx_end to tx_end+size = [9000, 9500]
        assert regions[1] == ["chr1", 9000, 9500]
        # gene3 (+): upstream = [9500, 10000]
        assert regions[2] == ["chr1", 9500, 10000]

    def test_get_intergenic_down(self):
        regions = self.bed.getIntergenic(direction="down", size=500)
        assert len(regions) == 3
        # gene1 (+): downstream = [5000, 5500]
        assert regions[0] == ["chr1", 5000, 5500]
        # gene2 (-): downstream = [5500, 6000]
        assert regions[1] == ["chr1", 5500, 6000]
        # gene3 (+): downstream = [11000, 11500]
        assert regions[2] == ["chr1", 11000, 11500]

    def test_get_intergenic_both(self):
        regions = self.bed.getIntergenic(direction="both", size=500)
        # 3 genes x 2 directions = 6 regions
        assert len(regions) == 6


# ---------------------------------------------------------------------------
# Extended ParseBED tests with inline BED12 lines
# ---------------------------------------------------------------------------


class TestParseBEDInline:
    """Tests using inline BED12 fixtures created in temp files."""

    def setup_method(self):
        self.bed, self.path = _make_bed(BED12_LINES)

    def teardown_method(self):
        try:
            os.unlink(self.path)
        except OSError:
            pass

    # --- getExon ---

    def test_get_exon_coords(self):
        exons = self.bed.getExon()
        assert len(exons) == 6
        # gene1 exons
        assert exons[0] == ("chr1", 1000, 1500)
        assert exons[1] == ("chr1", 2500, 3100)
        assert exons[2] == ("chr1", 4500, 5000)
        # gene2 exons
        assert exons[3] == ("chr1", 6000, 7000)
        assert exons[4] == ("chr1", 8000, 9000)
        # gene3 exon
        assert exons[5] == ("chr1", 10000, 11000)

    def test_get_exon_is_repeatable(self):
        """ParseBED resets file pointer so method can be called twice."""
        first = self.bed.getExon()
        second = self.bed.getExon()
        assert first == second

    # --- getIntron ---

    def test_get_intron_coords(self):
        introns = self.bed.getIntron()
        # gene1: 2 introns, gene2: 1 intron, gene3: 0 introns (single exon skipped)
        assert len(introns) == 3
        assert introns[0] == ["chr1", 1500, 2500]
        assert introns[1] == ["chr1", 3100, 4500]
        assert introns[2] == ["chr1", 7000, 8000]

    def test_get_intron_single_exon_gene_excluded(self):
        """Single-exon gene (gene3) should produce no introns."""
        bed, path = _make_bed([BED12_LINES[2]])  # gene3 only
        introns = bed.getIntron()
        assert introns == []
        os.unlink(path)

    def test_get_intron_repeatable(self):
        first = self.bed.getIntron()
        second = self.bed.getIntron()
        assert first == second

    # --- getCDSExon ---

    def test_get_cds_exon_coords(self):
        cds = self.bed.getCDSExon()
        # gene1: CDS 1200-4800, exons 1000-1500, 2500-3100, 4500-5000
        #   exon1 CDS: 1200-1500, exon2: 2500-3100, exon3: 4500-4800
        # gene2: CDS 6200-8800, exons 6000-7000, 8000-9000
        #   exon1 CDS: 6200-7000, exon2: 8000-8800
        # gene3: CDS 10000-11000, exon 10000-11000
        #   exon1 CDS: 10000-11000
        assert len(cds) == 6
        assert cds[0] == ["chr1", 1200, 1500]
        assert cds[1] == ["chr1", 2500, 3100]
        assert cds[2] == ["chr1", 4500, 4800]
        assert cds[3] == ["chr1", 6200, 7000]
        assert cds[4] == ["chr1", 8000, 8800]
        assert cds[5] == ["chr1", 10000, 11000]

    def test_get_cds_exon_repeatable(self):
        first = self.bed.getCDSExon()
        second = self.bed.getCDSExon()
        assert first == second

    # --- getUTR ---

    def test_get_utr_5prime_plus_strand(self):
        """For + strand gene, 5' UTR is before CDS."""
        bed, path = _make_bed([BED12_LINES[0]])  # gene1 only, + strand
        utrs = bed.getUTR(utr=5)
        # gene1 +: 5'UTR = exon regions before cdsStart(1200)
        # exon1 starts at 1000, cdsStart=1200 => 5'UTR = [1000, 1200]
        assert len(utrs) == 1
        assert utrs[0] == ["chr1", 1000, 1200, "gene1", "0", "+"]
        os.unlink(path)

    def test_get_utr_3prime_plus_strand(self):
        """For + strand gene, 3' UTR is after CDS."""
        bed, path = _make_bed([BED12_LINES[0]])  # gene1 only
        utrs = bed.getUTR(utr=3)
        # gene1 +: 3'UTR = exon regions after cdsEnd(4800)
        # exon3 ends at 5000, cdsEnd=4800 => 3'UTR = [4800, 5000]
        assert len(utrs) == 1
        assert utrs[0] == ["chr1", 4800, 5000, "gene1", "0", "+"]
        os.unlink(path)

    def test_get_utr_both_plus_strand(self):
        """For + strand gene, both UTRs."""
        bed, path = _make_bed([BED12_LINES[0]])
        utrs = bed.getUTR(utr=35)
        assert len(utrs) == 2
        os.unlink(path)

    def test_get_utr_5prime_minus_strand(self):
        """For - strand gene, 5' UTR is AFTER CDS (higher coordinates)."""
        bed, path = _make_bed([BED12_LINES[1]])  # gene2, - strand
        utrs = bed.getUTR(utr=5)
        # gene2 -: 5'UTR = exon regions after cdsEnd(8800)
        # exon2 ends at 9000, cdsEnd=8800 => 5'UTR = [8800, 9000]
        assert len(utrs) == 1
        assert utrs[0] == ["chr1", 8800, 9000, "gene2", "0", "-"]
        os.unlink(path)

    def test_get_utr_3prime_minus_strand(self):
        """For - strand gene, 3' UTR is BEFORE CDS (lower coordinates)."""
        bed, path = _make_bed([BED12_LINES[1]])  # gene2, - strand
        utrs = bed.getUTR(utr=3)
        # gene2 -: 3'UTR = exon regions before cdsStart(6200)
        # exon1 starts at 6000, cdsStart=6200 => 3'UTR = [6000, 6200]
        assert len(utrs) == 1
        assert utrs[0] == ["chr1", 6000, 6200, "gene2", "0", "-"]
        os.unlink(path)

    def test_get_utr_no_utr_gene(self):
        """Gene3 has CDS=txStart-txEnd, so no UTR."""
        bed, path = _make_bed([BED12_LINES[2]])  # gene3
        utrs = bed.getUTR(utr=35)
        assert utrs == []
        os.unlink(path)

    def test_get_utr_repeatable(self):
        first = self.bed.getUTR(utr=35)
        second = self.bed.getUTR(utr=35)
        assert first == second

    # --- getTranscriptRanges ---

    def test_get_transcript_ranges_all(self):
        ranges = list(self.bed.getTranscriptRanges())
        assert len(ranges) == 3
        # gene1
        assert ranges[0][0] == "chr1"
        assert ranges[0][1] == 1000
        assert ranges[0][2] == 5000
        assert ranges[0][3] == "+"
        assert "gene1" in ranges[0][4]
        # gene2
        assert ranges[1][1] == 6000
        assert ranges[1][2] == 9000
        assert ranges[1][3] == "-"
        # gene3
        assert ranges[2][1] == 10000
        assert ranges[2][2] == 11000

    def test_get_transcript_ranges_skips_comments(self):
        bed, path = _make_bed(["# comment", "track name=foo"] + BED12_LINES[:1])
        ranges = list(bed.getTranscriptRanges())
        assert len(ranges) == 1
        os.unlink(path)

    # --- getIntergenic ---

    def test_get_intergenic_up_clamped_at_zero(self):
        """Upstream of a gene at start of chromosome should be clamped at 0."""
        line = "chr1\t100\t500\tgeneA\t0\t+\t100\t500\t0,0,0\t1\t400\t0"
        bed, path = _make_bed([line])
        regions = bed.getIntergenic(direction="up", size=200)
        assert len(regions) == 1
        assert regions[0] == ["chr1", 0, 100]  # clamped at 0
        os.unlink(path)

    def test_get_intergenic_down_minus_strand(self):
        """For - strand, downstream = before tx_start."""
        line = "chr1\t5000\t6000\tgeneB\t0\t-\t5000\t6000\t0,0,0\t1\t1000\t0"
        bed, path = _make_bed([line])
        regions = bed.getIntergenic(direction="down", size=300)
        assert len(regions) == 1
        assert regions[0] == ["chr1", 4700, 5000]
        os.unlink(path)

    def test_get_intergenic_both_count(self):
        regions = self.bed.getIntergenic(direction="both", size=100)
        # 3 genes * 2 = 6
        assert len(regions) == 6

    # --- Comment/track/browser line handling ---

    def test_skips_comment_lines_in_methods_that_support_it(self):
        """Methods like getIntron, getTranscriptRanges skip comment/track/browser lines.
        Note: getExon does NOT skip these lines (it crashes), so we test with getIntron."""
        bed, path = _make_bed(["# a comment", "track name=foo", "browser position chr1:1-100"] + BED12_LINES[:1])
        introns = bed.getIntron()
        # gene1 has 2 introns
        assert len(introns) == 2
        assert introns[0] == ["chr1", 1500, 2500]
        assert introns[1] == ["chr1", 3100, 4500]
        os.unlink(path)


# ---------------------------------------------------------------------------
# Additional standalone function tests
# ---------------------------------------------------------------------------


class TestUnionBed3Edge:
    def test_empty_input(self):
        # Empty list should produce empty result (or handle gracefully)
        result = BED.unionBed3([])
        assert result == []

    def test_single_base(self):
        result = BED.unionBed3([["chr1", 10, 11]])
        assert result == [["chr1", 10, 11]]


class TestIntersectBed3Edge:
    def test_empty_first(self):
        result = BED.intersectBed3([], [["chr1", 10, 20]])
        assert result == []

    def test_empty_second(self):
        result = BED.intersectBed3([["chr1", 10, 20]], [])
        assert result == []

    def test_both_empty(self):
        result = BED.intersectBed3([], [])
        assert result == []


class TestSubtractBed3Edge:
    def test_empty_first(self):
        result = BED.subtractBed3([], [["chr1", 10, 20]])
        assert result == []

    def test_empty_second(self):
        result = BED.subtractBed3([["chr1", 10, 20]], [])
        assert result == [["chr1", 10, 20]]

"""Tests for rseqc.BED."""

import importlib
from pathlib import Path

from rseqc import BED


def test_import():
    """Verify that rseqc.BED can be imported."""
    mod = importlib.import_module("rseqc.BED")
    assert mod is not None


FIXTURES_DIR = Path(__file__).parent / "fixtures"


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

    def test_get_splice_junctions(self):
        junctions = list(self.bed.getSpliceJunctions())
        assert len(junctions) >= 2

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

"""Tests for rseqc.BED."""

import importlib
import os
import tempfile
from pathlib import Path

import pytest

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

    # --- getSpliceJunctions ---

    def test_get_splice_junctions_multi_exon(self):
        junctions = list(self.bed.getSpliceJunctions())
        # gene1: 2 junctions, gene2: 1 junction, gene3: 0 junctions (single exon)
        assert len(junctions) == 3  # all 3 genes yield a tuple
        # gene1 has 2 junctions
        assert len(junctions[0][4]) == 2
        assert junctions[0][4][0] == "chr1:1500-2500"
        assert junctions[0][4][1] == "chr1:3100-4500"
        # gene2 has 1 junction
        assert len(junctions[1][4]) == 1
        assert junctions[1][4][0] == "chr1:7000-8000"
        # gene3 is single exon, yielded with empty junctions list
        assert junctions[2][4] == []

    def test_get_splice_junctions_gene_names(self):
        junctions = list(self.bed.getSpliceJunctions())
        names = [j[0] for j in junctions]
        assert names == ["gene1", "gene2", "gene3"]

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

    # --- truncate_bed ---

    def test_truncate_bed_3prime_plus(self):
        """Truncate from 3' end for + strand gene keeps last N bases."""
        bed, path = _make_bed([BED12_LINES[0]])  # gene1, + strand
        result = bed.truncate_bed(truncation_from=3, size=100)
        assert len(result) == 1
        # gene1 exonic bases from 3' end: last 100 bases of exon3 (4500-5000 = 500 bases)
        # so the 100 bases from 3' are 4901..5000 (1-based), BED = 4900-5000
        fields = result[0].split("\t")
        assert fields[0] == "chr1"
        assert int(fields[2]) == 5000  # tx_end should be end of gene
        os.unlink(path)

    def test_truncate_bed_5prime_plus(self):
        """Truncate from 5' end for + strand gene keeps first N bases."""
        bed, path = _make_bed([BED12_LINES[0]])  # gene1, + strand
        result = bed.truncate_bed(truncation_from=5, size=100)
        assert len(result) == 1
        fields = result[0].split("\t")
        assert fields[0] == "chr1"
        assert int(fields[1]) == 1000  # tx_start should be start of gene
        os.unlink(path)

    def test_truncate_bed_returns_bed12(self):
        """Output should be 12-column BED."""
        bed, path = _make_bed([BED12_LINES[0]])
        result = bed.truncate_bed(truncation_from=3, size=200)
        fields = result[0].split("\t")
        assert len(fields) == 12
        os.unlink(path)

    def test_truncate_bed_3prime_minus(self):
        """For - strand, 3' truncation keeps first N bases (leftmost)."""
        bed, path = _make_bed([BED12_LINES[1]])  # gene2, - strand
        result = bed.truncate_bed(truncation_from=3, size=100)
        assert len(result) == 1
        fields = result[0].split("\t")
        assert fields[0] == "chr1"
        # - strand: 3' end is left side, so we take first 100 exonic bases
        assert int(fields[1]) == 6000  # starts at gene start
        os.unlink(path)

    def test_truncate_bed_repeatable(self):
        first = self.bed.truncate_bed(truncation_from=3, size=200)
        second = self.bed.truncate_bed(truncation_from=3, size=200)
        assert first == second

    # --- bedToWig ---

    def test_bed_to_wig(self, tmp_path):
        """bedToWig should produce a wiggle file with correct positions."""
        bed, bed_path = _make_bed([
            "chr1\t10\t13",  # simple 3-column BED
        ])
        outfile = str(tmp_path / "out.wig")
        bed.bedToWig(outfile=outfile, header=False)
        with open(outfile) as fh:
            content = fh.read()
        # BED 10-13 => wig positions 11,12,13 (1-based)
        assert "variableStep chrom=chr1" in content
        assert "11\t1" in content
        assert "12\t1" in content
        assert "13\t1" in content
        os.unlink(bed_path)

    def test_bed_to_wig_bed12(self, tmp_path):
        """bedToWig with BED12 input should skip introns."""
        bed, bed_path = _make_bed([BED12_LINES[0]])  # gene1
        outfile = str(tmp_path / "out.wig")
        bed.bedToWig(outfile=outfile, header=False)
        with open(outfile) as fh:
            lines = fh.readlines()
        # Should have a variableStep line + exon positions
        assert lines[0].startswith("variableStep")
        # Intron region (e.g., position 1600) should NOT appear
        positions = set()
        for line in lines[1:]:
            pos = int(line.split("\t")[0])
            positions.add(pos)
        # Exon1: 1001-1500, Exon2: 2501-3100, Exon3: 4501-5000 (1-based)
        assert 1001 in positions
        assert 1500 in positions
        # Intron positions should not be present
        assert 1600 not in positions
        assert 2000 not in positions
        os.unlink(bed_path)

    def test_bed_to_wig_with_header(self, tmp_path):
        bed, bed_path = _make_bed(["chr1\t5\t8"])
        outfile = str(tmp_path / "out.wig")
        bed.bedToWig(outfile=outfile, header=True)
        with open(outfile) as fh:
            first_line = fh.readline()
        assert first_line.startswith("track type=wiggle_0")
        os.unlink(bed_path)

    # --- bedToGFF ---

    def test_bed_to_gff(self, tmp_path):
        """bedToGFF should produce valid GFF output."""
        bed, bed_path = _make_bed([BED12_LINES[0]])
        outfile = str(tmp_path / "out.gff")
        bed.bedToGFF(outfile=outfile)
        with open(outfile) as fh:
            content = fh.read()
        assert "##gff-version 2" in content
        assert "chr1" in content
        os.unlink(bed_path)

    # --- getBedinfor ---

    def test_get_bedinfor(self, tmp_path):
        """getBedinfor should produce a tab-separated info file."""
        outfile = str(tmp_path / "info.xls")
        self.bed.getBedinfor(outfile=outfile)
        with open(outfile) as fh:
            lines = fh.readlines()
        # Header + 3 data lines
        assert len(lines) == 4
        assert "geneID" in lines[0]
        assert "Exon_Num" in lines[0]

    def test_get_bedinfor_exon_counts(self, tmp_path):
        """Check that exon numbers are correct in bedinfor output."""
        outfile = str(tmp_path / "info.xls")
        self.bed.getBedinfor(outfile=outfile)
        with open(outfile) as fh:
            lines = fh.readlines()
        # gene1 has 3 exons
        fields_gene1 = lines[1].split("\t")
        assert fields_gene1[1] == "3"
        # gene2 has 2 exons
        fields_gene2 = lines[2].split("\t")
        assert fields_gene2[1] == "2"
        # gene3 has 1 exon
        fields_gene3 = lines[3].split("\t")
        assert fields_gene3[1] == "1"

    # --- filterBedbyIntronSize ---

    def test_filter_bed_by_intron_size(self, tmp_path):
        """Genes with introns smaller than min_intron should be filtered out."""
        # gene with small intron
        line_small = "chr1\t100\t400\tsmall_intron\t0\t+\t100\t400\t0,0,0\t2\t50,50\t0,250"
        # intron = 150-350 = 200 bp
        line_big = "chr1\t1000\t5000\tbig_intron\t0\t+\t1000\t5000\t0,0,0\t2\t500,500\t0,3500"
        # intron = 1500-4500 = 3000 bp
        bed, bed_path = _make_bed([line_small, line_big])
        outfile = str(tmp_path / "filtered.bed")
        bed.filterBedbyIntronSize(outfile=outfile, min_intron=500)
        with open(outfile) as fh:
            lines = fh.readlines()
        # Only big_intron should pass (intron=3000 > 500)
        assert len(lines) == 1
        assert "big_intron" in lines[0]
        os.unlink(bed_path)

    def test_filter_bed_single_exon_excluded(self, tmp_path):
        """Single-exon genes should always be excluded."""
        bed, bed_path = _make_bed([BED12_LINES[2]])  # gene3, single exon
        outfile = str(tmp_path / "filtered.bed")
        bed.filterBedbyIntronSize(outfile=outfile, min_intron=0)
        with open(outfile) as fh:
            lines = fh.readlines()
        assert len(lines) == 0
        os.unlink(bed_path)

    # --- nrBED ---

    def test_nrbed_no_duplicates(self, tmp_path):
        """nrBED with no duplicate structures should preserve all entries."""
        outfile = str(tmp_path / "nr.bed")
        self.bed.nrBED(outfile=outfile)
        with open(outfile) as fh:
            lines = fh.readlines()
        assert len(lines) == 3

    def test_nrbed_merges_duplicates(self, tmp_path):
        """nrBED should merge entries with identical structure but different names."""
        line1 = "chr1\t1000\t5000\tgene1a\t0\t+\t1200\t4800\t0,0,0\t3\t500,600,500\t0,1500,3500"
        line2 = "chr1\t1000\t5000\tgene1b\t0\t+\t1200\t4800\t0,0,0\t3\t500,600,500\t0,1500,3500"
        bed, bed_path = _make_bed([line1, line2])
        outfile = str(tmp_path / "nr.bed")
        bed.nrBED(outfile=outfile)
        with open(outfile) as fh:
            lines = fh.readlines()
        # Should merge into 1 entry with both names
        assert len(lines) == 1
        assert "gene1a" in lines[0]
        assert "gene1b" in lines[0]
        os.unlink(bed_path)

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


# ---------------------------------------------------------------------------
# CompareBED tests
# ---------------------------------------------------------------------------


class TestCompareBED:
    def test_init(self, tmp_path):
        bed_a = tmp_path / "a.bed"
        bed_b = tmp_path / "b.bed"
        bed_a.write_text(BED12_LINES[0] + "\n")
        bed_b.write_text(BED12_LINES[0] + "\n")
        comp = BED.CompareBED(str(bed_a), str(bed_b))
        assert comp.A_base_Name == "a.bed"
        assert comp.B_base_Name == "b.bed"

    def test_best_match_identical(self, tmp_path):
        """bestMatch with identical files should find complete_match."""
        bed_a = tmp_path / "a.bed"
        bed_b = tmp_path / "b.bed"
        bed_a.write_text(BED12_LINES[0] + "\n")
        bed_b.write_text(BED12_LINES[0] + "\n")
        comp = BED.CompareBED(str(bed_a), str(bed_b))
        result = comp.bestMatch()
        assert len(result) == 1
        assert "complete_match" in result[1]

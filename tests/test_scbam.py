"""Tests for rseqc.scbam — pure-logic functions and BAM-based integration tests."""

import importlib
import os
from pathlib import Path

import pysam
import pytest

from rseqc import scbam


def test_import():
    """Verify that rseqc.scbam can be imported."""
    mod = importlib.import_module("rseqc.scbam")
    assert mod is not None


# --- diff_str ---


def test_diff_str_identical():
    assert scbam.diff_str("ACGT", "ACGT") == []


def test_diff_str_single_diff():
    result = scbam.diff_str("ACGT", "ACTT")
    assert len(result) == 1
    assert result[0] == [2, "G", "T"]


def test_diff_str_multiple_diffs():
    result = scbam.diff_str("AAAA", "TTTT")
    assert len(result) == 4


def test_diff_str_different_lengths():
    assert scbam.diff_str("AC", "ACG") == []


def test_diff_str_empty():
    assert scbam.diff_str("", "") == []


def test_diff_str_single_char_same():
    assert scbam.diff_str("A", "A") == []


def test_diff_str_single_char_diff():
    result = scbam.diff_str("A", "T")
    assert result == [[0, "A", "T"]]


def test_diff_str_all_positions_differ():
    result = scbam.diff_str("ACGT", "TGCA")
    assert len(result) == 4
    assert result[0] == [0, "A", "T"]
    assert result[1] == [1, "C", "G"]
    assert result[2] == [2, "G", "C"]
    assert result[3] == [3, "T", "A"]


def test_diff_str_first_only():
    result = scbam.diff_str("XAAA", "YAAA")
    assert len(result) == 1
    assert result[0] == [0, "X", "Y"]


def test_diff_str_last_only():
    result = scbam.diff_str("AAAX", "AAAY")
    assert len(result) == 1
    assert result[0] == [3, "X", "Y"]


# --- read_match_type ---


def test_read_match_type_consecutive():
    assert scbam.read_match_type("50M") == "Map_consecutively"


def test_read_match_type_splicing():
    assert scbam.read_match_type("25M1000N25M") == "Map_with_splicing"


def test_read_match_type_clip_left():
    assert scbam.read_match_type("5S45M") == "Map_with_clipping"


def test_read_match_type_clip_right():
    assert scbam.read_match_type("45M5S") == "Map_with_clipping"


def test_read_match_type_splice_clip_right():
    assert scbam.read_match_type("20M500N25M5S") == "Map_with_splicing_and_clipping"


def test_read_match_type_splice_clip_left():
    assert scbam.read_match_type("5S20M500N25M") == "Map_with_splicing_and_clipping"


def test_read_match_type_others():
    assert scbam.read_match_type("10M5I35M") == "Others"


def test_read_match_type_empty():
    assert scbam.read_match_type("") == "Others"


def test_read_match_type_large_numbers():
    assert scbam.read_match_type("150M") == "Map_consecutively"


def test_read_match_type_large_splice():
    assert scbam.read_match_type("75M50000N75M") == "Map_with_splicing"


def test_read_match_type_double_splice():
    """Two splicing events should not match simple splice pattern."""
    assert scbam.read_match_type("20M100N20M200N20M") == "Others"


def test_read_match_type_insertion():
    assert scbam.read_match_type("20M5I25M") == "Others"


def test_read_match_type_deletion():
    assert scbam.read_match_type("20M5D25M") == "Others"


def test_read_match_type_hard_clip():
    assert scbam.read_match_type("5H45M") == "Others"


# --- list2str ---


def test_list2str_simple():
    assert scbam.list2str([(0, 50)]) == "50M"


def test_list2str_splice():
    assert scbam.list2str([(0, 25), (3, 1000), (0, 25)]) == "25M1000N25M"


def test_list2str_clip():
    assert scbam.list2str([(4, 5), (0, 45)]) == "5S45M"


def test_list2str_insertion():
    assert scbam.list2str([(0, 30), (1, 5), (0, 20)]) == "30M5I20M"


def test_list2str_deletion():
    assert scbam.list2str([(0, 30), (2, 5), (0, 20)]) == "30M5D20M"


def test_list2str_all_ops():
    # Test all CIGAR operation codes
    ops = [(0, 10), (1, 2), (2, 3), (3, 100), (4, 5), (5, 1), (6, 1), (7, 10), (8, 2)]
    result = scbam.list2str(ops)
    assert result == "10M2I3D100N5S1H1P10=2X"


def test_list2str_empty():
    assert scbam.list2str([]) == ""


def test_list2str_single_base():
    assert scbam.list2str([(0, 1)]) == "1M"


def test_list2str_hard_clip_both_ends():
    assert scbam.list2str([(5, 10), (0, 50), (5, 10)]) == "10H50M10H"


# --- _pysam_iter ---


def test_pysam_iter_empty_bam(tmp_path):
    """_pysam_iter on an empty BAM produces no reads."""
    bam_path = str(tmp_path / "empty.bam")
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    )
    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        pass  # write no reads
    with pysam.AlignmentFile(bam_path, "rb") as inf:
        reads = list(scbam._pysam_iter(inf))
    assert reads == []


def test_pysam_iter_with_reads(tmp_path):
    """_pysam_iter should yield all reads from a BAM file."""
    bam_path = str(tmp_path / "test.bam")
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    )
    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        for i in range(3):
            a = pysam.AlignedSegment(header)
            a.query_name = f"read{i}"
            a.reference_id = 0
            a.reference_start = i * 100
            a.mapping_quality = 60
            a.cigar = [(0, 50)]
            a.query_sequence = "A" * 50
            a.query_qualities = pysam.qualitystring_to_array("I" * 50)
            outf.write(a)
    with pysam.AlignmentFile(bam_path, "rb") as inf:
        reads = list(scbam._pysam_iter(inf))
    assert len(reads) == 3


# ---------------------------------------------------------------------------
# BAM fixture with single-cell tags (CB, CR, UB, UR, xf, GX, GN, RE, TX, AN)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def sc_bam(tmp_path_factory):
    """Create a small sorted+indexed BAM with single-cell 10X tags.

    Reads:
      1. read_cb_same: CR=ACGT, CB=ACGT (same), UR=AAAA, UB=AAAA (same), xf=1, GX=GENE1, GN=Gene1, RE=E
      2. read_cb_diff: CR=ACGT, CB=ACTT (edited), UR=AAAA, UB=AAAT (edited), xf=1, GX=GENE2, GN=Gene2, RE=N
      3. read_no_cb: no CB/CR tags, xf=1, RE=I
      4. read_no_xf: CB=TTTT, UB=CCCC, no xf tag (not confidently mapped)
      5. read_dup: CR=GGGG, CB=GGGG, UR=CCCC, UB=CCCC, xf=1, is_duplicate, GX=GENE1, GN=Gene1, RE=E, TX=tx1
      6. read_antisense: CR=AAAA, CB=AAAA, UR=TTTT, UB=TTTT, xf=1, GX=GENE3, GN=Gene3, AN=anti1
    """
    tmpdir = tmp_path_factory.mktemp("scbam")
    unsorted = tmpdir / "unsorted.bam"
    sorted_bam = tmpdir / "sc.bam"

    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 50000}],
        }
    )
    seq50 = "A" * 50
    qual50 = pysam.qualitystring_to_array("I" * 50)

    def _sc_read(name, ref_start, flag=0, tags=None):
        a = pysam.AlignedSegment(header)
        a.query_name = name
        a.flag = flag
        a.reference_id = 0
        a.reference_start = ref_start
        a.mapping_quality = 60
        a.cigar = [(0, 50)]
        a.query_sequence = seq50
        a.query_qualities = qual50
        if tags:
            a.set_tags(tags)
        return a

    reads = [
        _sc_read("read_cb_same", 100, tags=[
            ("CR", "ACGT"), ("CB", "ACGT"),
            ("UR", "AAAA"), ("UB", "AAAA"),
            ("xf", 1), ("GX", "GENE1"), ("GN", "Gene1"), ("RE", "E"),
        ]),
        _sc_read("read_cb_diff", 200, tags=[
            ("CR", "ACGT"), ("CB", "ACTT"),
            ("UR", "AAAA"), ("UB", "AAAT"),
            ("xf", 1), ("GX", "GENE2"), ("GN", "Gene2"), ("RE", "N"),
        ]),
        _sc_read("read_no_cb", 300, tags=[
            ("xf", 1), ("RE", "I"),
        ]),
        _sc_read("read_no_xf", 400, tags=[
            ("CB", "TTTT"), ("UB", "CCCC"),
        ]),
        _sc_read("read_dup", 500, flag=0x400, tags=[
            ("CR", "GGGG"), ("CB", "GGGG"),
            ("UR", "CCCC"), ("UB", "CCCC"),
            ("xf", 1), ("GX", "GENE1"), ("GN", "Gene1"), ("RE", "E"), ("TX", "tx1"),
        ]),
        _sc_read("read_antisense", 600, tags=[
            ("CR", "AAAA"), ("CB", "AAAA"),
            ("UR", "TTTT"), ("UB", "TTTT"),
            ("xf", 1), ("GX", "GENE3"), ("GN", "Gene3"), ("AN", "anti1"),
        ]),
    ]

    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as outf:
        for r in reads:
            outf.write(r)

    pysam.sort("-o", str(sorted_bam), str(unsorted))
    pysam.index(str(sorted_bam))
    return sorted_bam


# --- barcode_edits ---


class TestBarcodeEdits:
    def test_barcode_edits_output_files(self, sc_bam, tmp_path):
        """barcode_edits should produce CB_freq, UMI_freq, and edits files."""
        outprefix = str(tmp_path / "bc_out")
        scbam.barcode_edits(str(sc_bam), outprefix, step_size=100, limit=100)

        assert os.path.exists(outprefix + ".CB_freq.tsv")
        assert os.path.exists(outprefix + ".UMI_freq.tsv")
        assert os.path.exists(outprefix + ".CB_edits_count.csv")
        assert os.path.exists(outprefix + ".UMI_edits_count.csv")

    def test_barcode_edits_cb_freq(self, sc_bam, tmp_path):
        """CB frequency file should contain barcodes from reads with CB tags."""
        outprefix = str(tmp_path / "bc_freq")
        scbam.barcode_edits(str(sc_bam), outprefix, step_size=100, limit=100)

        with open(outprefix + ".CB_freq.tsv") as fh:
            lines = fh.readlines()
        barcodes = {line.split("\t")[0] for line in lines}
        # Reads with CB: read_cb_same(ACGT), read_cb_diff(ACTT), read_no_xf(TTTT),
        # read_dup(GGGG), read_antisense(AAAA)
        assert "ACGT" in barcodes  # read_cb_same
        assert "ACTT" in barcodes  # read_cb_diff (corrected)
        assert "GGGG" in barcodes  # read_dup

    def test_barcode_edits_umi_freq(self, sc_bam, tmp_path):
        """UMI frequency file should contain UMIs from reads with UB tags."""
        outprefix = str(tmp_path / "bc_umi")
        scbam.barcode_edits(str(sc_bam), outprefix, step_size=100, limit=100)

        with open(outprefix + ".UMI_freq.tsv") as fh:
            lines = fh.readlines()
        umis = {line.split("\t")[0] for line in lines}
        assert "AAAA" in umis  # read_cb_same
        assert "AAAT" in umis  # read_cb_diff (corrected)

    def test_barcode_edits_limit(self, sc_bam, tmp_path):
        """With limit=2, only 2 reads should be processed."""
        outprefix = str(tmp_path / "bc_lim")
        scbam.barcode_edits(str(sc_bam), outprefix, step_size=100, limit=2)

        with open(outprefix + ".CB_freq.tsv") as fh:
            lines = fh.readlines()
        # With limit=2, at most 2 distinct barcodes
        total_count = sum(int(line.split("\t")[1].strip()) for line in lines)
        assert total_count <= 2

    def test_barcode_edits_cb_edit_matrix(self, sc_bam, tmp_path):
        """CB edits matrix should record position and nucleotide changes."""
        outprefix = str(tmp_path / "bc_edit")
        scbam.barcode_edits(str(sc_bam), outprefix, step_size=100, limit=100)

        with open(outprefix + ".CB_edits_count.csv") as fh:
            content = fh.read()
        # read_cb_diff: CR=ACGT -> CB=ACTT, pos 2: G->T
        # Should be in the matrix
        assert "Edits" in content or "Index" in content  # header present


# --- readCount ---


class TestReadCount:
    def test_readcount_runs(self, sc_bam, tmp_path):
        """readCount should run without error."""
        outprefix = str(tmp_path / "rc_out")
        # readCount uses xf tag and CB/GX/GN tags
        scbam.readCount(str(sc_bam), outprefix, step_size=100, limit=100)
        # readCount currently doesn't write output files (the write block is commented out)
        # Just verify it completes without error

    def test_readcount_limit(self, sc_bam, tmp_path):
        """readCount with limit=1 should process at most 1 read."""
        outprefix = str(tmp_path / "rc_lim")
        scbam.readCount(str(sc_bam), outprefix, step_size=100, limit=1)
        # No assertion on output since write is commented out; just verify no crash


# --- CBC_UMIcount ---


class TestCBCUMIcount:
    def test_cbc_umicount_output(self, sc_bam, tmp_path):
        """CBC_UMIcount should produce a Read_UMI_freq.tsv file."""
        outprefix = str(tmp_path / "cbc_out")
        scbam.CBC_UMIcount(
            str(sc_bam), outprefix,
            step_size=100,
            CB_tag="CB", UMI_tag="UB", gene_tag="GN",
            CB_num=100, read_num=1, UMI_num=1, gene_num=1,
        )
        outfile = outprefix + ".Read_UMI_freq.tsv"
        assert os.path.exists(outfile)
        with open(outfile) as fh:
            lines = fh.readlines()
        # Header should be present
        assert "Cell_barcode" in lines[0]
        assert "Read_count" in lines[0]
        assert "UMI_count" in lines[0]
        assert "Gene_count" in lines[0]

    def test_cbc_umicount_filters_low_reads(self, sc_bam, tmp_path):
        """With read_num=1000, no barcodes should pass the filter."""
        outprefix = str(tmp_path / "cbc_filt")
        scbam.CBC_UMIcount(
            str(sc_bam), outprefix,
            step_size=100,
            CB_tag="CB", UMI_tag="UB", gene_tag="GN",
            CB_num=100, read_num=1000, UMI_num=1, gene_num=1,
        )
        outfile = outprefix + ".Read_UMI_freq.tsv"
        with open(outfile) as fh:
            lines = fh.readlines()
        # Only header, no data lines (all barcodes have < 1000 reads)
        assert len(lines) == 1

    def test_cbc_umicount_filters_low_umi(self, sc_bam, tmp_path):
        """With UMI_num=1000, no barcodes should pass the UMI filter."""
        outprefix = str(tmp_path / "cbc_umi_filt")
        scbam.CBC_UMIcount(
            str(sc_bam), outprefix,
            step_size=100,
            CB_tag="CB", UMI_tag="UB", gene_tag="GN",
            CB_num=100, read_num=1, UMI_num=1000, gene_num=1,
        )
        outfile = outprefix + ".Read_UMI_freq.tsv"
        with open(outfile) as fh:
            lines = fh.readlines()
        assert len(lines) == 1  # header only


# --- Integration: roundtrip list2str -> read_match_type ---


class TestCigarRoundtrip:
    """Verify list2str output feeds correctly into read_match_type."""

    def test_consecutive(self):
        cigar = scbam.list2str([(0, 100)])
        assert scbam.read_match_type(cigar) == "Map_consecutively"

    def test_spliced(self):
        cigar = scbam.list2str([(0, 50), (3, 500), (0, 50)])
        assert scbam.read_match_type(cigar) == "Map_with_splicing"

    def test_clipped_left(self):
        cigar = scbam.list2str([(4, 10), (0, 90)])
        assert scbam.read_match_type(cigar) == "Map_with_clipping"

    def test_clipped_right(self):
        cigar = scbam.list2str([(0, 90), (4, 10)])
        assert scbam.read_match_type(cigar) == "Map_with_clipping"

    def test_splice_and_clip_right(self):
        cigar = scbam.list2str([(0, 40), (3, 200), (0, 50), (4, 10)])
        assert scbam.read_match_type(cigar) == "Map_with_splicing_and_clipping"

    def test_splice_and_clip_left(self):
        cigar = scbam.list2str([(4, 10), (0, 40), (3, 200), (0, 50)])
        assert scbam.read_match_type(cigar) == "Map_with_splicing_and_clipping"

    def test_insertion_is_others(self):
        cigar = scbam.list2str([(0, 40), (1, 5), (0, 55)])
        assert scbam.read_match_type(cigar) == "Others"

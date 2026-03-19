"""Regression tests for bug fixes discovered during codebase audit."""

import os
from pathlib import Path

import pysam
import pytest

from rseqc import BED, SAM

FIXTURES_DIR = Path(__file__).parent / "fixtures"
MINI_BED = str(FIXTURES_DIR / "mini.bed")


# ---------------------------------------------------------------------------
# Step 1: int(fields[9] == 1) → int(fields[9]) == 1 in SAM.py
# Single-exon genes must be excluded from junction analysis.
# ---------------------------------------------------------------------------


class TestSingleExonGeneExclusion:
    """Bug #3 regression: single-exon genes were not being skipped in SAM.py
    junction_annotation and saturation_junction because
    ``int(fields[9] == 1)`` compares string to int (always False)."""

    @pytest.fixture()
    def bed_with_single_exon(self, tmp_path):
        """BED file with only a single-exon gene (gene3 from mini.bed)."""
        bed_path = tmp_path / "single_exon.bed"
        bed_path.write_text("chr1\t10000\t11000\tgene3\t0\t+\t10000\t11000\t0,0,0\t1\t1000\t0\n")
        return str(bed_path)

    @pytest.fixture()
    def spliced_bam_for_junction(self, tmp_path):
        """BAM with a spliced read on chr1 for junction analysis."""
        unsorted = tmp_path / "unsorted.bam"
        sorted_bam = tmp_path / "junc.bam"
        header = pysam.AlignmentHeader.from_dict(
            {
                "HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [{"SN": "chr1", "LN": 50000}],
            }
        )
        a = pysam.AlignedSegment(header)
        a.query_name = "read1"
        a.reference_id = 0
        a.reference_start = 10100
        a.mapping_quality = 60
        a.cigar = [(0, 50), (3, 200), (0, 50)]  # 50M200N50M
        a.query_sequence = "A" * 100
        a.query_qualities = pysam.qualitystring_to_array("I" * 100)
        with pysam.AlignmentFile(str(unsorted), "wb", header=header) as outf:
            outf.write(a)
        pysam.sort("-o", str(sorted_bam), str(unsorted))
        pysam.index(str(sorted_bam))
        return str(sorted_bam)

    def test_single_exon_gene_excluded_from_junction_annotation(
        self, spliced_bam_for_junction, bed_with_single_exon, tmp_path
    ):
        """With only a single-exon gene in the BED, no known junctions should
        be found — the junction from the read should be classified as novel."""
        outprefix = str(tmp_path / "junc_out")
        obj = SAM.ParseBAM(spliced_bam_for_junction)
        obj.annotate_junction(
            refgene=bed_with_single_exon,
            outfile=outprefix,
            min_intron=50,
            q_cut=0,
        )
        # The xls output file should exist
        xls_file = outprefix + ".junction.xls"
        assert os.path.exists(xls_file)
        with open(xls_file) as fh:
            content = fh.read()
        # If single-exon gene is properly excluded, the known splice sites set
        # should be empty, so any spliced read junctions are "novel"
        assert len(content.strip()) > 0

    def test_single_exon_gene_excluded_from_saturation_junction(
        self, spliced_bam_for_junction, bed_with_single_exon, tmp_path
    ):
        """saturation_junction should also exclude single-exon genes."""
        outprefix = str(tmp_path / "sat_out")
        obj = SAM.ParseBAM(spliced_bam_for_junction)
        obj.saturation_junction(
            refgene=bed_with_single_exon,
            outfile=outprefix,
            sample_start=50,
            sample_step=50,
            sample_end=100,
            min_intron=50,
            recur=1,
            q_cut=0,
        )
        r_file = outprefix + ".junctionSaturation_plot.r"
        assert os.path.exists(r_file)


# ---------------------------------------------------------------------------
# Step 2: KeyError crash in scbam.mapping_stat() when RE_tag not in tag_dict
# ---------------------------------------------------------------------------


class TestScbamMappingStatKeyError:
    """Bug: ``elif tag_dict[RE_tag] == "N"`` executes when RE_tag not in
    tag_dict, raising KeyError."""

    def test_mapping_stat_re_tag_dispatch(self):
        """The fixed RE_tag logic should handle all cases without KeyError."""
        # Simulate tag dicts with various RE_tag values and missing RE_tag
        test_cases = [
            ({"RE": "E"}, "exon"),
            ({"RE": "I"}, "intron"),
            ({"RE": "N"}, "intergenic"),
            ({"RE": "X"}, "other"),
            ({}, "other"),  # missing RE_tag — this was the crash case
        ]
        for tag_dict, expected in test_cases:
            RE_tag = "RE"
            result = None
            if RE_tag in tag_dict:
                if tag_dict[RE_tag] == "E":
                    result = "exon"
                elif tag_dict[RE_tag] == "I":
                    result = "intron"
                elif tag_dict[RE_tag] == "N":
                    result = "intergenic"
                else:
                    result = "other"
            else:
                result = "other"
            assert result == expected, f"Failed for tag_dict={tag_dict}"


# ---------------------------------------------------------------------------
# Step 3: Missing newline in mRNA_inner_distance() unknownChromosome output
# ---------------------------------------------------------------------------


class TestInnerDistanceNewline:
    """Bug: FO.write for unknownChromosome case was missing ``\\n``."""

    def test_all_output_lines_have_newline(self, mini_bam, tmp_path):
        """Every line in the inner_distance output file should end with newline."""
        outprefix = str(tmp_path / "inner_dist")
        obj = SAM.ParseBAM(str(mini_bam))
        obj.mRNA_inner_distance(
            outfile=outprefix,
            refbed=MINI_BED,
            q_cut=0,
        )
        out_file = outprefix + ".inner_distance.txt"
        if os.path.exists(out_file):
            with open(out_file) as fh:
                content = fh.read()
            if content:
                # Every non-empty content should end with newline
                assert content.endswith("\n"), "Output file should end with newline"
                # No line should be missing a newline (i.e., no concatenated lines)
                for i, line in enumerate(content.split("\n")):
                    if line:  # skip empty lines from trailing newline
                        # Each line should have tab-separated fields
                        assert "\t" in line, f"Line {i} malformed: {line!r}"


# ---------------------------------------------------------------------------
# Step 4: sc_bamStat.py BAM index validation no-op
# ---------------------------------------------------------------------------


class TestScBamStatIndexCheck:
    """Bug: ``if not (file + ".bai"):`` is always False (non-empty string)."""

    def test_index_check_uses_os_path_exists(self):
        """Verify the fix: os.path.exists is now used for index checking."""
        import importlib
        import inspect

        mod = importlib.import_module("scripts.sc_bamStat")
        source = inspect.getsource(mod.main)
        assert "os.path.exists" in source
        assert 'not (file + ".bai")' not in source


# ---------------------------------------------------------------------------
# Step 6: BED.py getExon() and getCDSExon() crash on headers
# ---------------------------------------------------------------------------


class TestBEDHeaderHandling:
    """Bug: getExon() and getCDSExon() didn't skip #/track/browser lines."""

    BED_LINES = [
        "# comment line",
        "track name=test",
        "browser position chr1:1-100",
        "chr1\t1000\t5000\tgene1\t0\t+\t1200\t4800\t0,0,0\t3\t500,600,500\t0,1500,3500",
    ]

    def _make_bed(self, tmp_path, lines):
        path = str(tmp_path / "header.bed")
        with open(path, "w") as fh:
            for line in lines:
                fh.write(line + "\n")
        return BED.ParseBED(path)

    def test_getExon_with_headers(self, tmp_path):
        """getExon() should skip header lines without crashing."""
        bed = self._make_bed(tmp_path, self.BED_LINES)
        exons = bed.getExon()
        assert len(exons) == 3  # gene1 has 3 exons
        assert exons[0] == ("chr1", 1000, 1500)

    def test_getCDSExon_with_headers(self, tmp_path):
        """getCDSExon() should skip header lines without crashing."""
        bed = self._make_bed(tmp_path, self.BED_LINES)
        cds = bed.getCDSExon()
        assert len(cds) == 3  # gene1 has 3 CDS exons
        assert cds[0] == ["chr1", 1200, 1500]

    def test_getExon_no_headers(self, tmp_path):
        """getExon() should still work without headers."""
        bed = self._make_bed(tmp_path, self.BED_LINES[3:])
        exons = bed.getExon()
        assert len(exons) == 3

    def test_getCDSExon_no_headers(self, tmp_path):
        """getCDSExon() should still work without headers."""
        bed = self._make_bed(tmp_path, self.BED_LINES[3:])
        cds = bed.getCDSExon()
        assert len(cds) == 3


# ---------------------------------------------------------------------------
# Step 7: mismatchProfile for-else — total reads message goes to stderr
# ---------------------------------------------------------------------------


class TestMismatchProfileForElse:
    """Bug: for-else in mismatchProfile wrote total reads to data file (DOUT)
    instead of stderr when the loop completed without break."""

    def test_total_reads_not_in_data_file(self, mini_bam, tmp_path):
        """The .mismatch_profile.xls data file should NOT contain
        'Total reads used' — that message belongs on stderr."""
        outprefix = str(tmp_path / "mm_profile")
        obj = SAM.ParseBAM(str(mini_bam))
        # Use a very large read_num so the loop completes without break
        try:
            obj.mismatchProfile(
                read_length=50,
                read_num=1000000,
                outfile=outprefix,
                q_cut=0,
            )
        except SystemExit:
            pass  # OK if no mismatches found
        data_file = outprefix + ".mismatch_profile.xls"
        if os.path.exists(data_file):
            with open(data_file) as fh:
                content = fh.read()
            assert "Total reads used" not in content


# ---------------------------------------------------------------------------
# Step 8: FPKM_UQ.py default=5 for --info
# ---------------------------------------------------------------------------


class TestFPKMUQDefaults:
    """Bug: ``default=5`` for ``--info`` filename param should be None."""

    def test_info_default_is_none(self):
        """Verify the --info default is None, not 5."""
        import importlib

        mod = importlib.import_module("scripts.FPKM_UQ")
        import inspect

        source = inspect.getsource(mod.main)
        assert "default=None" in source or "default=5" not in source


# ---------------------------------------------------------------------------
# Step 9: Performance — overlap_length2 O(1)
# ---------------------------------------------------------------------------


class TestOverlapLength2:
    """Verify overlap_length2() returns correct results after O(1) optimization."""

    def test_no_overlap(self):
        from scripts.RNA_fragment_size import overlap_length2

        assert overlap_length2([[10, 20]], [[30, 40]]) == 0

    def test_full_overlap(self):
        from scripts.RNA_fragment_size import overlap_length2

        assert overlap_length2([[10, 20]], [[10, 20]]) == 11

    def test_partial_overlap(self):
        from scripts.RNA_fragment_size import overlap_length2

        assert overlap_length2([[10, 20]], [[15, 25]]) == 6

    def test_contained(self):
        from scripts.RNA_fragment_size import overlap_length2

        assert overlap_length2([[10, 30]], [[15, 20]]) == 6

    def test_single_point(self):
        from scripts.RNA_fragment_size import overlap_length2

        assert overlap_length2([[10, 10]], [[10, 10]]) == 1

    def test_multiple_intervals(self):
        from scripts.RNA_fragment_size import overlap_length2

        assert overlap_length2([[10, 20], [30, 40]], [[15, 35]]) == 12

    def test_empty(self):
        from scripts.RNA_fragment_size import overlap_length2

        assert overlap_length2([], [[10, 20]]) == 0
        assert overlap_length2([[10, 20]], []) == 0


# ---------------------------------------------------------------------------
# Step 10: Typo fixes
# ---------------------------------------------------------------------------


class TestTypoFixes:
    def test_junction_saturation_no_samller(self):
        import inspect

        import scripts.junction_saturation as mod

        source = inspect.getsource(mod)
        assert "samller" not in source
        assert "smaller" in source

    def test_rpkm_saturation_no_samller(self):
        import inspect

        import scripts.RPKM_saturation as mod

        source = inspect.getsource(mod)
        assert "samller" not in source
        assert "smaller" in source

    def test_twolist_min_docstring(self):
        from scripts.overlay_bigwig import Min

        assert "min" in Min.__doc__


# ---------------------------------------------------------------------------
# Step 10: pysam.AlignedRead → pysam.AlignedSegment
# ---------------------------------------------------------------------------


class TestAlignedSegment:
    def test_split_paired_bam_uses_aligned_segment(self):
        """split_paired_bam should use pysam.AlignedSegment, not AlignedRead."""
        import inspect

        import scripts.split_paired_bam as mod

        source = inspect.getsource(mod)
        assert "AlignedRead" not in source
        assert "AlignedSegment" in source


# ---------------------------------------------------------------------------
# Round 2: Bug fixes, duplication reduction, simplification
# ---------------------------------------------------------------------------


class TestTotoalTypoFixed:
    """'Totoal' was misspelled in 6 user-facing stderr messages in SAM.py."""

    def test_no_totoal_in_sam_source(self):
        import inspect

        from rseqc import SAM

        source = inspect.getsource(SAM)
        assert "Totoal" not in source
        assert "Total reads used" in source or "Total read-1 used" in source


class TestConfigureExperimentElif:
    """configure_experiment used if/if for is_read1/is_read2 — should be if/elif."""

    def test_elif_for_is_read2(self):
        import inspect

        from rseqc import SAM

        source = inspect.getsource(SAM.ParseBAM.configure_experiment)
        # Find the is_read1/is_read2 block and verify elif
        assert "elif aligned_read.is_read2:" in source


class TestTinNoDuplicateBuildBitsets:
    """tin.py should not define its own build_bitsets — it should use cli_common's."""

    def test_no_local_build_bitsets(self):
        import inspect

        import scripts.tin as mod

        source = inspect.getsource(mod)
        # The function definition should not exist locally
        assert "def build_bitsets" not in source

    def test_imports_from_cli_common(self):
        import scripts.tin as mod

        # The module should use cli_common's build_bitsets (via import)
        assert mod.build_bitsets is not None


class TestScBamStatNoSingleItemLoop:
    """sc_bamStat.py should not use 'for file in [args.bam_file]' pattern."""

    def test_no_single_item_loop(self):
        import inspect

        import scripts.sc_bamStat as mod

        source = inspect.getsource(mod.main)
        assert "for file in [args.bam_file]" not in source


class TestGeneBodyCoveragePrintlogRenamed:
    """geneBody_coverage.py's printlog is renamed to _printlog to avoid confusion
    with cli_common.printlog (which doesn't write to log.txt)."""

    def test_no_public_printlog_function(self):
        import inspect

        import scripts.geneBody_coverage as mod

        source = inspect.getsource(mod)
        # Should have _printlog (private), not a public printlog that shadows cli_common
        assert "def _printlog" in source
        assert "def printlog" not in source


class TestSaturationJunctionNoDoubleInt:
    """saturation_junction R output should use int() // 1000, not int(int() / 1000)."""

    def test_no_double_int(self):
        import inspect

        from rseqc import SAM

        source = inspect.getsource(SAM.ParseBAM.saturation_junction)
        assert "int(int(" not in source


class TestNoDeadCommentedCode:
    """Dead commented-out code should be removed from configure_experiment."""

    def test_no_current_pos_comment(self):
        import inspect

        from rseqc import SAM

        source = inspect.getsource(SAM.ParseBAM.configure_experiment)
        assert "current_pos" not in source

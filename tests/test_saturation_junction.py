"""Tests for SAM.ParseBAM.saturation_junction() incremental counting
and configure_experiment() fast-path optimization."""

import collections
import io
from pathlib import Path
from unittest.mock import patch

import pysam
import pytest

from rseqc import SAM

FIXTURES_DIR = Path(__file__).parent / "fixtures"
MINI_BED = str(FIXTURES_DIR / "mini.bed")


# ---------------------------------------------------------------------------
# Fixture: BAM with multiple splice junctions for saturation testing
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def spliced_bam(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Create a BAM with many spliced reads hitting known and novel junctions.

    mini.bed gene1 introns: 1500-2500, 3100-4500
    mini.bed gene2 introns: 7000-8000
    We create reads with known splice junctions (matching these introns)
    and novel junctions (not in the BED file).
    """
    tmpdir = tmp_path_factory.mktemp("spliced_bam")
    unsorted_bam = tmpdir / "unsorted.bam"
    sorted_bam = tmpdir / "spliced.bam"

    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [
                {"SN": "chr1", "LN": 50000},
                {"SN": "chr2", "LN": 50000},
            ],
        }
    )

    def _make_spliced_read(
        name: str,
        ref_id: int,
        ref_start: int,
        intron_start: int,
        intron_len: int,
        mapq: int = 60,
    ) -> pysam.AlignedSegment:
        """Create a read with one splice junction: exon1 + N + exon2."""
        exon1_len = intron_start - ref_start
        exon2_len = 30
        cigar = [(0, exon1_len), (3, intron_len), (0, exon2_len)]
        seq_len = exon1_len + exon2_len

        a = pysam.AlignedSegment(header)
        a.query_name = name
        a.flag = 0
        a.reference_id = ref_id
        a.reference_start = ref_start
        a.mapping_quality = mapq
        a.cigar = cigar
        a.query_sequence = "A" * seq_len
        a.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
        return a

    reads = []
    # Known junctions (match mini.bed introns):
    # gene1 intron1: 1500-2500 (exon_ends[0]=1500, exon_starts[1]=2500)
    for i in range(10):
        reads.append(_make_spliced_read(f"known1_{i}", 0, 1450, 1500, 1000))
    # gene1 intron2: 3100-4500
    for i in range(5):
        reads.append(_make_spliced_read(f"known2_{i}", 0, 3050, 3100, 1400))
    # gene2 intron: 7000-8000
    for i in range(8):
        reads.append(_make_spliced_read(f"known3_{i}", 0, 6950, 7000, 1000))
    # Novel junctions (not in mini.bed):
    for i in range(6):
        reads.append(_make_spliced_read(f"novel1_{i}", 0, 12000, 12050, 500))
    for i in range(4):
        reads.append(_make_spliced_read(f"novel2_{i}", 0, 20000, 20050, 800))
    # Novel junction on chr2
    for i in range(3):
        reads.append(_make_spliced_read(f"novel3_{i}", 1, 1000, 1050, 600))

    with pysam.AlignmentFile(str(unsorted_bam), "wb", header=header) as outf:
        for read in reads:
            outf.write(read)

    pysam.sort("-o", str(sorted_bam), str(unsorted_bam))
    pysam.index(str(sorted_bam))

    return sorted_bam


# ---------------------------------------------------------------------------
# saturation_junction tests
# ---------------------------------------------------------------------------


class TestSaturationJunction:
    def test_produces_r_script_file(self, spliced_bam, tmp_path):
        """saturation_junction creates the .junctionSaturation_plot.r file."""
        outprefix = str(tmp_path / "sat_test")
        obj = SAM.ParseBAM(str(spliced_bam))
        with patch("sys.stderr", new_callable=io.StringIO):
            obj.saturation_junction(
                refgene=MINI_BED,
                outfile=outprefix,
                sample_start=50,
                sample_step=50,
                sample_end=100,
                min_intron=50,
                recur=1,
                q_cut=0,
            )
        r_file = Path(outprefix + ".junctionSaturation_plot.r")
        assert r_file.exists()
        assert r_file.stat().st_size > 0

    def test_r_script_contains_data_vectors(self, spliced_bam, tmp_path):
        """R script contains x, y (known), z (all), w (novel) vectors."""
        outprefix = str(tmp_path / "sat_vec")
        obj = SAM.ParseBAM(str(spliced_bam))
        with patch("sys.stderr", new_callable=io.StringIO):
            obj.saturation_junction(
                refgene=MINI_BED,
                outfile=outprefix,
                sample_start=50,
                sample_step=50,
                sample_end=100,
                min_intron=50,
                recur=1,
                q_cut=0,
            )
        content = Path(outprefix + ".junctionSaturation_plot.r").read_text()
        assert "x=c(" in content
        assert "y=c(" in content  # known
        assert "z=c(" in content  # all
        assert "w=c(" in content  # novel

    def test_known_plus_novel_equals_all(self, spliced_bam, tmp_path):
        """At each sampling point: known + novel == all junctions."""
        outprefix = str(tmp_path / "sat_sum")
        obj = SAM.ParseBAM(str(spliced_bam))
        with patch("sys.stderr", new_callable=io.StringIO):
            obj.saturation_junction(
                refgene=MINI_BED,
                outfile=outprefix,
                sample_start=25,
                sample_step=25,
                sample_end=100,
                min_intron=50,
                recur=1,
                q_cut=0,
            )
        content = Path(outprefix + ".junctionSaturation_plot.r").read_text()
        # Parse the vectors from R script
        known = _parse_r_vector(content, "y")
        all_j = _parse_r_vector(content, "z")
        novel = _parse_r_vector(content, "w")
        assert len(known) == len(all_j) == len(novel)
        for k, a, n in zip(known, all_j, novel):
            assert k + n == a, f"known({k}) + novel({n}) != all({a})"

    def test_counts_are_nondecreasing(self, spliced_bam, tmp_path):
        """known, all, and novel counts are non-decreasing across percentiles."""
        outprefix = str(tmp_path / "sat_mono")
        obj = SAM.ParseBAM(str(spliced_bam))
        with patch("sys.stderr", new_callable=io.StringIO):
            obj.saturation_junction(
                refgene=MINI_BED,
                outfile=outprefix,
                sample_start=25,
                sample_step=25,
                sample_end=100,
                min_intron=50,
                recur=1,
                q_cut=0,
            )
        content = Path(outprefix + ".junctionSaturation_plot.r").read_text()
        for label in ("y", "z", "w"):
            vals = _parse_r_vector(content, label)
            for i in range(1, len(vals)):
                assert vals[i] >= vals[i - 1], f"Vector {label} decreased: {vals[i - 1]} -> {vals[i]}"

    def test_recur_threshold_delays_known_count(self, spliced_bam, tmp_path):
        """With recur>1, a junction needs multiple observations to be 'known'."""
        # Run with recur=1
        out1 = str(tmp_path / "sat_r1")
        obj1 = SAM.ParseBAM(str(spliced_bam))
        with patch("sys.stderr", new_callable=io.StringIO):
            with patch("random.shuffle"):  # deterministic
                obj1.saturation_junction(
                    refgene=MINI_BED,
                    outfile=out1,
                    sample_start=50,
                    sample_step=50,
                    sample_end=100,
                    min_intron=50,
                    recur=1,
                    q_cut=0,
                )
        # Run with recur=3
        out2 = str(tmp_path / "sat_r3")
        obj2 = SAM.ParseBAM(str(spliced_bam))
        with patch("sys.stderr", new_callable=io.StringIO):
            with patch("random.shuffle"):  # deterministic
                obj2.saturation_junction(
                    refgene=MINI_BED,
                    outfile=out2,
                    sample_start=50,
                    sample_step=50,
                    sample_end=100,
                    min_intron=50,
                    recur=3,
                    q_cut=0,
                )
        content1 = Path(out1 + ".junctionSaturation_plot.r").read_text()
        content2 = Path(out2 + ".junctionSaturation_plot.r").read_text()
        known_r1 = _parse_r_vector(content1, "y")
        known_r3 = _parse_r_vector(content2, "y")
        # With higher recur, known count should be <= recur=1 at each step
        for k1, k3 in zip(known_r1, known_r3):
            assert k3 <= k1, f"recur=3 known({k3}) > recur=1 known({k1})"

    def test_novel_counts_independent_of_recur(self, spliced_bam, tmp_path):
        """Novel junction counts don't depend on recur threshold."""
        out1 = str(tmp_path / "sat_novel_r1")
        obj1 = SAM.ParseBAM(str(spliced_bam))
        with patch("sys.stderr", new_callable=io.StringIO):
            with patch("random.shuffle"):
                obj1.saturation_junction(
                    refgene=MINI_BED,
                    outfile=out1,
                    sample_start=50,
                    sample_step=50,
                    sample_end=100,
                    min_intron=50,
                    recur=1,
                    q_cut=0,
                )
        out2 = str(tmp_path / "sat_novel_r5")
        obj2 = SAM.ParseBAM(str(spliced_bam))
        with patch("sys.stderr", new_callable=io.StringIO):
            with patch("random.shuffle"):
                obj2.saturation_junction(
                    refgene=MINI_BED,
                    outfile=out2,
                    sample_start=50,
                    sample_step=50,
                    sample_end=100,
                    min_intron=50,
                    recur=5,
                    q_cut=0,
                )
        content1 = Path(out1 + ".junctionSaturation_plot.r").read_text()
        content2 = Path(out2 + ".junctionSaturation_plot.r").read_text()
        novel_r1 = _parse_r_vector(content1, "w")
        novel_r5 = _parse_r_vector(content2, "w")
        assert novel_r1 == novel_r5

    def test_incremental_matches_reference(self, spliced_bam, tmp_path):
        """Verify incremental counting produces identical results to a naive rescan.

        This is the key correctness test: we run saturation_junction and then
        independently verify counts by doing a full rescan of the accumulated dict.
        """
        outprefix = str(tmp_path / "sat_ref")
        obj = SAM.ParseBAM(str(spliced_bam))

        # Capture the R script output
        with patch("sys.stderr", new_callable=io.StringIO):
            with patch("random.shuffle"):  # deterministic order
                obj.saturation_junction(
                    refgene=MINI_BED,
                    outfile=outprefix,
                    sample_start=25,
                    sample_step=25,
                    sample_end=100,
                    min_intron=50,
                    recur=2,
                    q_cut=0,
                )
        content = Path(outprefix + ".junctionSaturation_plot.r").read_text()
        actual_known = _parse_r_vector(content, "y")
        actual_all = _parse_r_vector(content, "z")
        actual_novel = _parse_r_vector(content, "w")

        # Now do the same computation with a naive reference implementation
        from rseqc import bam_cigar
        from rseqc.SAM import _passes_qc, _pysam_iter

        # Parse known splice sites from BED
        knownSpliceSites = set()
        chrom_list = set()
        with open(MINI_BED) as fh:
            for line in fh:
                if line.startswith(("#", "track", "browser")):
                    continue
                fields = line.split()
                if len(fields) < 12:
                    continue
                chrom = fields[0].upper()
                chrom_list.add(chrom)
                tx_start = int(fields[1])
                if int(fields[9] == 1):
                    continue
                exon_starts = list(map(int, fields[11].rstrip(",\n").split(",")))
                exon_starts = [x + tx_start for x in exon_starts]
                exon_ends = list(map(int, fields[10].rstrip(",\n").split(",")))
                exon_ends = [x + y for x, y in zip(exon_starts, exon_ends)]
                intron_start = exon_ends[:-1]
                intron_end = exon_starts[1:]
                for st, end in zip(intron_start, intron_end):
                    knownSpliceSites.add(chrom + ":" + str(st) + "-" + str(end))

        # Collect splice sites from BAM
        obj2 = SAM.ParseBAM(str(spliced_bam))
        samSpliceSites = []
        for aligned_read in _pysam_iter(obj2.samfile):
            if not _passes_qc(aligned_read, 0):
                continue
            try:
                chrom = obj2.samfile.getrname(aligned_read.tid).upper()
            except (ValueError, KeyError, AttributeError):
                continue
            if chrom not in chrom_list:
                continue
            hit_st = aligned_read.pos
            intron_blocks = bam_cigar.fetch_intron(hit_st, aligned_read.cigar)
            for intrn in intron_blocks:
                if intrn[1] - intrn[0] < 50:
                    continue
                samSpliceSites.append(chrom + ":" + str(intrn[0]) + "-" + str(intrn[1]))

        # Don't shuffle (matches patch("random.shuffle"))
        # Naive rescan at each percentile
        uniqSpliceSites: dict[str, int] = collections.defaultdict(int)
        ref_known = []
        ref_all = []
        ref_novel = []
        recur = 2
        tmp_steps = list(range(25, 100, 25))
        tmp_steps.append(100)
        SR_num = len(samSpliceSites)
        for pertl in tmp_steps:
            index_st = int(SR_num * ((pertl - 25) / 100.0))
            index_end = int(SR_num * (pertl / 100.0))
            if index_st < 0:
                index_st = 0
            for i in range(index_st, index_end):
                uniqSpliceSites[samSpliceSites[i]] += 1

            ref_all.append(len(uniqSpliceSites))
            known_n = sum(1 for sj in uniqSpliceSites if sj in knownSpliceSites and uniqSpliceSites[sj] >= recur)
            ref_known.append(known_n)
            novel_n = sum(1 for sj in uniqSpliceSites if sj not in knownSpliceSites)
            ref_novel.append(novel_n)

        assert actual_known == ref_known
        assert actual_all == ref_all
        assert actual_novel == ref_novel


# ---------------------------------------------------------------------------
# configure_experiment tests
# ---------------------------------------------------------------------------


class TestConfigureExperimentOptimized:
    def test_single_end_fractions_stable(self, tmp_path):
        """Single-end reads produce consistent spec1+spec2+other≈1."""
        bam_path = _make_single_end_bam(tmp_path)
        obj = SAM.ParseBAM(str(bam_path))
        with patch("sys.stderr", new_callable=io.StringIO):
            result = obj.configure_experiment(refbed=MINI_BED, sample_size=100, q_cut=0)
        assert result[0] == "SingleEnd"
        total = result[1] + result[2] + result[3]
        assert abs(total - 1.0) < 0.01

    def test_paired_end_fractions_stable(self, mini_bam):
        """Paired-end reads produce consistent spec1+spec2+other≈1."""
        obj = SAM.ParseBAM(str(mini_bam))
        with patch("sys.stderr", new_callable=io.StringIO):
            result = obj.configure_experiment(refbed=MINI_BED, sample_size=100, q_cut=0)
        if result[0] == "PairEnd":
            total = result[1] + result[2] + result[3]
            assert abs(total - 1.0) < 0.01

    def test_multi_gene_overlap_deterministic(self, tmp_path):
        """Reads overlapping genes on both strands produce deterministic results.

        The optimized code uses sorted(set(hits)) for multi-hit reads, which
        should be deterministic unlike the original set() iteration order.
        """
        bam_path = _make_single_end_bam(tmp_path)
        # Run twice and verify identical results
        results = []
        for _ in range(2):
            obj = SAM.ParseBAM(str(bam_path))
            with patch("sys.stderr", new_callable=io.StringIO):
                result = obj.configure_experiment(refbed=MINI_BED, sample_size=100, q_cut=0)
            results.append(result)
        assert results[0] == results[1]

    def test_returns_four_elements(self, mini_bam):
        """configure_experiment returns [protocol, spec1, spec2, other]."""
        obj = SAM.ParseBAM(str(mini_bam))
        with patch("sys.stderr", new_callable=io.StringIO):
            result = obj.configure_experiment(refbed=MINI_BED, sample_size=100, q_cut=0)
        assert isinstance(result, list)
        assert len(result) == 4
        assert result[0] in ("PairEnd", "SingleEnd", "Mixture")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _parse_r_vector(content: str, var_name: str) -> list[int]:
    """Extract integer values from an R vector assignment like 'y=c(1,2,3)'."""
    for line in content.splitlines():
        if line.startswith(var_name + "=c("):
            inner = line.split("=c(", 1)[1].rstrip(")")
            return [int(x) for x in inner.split(",")]
    raise ValueError(f"Vector {var_name} not found in R script")


def _make_single_end_bam(tmp_path: Path) -> Path:
    """Create a BAM with single-end reads only (no paired flags)."""
    unsorted = tmp_path / "se_unsorted.bam"
    sorted_bam = tmp_path / "se.bam"

    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 50000}],
        }
    )

    reads = []
    # Forward-strand reads in gene1 region
    for i in range(20):
        a = pysam.AlignedSegment(header)
        a.query_name = f"se_fwd_{i}"
        a.flag = 0  # forward, single-end
        a.reference_id = 0
        a.reference_start = 1050 + i * 10
        a.mapping_quality = 60
        a.cigar = [(0, 50)]
        a.query_sequence = "A" * 50
        a.query_qualities = pysam.qualitystring_to_array("I" * 50)
        reads.append(a)
    # Reverse-strand reads in gene2 region
    for i in range(20):
        a = pysam.AlignedSegment(header)
        a.query_name = f"se_rev_{i}"
        a.flag = 0x10  # reverse strand, single-end
        a.reference_id = 0
        a.reference_start = 6050 + i * 10
        a.mapping_quality = 60
        a.cigar = [(0, 50)]
        a.query_sequence = "A" * 50
        a.query_qualities = pysam.qualitystring_to_array("I" * 50)
        reads.append(a)

    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as outf:
        for read in reads:
            outf.write(read)

    pysam.sort("-o", str(sorted_bam), str(unsorted))
    pysam.index(str(sorted_bam))
    return sorted_bam

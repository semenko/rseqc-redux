"""Tests for scripts/split_bam.py — verify read routing before refactoring."""

import subprocess
import sys
from pathlib import Path

import pysam
import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture()
def split_bam_fixture(tmp_path: Path) -> tuple[Path, Path]:
    """Create a BAM and BED file exercising all split_bam code paths.

    BED defines exon regions on chr1: 1000-1500, 2500-3100 (gene1).
    No exon regions on chr2.

    Reads cover:
      - qcfail → junk
      - unmapped → junk
      - mapped, chrom not in exon_ranges → ex
      - mate_is_unmapped, position in exon → in
      - mate_is_unmapped, position not in exon → ex
      - both mapped, read_start in exon → in
      - both mapped, only mate_start in exon → in
      - both mapped, neither in exon → ex
    """
    bed_file = tmp_path / "test.bed"
    bed_file.write_text("chr1\t1000\t3100\tgene1\t0\t+\t1200\t3000\t0,0,0\t2\t500,600\t0,1500\n")
    # Exon regions: chr1:1000-1500, chr1:2500-3100

    bam_unsorted = tmp_path / "unsorted.bam"
    bam_sorted = tmp_path / "test.bam"

    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [
                {"SN": "chr1", "LN": 50000},
                {"SN": "chr2", "LN": 50000},
            ],
        }
    )

    seq50 = "A" * 50
    qual50 = pysam.qualitystring_to_array("I" * 50)

    def _make_read(
        name: str,
        ref_id: int = 0,
        ref_start: int = 0,
        mapq: int = 60,
        flag: int = 0,
        mate_ref_id: int = -1,
        mate_ref_start: int = 0,
    ) -> pysam.AlignedSegment:
        a = pysam.AlignedSegment(header)
        a.query_name = name
        a.flag = flag
        if flag & 0x4:  # unmapped
            a.reference_id = -1
            a.reference_start = 0
            a.mapping_quality = 0
            a.cigar = None
        else:
            a.reference_id = ref_id
            a.reference_start = ref_start
            a.mapping_quality = mapq
            a.cigar = [(0, 50)]  # 50M
        a.query_sequence = seq50
        a.query_qualities = qual50
        if mate_ref_id >= 0:
            a.next_reference_id = mate_ref_id
            a.next_reference_start = mate_ref_start
        return a

    reads = [
        # 1. qcfail → junk
        _make_read("read_qcfail", ref_start=1050, flag=0x200),
        # 2. unmapped → junk
        _make_read("read_unmapped", flag=0x4),
        # 3. chrom not in exon_ranges → ex
        _make_read("read_chr2", ref_id=1, ref_start=500),
        # 4. mate_is_unmapped, position in exon → in
        _make_read(
            "read_mate_unmapped_in",
            ref_start=1050,
            flag=0x1 | 0x8,
        ),
        # 5. mate_is_unmapped, position not in exon → ex
        _make_read(
            "read_mate_unmapped_out",
            ref_start=3200,
            flag=0x1 | 0x8,
        ),
        # 6. both mapped, read_start in exon → in
        _make_read(
            "read_pair_read_in",
            ref_start=1050,
            flag=0x1 | 0x2 | 0x40,
            mate_ref_id=0,
            mate_ref_start=3200,
        ),
        # 7. both mapped, only mate_start in exon → in
        _make_read(
            "read_pair_mate_in",
            ref_start=3200,
            flag=0x1 | 0x2 | 0x40,
            mate_ref_id=0,
            mate_ref_start=1050,
        ),
        # 8. both mapped, neither in exon → ex
        _make_read(
            "read_pair_neither",
            ref_start=3200,
            flag=0x1 | 0x2 | 0x40,
            mate_ref_id=0,
            mate_ref_start=3300,
        ),
    ]

    with pysam.AlignmentFile(str(bam_unsorted), "wb", header=header) as outf:
        for read in reads:
            outf.write(read)

    pysam.sort("-o", str(bam_sorted), str(bam_unsorted))
    pysam.index(str(bam_sorted))

    return bam_sorted, bed_file


def _count_reads(bam_path: str) -> int:
    """Count total reads in a BAM file."""
    count = 0
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as f:
        for _ in f:
            count += 1
    return count


def _read_names(bam_path: str) -> set[str]:
    """Get set of read names from a BAM file."""
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as f:
        return {read.query_name for read in f}


def test_split_bam_routing(
    split_bam_fixture: tuple[Path, Path],
    tmp_path: Path,
) -> None:
    """Verify each read is routed to the correct output file."""
    bam_path, bed_path = split_bam_fixture
    outprefix = str(tmp_path / "out")

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.split_bam",
            "-i",
            str(bam_path),
            "-r",
            str(bed_path),
            "-o",
            outprefix,
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, f"split_bam failed: {result.stderr}"

    in_bam = outprefix + ".in.bam"
    ex_bam = outprefix + ".ex.bam"
    junk_bam = outprefix + ".junk.bam"

    # Verify output files exist
    assert Path(in_bam).exists()
    assert Path(ex_bam).exists()
    assert Path(junk_bam).exists()

    # Verify read counts
    assert _count_reads(junk_bam) == 2  # qcfail + unmapped
    assert _count_reads(in_bam) == 3  # mate_unmapped_in + pair_read_in + pair_mate_in
    assert _count_reads(ex_bam) == 3  # chr2 + mate_unmapped_out + pair_neither

    # Verify specific reads in each file
    junk_names = _read_names(junk_bam)
    assert "read_qcfail" in junk_names
    assert "read_unmapped" in junk_names

    in_names = _read_names(in_bam)
    assert "read_mate_unmapped_in" in in_names
    assert "read_pair_read_in" in in_names
    assert "read_pair_mate_in" in in_names

    ex_names = _read_names(ex_bam)
    assert "read_chr2" in ex_names
    assert "read_mate_unmapped_out" in ex_names
    assert "read_pair_neither" in ex_names


def test_split_bam_stdout_counts(
    split_bam_fixture: tuple[Path, Path],
    tmp_path: Path,
) -> None:
    """Verify the summary counts printed to stdout."""
    bam_path, bed_path = split_bam_fixture
    outprefix = str(tmp_path / "out")

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.split_bam",
            "-i",
            str(bam_path),
            "-r",
            str(bed_path),
            "-o",
            outprefix,
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0

    # Total should be 8
    assert "Total records:" in result.stdout
    lines = result.stdout.strip().split("\n")
    # Parse the counts from the output — format is "%-55s%d", so extract
    # trailing digits from each line.
    import re

    counts = {}
    for line in lines:
        m = re.search(r"(\d+)\s*$", line)
        if not m:
            continue
        val = int(m.group(1))
        if "Total records:" in line:
            counts["total"] = val
        elif ".in.bam" in line:
            counts["in"] = val
        elif ".ex.bam" in line:
            counts["ex"] = val
        elif ".junk.bam" in line:
            counts["junk"] = val

    assert counts["total"] == 8
    assert counts["in"] == 3
    assert counts["ex"] == 3
    assert counts["junk"] == 2


def test_split_bam_with_mini_fixture(mini_bam: Path, tmp_path: Path) -> None:
    """Integration test with the shared mini_bam fixture."""
    bed_file = str(FIXTURES_DIR / "mini.bed")
    outprefix = str(tmp_path / "mini_split")

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.split_bam",
            "-i",
            str(mini_bam),
            "-r",
            bed_file,
            "-o",
            outprefix,
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, f"split_bam failed: {result.stderr}"

    in_bam = outprefix + ".in.bam"
    ex_bam = outprefix + ".ex.bam"
    junk_bam = outprefix + ".junk.bam"

    # mini_bam has 11 reads total:
    #   - 1 unmapped → junk
    #   - 1 on chr2 (no exons) → ex
    #   - 9 on chr1 in exon regions → in
    in_count = _count_reads(in_bam)
    ex_count = _count_reads(ex_bam)
    junk_count = _count_reads(junk_bam)

    assert junk_count == 1  # read_unmapped
    assert ex_count == 1  # read_chr2
    assert in_count == 9  # all others hit exon regions
    assert in_count + ex_count + junk_count == 11

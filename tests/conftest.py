from pathlib import Path

import pysam
import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def fixtures_dir() -> Path:
    """Return the path to the test fixtures directory."""
    return FIXTURES_DIR


@pytest.fixture(scope="session")
def mini_bam(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Create a small sorted+indexed BAM file matching mini.bed genes.

    Gene layout (from mini.bed):
      chr1 gene1 (+, 3 exons: 1000-1500, 2500-3100, 4500-5000)
      chr1 gene2 (-, 2 exons: 6000-7000, 8000-9000)
      chr1 gene3 (+, 1 exon: 10000-11000)

    Reads:
      1. read_unique1: 50M at chr1:1050, MAPQ=60, unique mapped (gene1 exon1)
      2. read_unique2: 50M at chr1:1100, MAPQ=60, unique mapped (gene1 exon1)
      3. read_spliced: 50M1400N50M at chr1:1450, MAPQ=60, spliced across gene1 intron1
      4. read_reverse: 50M at chr1:6500, MAPQ=60, reverse strand (gene2)
      5. read_chr2: 50M at chr2:500, MAPQ=60, different chromosome
      6. read_unmapped: unmapped (flag 0x4)
      7. read_duplicate: 50M at chr1:1050, MAPQ=60, PCR duplicate (flag 0x400)
      8. read_secondary: 50M at chr1:1050, MAPQ=60, secondary alignment (flag 0x100)
      9. read_lowmapq: 50M at chr1:1050, MAPQ=2, low mapping quality
     10. read_pair_r1: 50M at chr1:1050, MAPQ=60, proper pair R1 (flags 0x1|0x2|0x40)
         read_pair_r2: 50M at chr1:1200, MAPQ=60, proper pair R2 (flags 0x1|0x2|0x80)
    """
    tmpdir = tmp_path_factory.mktemp("bam")
    unsorted_bam = tmpdir / "unsorted.bam"
    sorted_bam = tmpdir / "mini.bam"

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
        cigar: list[tuple[int, int]] | None = None,
        mapq: int = 60,
        flag: int = 0,
        mate_ref_id: int = -1,
        mate_ref_start: int = 0,
        template_length: int = 0,
    ) -> pysam.AlignedSegment:
        a = pysam.AlignedSegment(header)
        a.query_name = name
        a.flag = flag
        if flag & 0x4:  # unmapped
            a.reference_id = -1
            a.reference_start = 0
            a.mapping_quality = 0
            a.cigar = None
            a.query_sequence = seq50
            a.query_qualities = qual50
        else:
            a.reference_id = ref_id
            a.reference_start = ref_start
            a.mapping_quality = mapq
            a.cigar = cigar or [(0, 50)]  # 50M
            a.query_sequence = "A" * sum(ln for op, ln in a.cigar if op in (0, 1, 4))
            a.query_qualities = pysam.qualitystring_to_array("I" * len(a.query_sequence))
        if mate_ref_id >= 0:
            a.next_reference_id = mate_ref_id
            a.next_reference_start = mate_ref_start
            a.template_length = template_length
        return a

    reads = [
        # 1. unique mapped gene1 exon1
        _make_read("read_unique1", ref_start=1050),
        # 2. unique mapped gene1 exon1
        _make_read("read_unique2", ref_start=1100),
        # 3. spliced: 50M 1400N 50M (spans gene1 intron1: 1500-2500 + partial exon2)
        _make_read(
            "read_spliced",
            ref_start=1450,
            cigar=[(0, 50), (3, 1400), (0, 50)],
        ),
        # 4. reverse strand, gene2
        _make_read("read_reverse", ref_start=6500, flag=0x10),
        # 5. different chromosome
        _make_read("read_chr2", ref_id=1, ref_start=500),
        # 6. unmapped
        _make_read("read_unmapped", flag=0x4),
        # 7. PCR duplicate
        _make_read("read_duplicate", ref_start=1050, flag=0x400),
        # 8. secondary alignment
        _make_read("read_secondary", ref_start=1050, flag=0x100),
        # 9. low MAPQ
        _make_read("read_lowmapq", ref_start=1050, mapq=2),
        # 10. proper pair R1
        _make_read(
            "read_pair",
            ref_start=1050,
            flag=0x1 | 0x2 | 0x40,
            mate_ref_id=0,
            mate_ref_start=1200,
            template_length=200,
        ),
        # 10. proper pair R2
        _make_read(
            "read_pair",
            ref_start=1200,
            flag=0x1 | 0x2 | 0x80,
            mate_ref_id=0,
            mate_ref_start=1050,
            template_length=-200,
        ),
    ]

    with pysam.AlignmentFile(str(unsorted_bam), "wb", header=header) as outf:
        for read in reads:
            outf.write(read)

    pysam.sort("-o", str(sorted_bam), str(unsorted_bam))
    pysam.index(str(sorted_bam))

    return sorted_bam

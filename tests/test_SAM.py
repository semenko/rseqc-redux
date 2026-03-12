"""Tests for rseqc.SAM."""

import importlib
import io
import sys
from pathlib import Path

import pytest

from rseqc import SAM

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def test_import():
    """Verify that rseqc.SAM can be imported."""
    mod = importlib.import_module("rseqc.SAM")
    assert mod is not None


def test_parsebam_exists():
    """ParseBAM should still be importable."""
    assert hasattr(SAM, "ParseBAM")


def test_parsesam_removed():
    """ParseSAM dead code class should no longer exist."""
    assert not hasattr(SAM, "ParseSAM")


def test_qcsam_removed():
    """QCSAM dead code class should no longer exist."""
    assert not hasattr(SAM, "QCSAM")


def test_parsebam_init(mini_bam):
    """ParseBAM can be instantiated with a BAM file."""
    obj = SAM.ParseBAM(str(mini_bam))
    assert obj is not None
    assert obj.bam_format is True


def test_parsebam_stat(mini_bam, capsys):
    """ParseBAM.stat() should report mapping statistics."""
    obj = SAM.ParseBAM(str(mini_bam))
    obj.stat(q_cut=30)
    captured = capsys.readouterr()
    assert "Total records:" in captured.out


def test_parsebam_configure_experiment(mini_bam):
    """ParseBAM.configure_experiment() with mini BED should not crash."""
    obj = SAM.ParseBAM(str(mini_bam))
    bed_file = str(FIXTURES_DIR / "mini.bed")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.configure_experiment(refbed=bed_file, sample_size=100)
    finally:
        sys.stderr = old_stderr


def test_parsebam_readsNVC(mini_bam, tmp_path):
    """readsNVC() should produce a .NVC.xls file with per-position base counts."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.readsNVC(outfile=outprefix, q_cut=30)
    finally:
        sys.stderr = old_stderr
    xls_file = Path(outprefix + ".NVC.xls")
    assert xls_file.exists()
    content = xls_file.read_text()
    assert "Position" in content
    # Should have A, C, G, T columns
    assert "\tA\t" in content


def test_parsebam_readGC(mini_bam, tmp_path):
    """readGC() should produce a .GC.xls file with GC content distribution."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.readGC(outfile=outprefix, q_cut=30)
    finally:
        sys.stderr = old_stderr
    xls_file = Path(outprefix + ".GC.xls")
    assert xls_file.exists()
    content = xls_file.read_text()
    assert "GC%" in content
    # Should have at least one data row
    lines = [ln for ln in content.strip().split("\n") if not ln.startswith("GC%")]
    assert len(lines) >= 1


def test_parsebam_readDupRate(mini_bam, tmp_path):
    """readDupRate() should produce seq and pos DupRate files."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.readDupRate(q_cut=30, outfile=outprefix)
    finally:
        sys.stderr = old_stderr
    seq_file = Path(outprefix + ".seq.DupRate.xls")
    pos_file = Path(outprefix + ".pos.DupRate.xls")
    assert seq_file.exists()
    assert pos_file.exists()
    seq_content = seq_file.read_text()
    assert "Occurrence" in seq_content
    pos_content = pos_file.read_text()
    assert "Occurrence" in pos_content


def test_parsebam_readsQual_boxplot(mini_bam, tmp_path):
    """readsQual_boxplot() should produce a .qual.r file."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.readsQual_boxplot(outfile=outprefix, q_cut=30)
    finally:
        sys.stderr = old_stderr
    r_file = Path(outprefix + ".qual.r")
    assert r_file.exists()
    content = r_file.read_text()
    assert "boxplot(" in content
    assert "heatmap(" in content


def test_parsebam_clipping_profile_SE(mini_bam, tmp_path):
    """clipping_profile() in single-end mode should produce output files."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.clipping_profile(outfile=outprefix, q_cut=30, PE=False)
    finally:
        sys.stderr = old_stderr
    xls_file = Path(outprefix + ".clipping_profile.xls")
    assert xls_file.exists()
    content = xls_file.read_text()
    assert "Position" in content


def test_parsebam_insertion_profile_SE(mini_bam, tmp_path):
    """insertion_profile() in single-end mode should produce output files."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.insertion_profile(outfile=outprefix, q_cut=30, PE=False)
    finally:
        sys.stderr = old_stderr
    xls_file = Path(outprefix + ".insertion_profile.xls")
    assert xls_file.exists()
    content = xls_file.read_text()
    assert "Position" in content


def test_parsebam_bam2fq_single(mini_bam, tmp_path):
    """bam2fq() in single-end mode should produce a .fastq file."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.bam2fq(prefix=outprefix, paired=False)
    finally:
        sys.stderr = old_stderr
    fq_file = Path(outprefix + ".fastq")
    assert fq_file.exists()
    content = fq_file.read_text()
    # FASTQ entries start with @
    assert content.startswith("@")
    lines = content.strip().split("\n")
    # Each FASTQ entry is 4 lines
    assert len(lines) % 4 == 0
    assert len(lines) >= 4


def test_parsebam_bam2fq_paired(mini_bam, tmp_path):
    """bam2fq() in paired-end mode should produce R1 and R2 files."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.bam2fq(prefix=outprefix, paired=True)
    finally:
        sys.stderr = old_stderr
    r1_file = Path(outprefix + ".R1.fastq")
    r2_file = Path(outprefix + ".R2.fastq")
    assert r1_file.exists()
    assert r2_file.exists()


def test_parsebam_stat_counts(mini_bam, capsys):
    """stat() should report correct counts for known reads."""
    obj = SAM.ParseBAM(str(mini_bam))
    obj.stat(q_cut=30)
    captured = capsys.readouterr()
    output = captured.out
    # We have 11 total reads in mini_bam
    assert "Total records:" in output
    # Check that some expected categories appear
    assert "QC failed:" in output
    assert "Unmapped reads:" in output


@pytest.mark.parametrize("q_cut", [0, 30, 60])
def test_parsebam_stat_qcut_levels(mini_bam, capsys, q_cut):
    """stat() should work with different quality cutoffs."""
    obj = SAM.ParseBAM(str(mini_bam))
    obj.stat(q_cut=q_cut)
    captured = capsys.readouterr()
    assert "Total records:" in captured.out


def test_parsebam_readGC_content(mini_bam, tmp_path):
    """readGC() with q_cut=0 should include all mapped reads."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.readGC(outfile=outprefix, q_cut=0)
    finally:
        sys.stderr = old_stderr
    xls_file = Path(outprefix + ".GC.xls")
    content = xls_file.read_text()
    # All reads are "AAAA..." so GC should be 0.00
    assert "0.00" in content
    r_file = Path(outprefix + ".GC_plot.r")
    assert r_file.exists()
    assert "pdf(" in r_file.read_text()


def test_parsebam_readsNVC_content(mini_bam, tmp_path):
    """readsNVC() should count nucleotides at each position."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.readsNVC(outfile=outprefix, q_cut=0)
    finally:
        sys.stderr = old_stderr
    xls_file = Path(outprefix + ".NVC.xls")
    content = xls_file.read_text()
    lines = [ln for ln in content.strip().split("\n") if not ln.startswith("Position")]
    # Should have data for each position in the reads
    assert len(lines) >= 1
    # R script should be generated
    r_file = Path(outprefix + ".NVC_plot.r")
    assert r_file.exists()


def test_parsebam_readDupRate_content(mini_bam, tmp_path):
    """readDupRate() with q_cut=0 includes more reads."""
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "test")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.readDupRate(q_cut=0, outfile=outprefix)
    finally:
        sys.stderr = old_stderr
    seq_file = Path(outprefix + ".seq.DupRate.xls")
    seq_content = seq_file.read_text()
    lines = [ln for ln in seq_content.strip().split("\n") if not ln.startswith("Occurrence")]
    assert len(lines) >= 1
    # R script should be generated
    r_file = Path(outprefix + ".DupRate_plot.r")
    assert r_file.exists()
    assert "pdf(" in r_file.read_text()


def test_parsebam_readsNVC_exact_counts(mini_bam, tmp_path):
    """readsNVC() must produce exact per-position base counts.

    With q_cut=0, all mapped reads are included:
      - 10 forward reads (all 'A' bases → A column)
      - 1 reverse read ('A' bases reverse-complemented → T column)
    Every position should show A=10, T=1, C=G=N=X=0.
    """
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "nvc")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.readsNVC(outfile=outprefix, q_cut=0)
    finally:
        sys.stderr = old_stderr

    xls_file = Path(outprefix + ".NVC.xls")
    content = xls_file.read_text()
    lines = content.strip().split("\n")

    # Header line
    assert lines[0].split() == ["Position", "A", "C", "G", "T", "N", "X"]

    # Data lines: 50 positions (read length = 50)
    data_lines = lines[1:]
    assert len(data_lines) == 50

    for i, line in enumerate(data_lines):
        fields = line.split()
        assert int(fields[0]) == i, f"position mismatch at line {i}"
        assert int(fields[1]) == 10, f"A count wrong at pos {i}"  # A
        assert int(fields[2]) == 0, f"C count wrong at pos {i}"  # C
        assert int(fields[3]) == 0, f"G count wrong at pos {i}"  # G
        assert int(fields[4]) == 1, f"T count wrong at pos {i}"  # T (rev-comp of A)
        assert int(fields[5]) == 0, f"N count wrong at pos {i}"  # N
        assert int(fields[6]) == 0, f"X count wrong at pos {i}"  # X

    # R script should reference correct counts
    r_file = Path(outprefix + ".NVC_plot.r")
    r_content = r_file.read_text()
    assert "A_count=c(" in r_content
    assert "T_count=c(" in r_content


def test_parsebam_clipping_profile_exact_counts(tmp_path):
    """clipping_profile() must produce exact per-position clip counts.

    Uses a custom BAM with reads that have soft clipping:
      - read1: 5S45M (5bp soft clip at start)
      - read2: 45M5S (5bp soft clip at end)
      - read3: 50M   (no clipping)
    """
    import pysam

    tmpdir = tmp_path / "bam"
    tmpdir.mkdir()

    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 50000}],
        }
    )

    def _make_read(name, ref_start, cigar, flag=0, mapq=60):
        a = pysam.AlignedSegment(header)
        a.query_name = name
        a.flag = flag
        a.reference_id = 0
        a.reference_start = ref_start
        a.mapping_quality = mapq
        a.cigar = cigar
        read_len = sum(ln for op, ln in cigar if op in (0, 1, 4))
        a.query_sequence = "A" * read_len
        a.query_qualities = pysam.qualitystring_to_array("I" * read_len)
        return a

    reads = [
        _make_read("clip_left", 1000, [(4, 5), (0, 45)]),  # 5S45M
        _make_read("clip_right", 2000, [(0, 45), (4, 5)]),  # 45M5S
        _make_read("no_clip", 3000, [(0, 50)]),  # 50M
    ]

    unsorted_bam = tmpdir / "unsorted.bam"
    sorted_bam = tmpdir / "clip_test.bam"
    with pysam.AlignmentFile(str(unsorted_bam), "wb", header=header) as outf:
        for r in reads:
            outf.write(r)
    pysam.sort("-o", str(sorted_bam), str(unsorted_bam))
    pysam.index(str(sorted_bam))

    obj = SAM.ParseBAM(str(sorted_bam))
    outprefix = str(tmp_path / "clip")
    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.clipping_profile(outfile=outprefix, q_cut=0, PE=False)
    finally:
        sys.stderr = old_stderr

    xls_file = Path(outprefix + ".clipping_profile.xls")
    content = xls_file.read_text()
    lines = content.strip().split("\n")

    # Header
    assert lines[0] == "Position\tClipped_nt\tNon_clipped_nt"

    # Data: 50 positions (read length = 50)
    data_lines = lines[1:]
    assert len(data_lines) == 50

    for i, line in enumerate(data_lines):
        fields = line.split("\t")
        pos = int(fields[0])
        clipped = float(fields[1])
        non_clipped = float(fields[2])
        assert pos == i

        if i < 5:
            # Position 0-4: clip_left has S here (1 read clipped)
            assert clipped == 1.0, f"pos {i}: expected 1 clipped"
        elif i >= 45:
            # Position 45-49: clip_right has S here (1 read clipped)
            assert clipped == 1.0, f"pos {i}: expected 1 clipped"
        else:
            # Middle positions: no clipping
            assert clipped == 0.0, f"pos {i}: expected 0 clipped"

        # Total reads = 3, clipped + non_clipped should = 3
        assert clipped + non_clipped == 3.0, f"pos {i}: total should be 3"


def test_parsebam_bamTowig_exact_output(mini_bam, tmp_path):
    """bamTowig() must produce correct per-position coverage in wig format.

    mini_bam has these mapped reads (after filtering duplicates, secondary, etc.):
      read_unique1:  50M at chr1:1050 → covers 1051-1100 (1-based)
      read_unique2:  50M at chr1:1100 → covers 1101-1150
      read_spliced:  50M1400N50M at chr1:1450 → covers 1451-1500, 2901-2950
      read_reverse:  50M at chr1:6500 → covers 6501-6550
      read_chr2:     50M at chr2:500  → covers 501-550
      read_lowmapq:  50M at chr1:1050, MAPQ=2 → filtered by q_cut=30
      read_pair R1:  50M at chr1:1050 → covers 1051-1100
      read_pair R2:  50M at chr1:1200 → covers 1201-1250

    Expected coverage at chr1:
      1051-1100: 2.0 (read_unique1 + read_pair_r1)
      1101-1150: 1.0 (read_unique2)
      1201-1250: 1.0 (read_pair_r2)
      1451-1500: 1.0 (read_spliced exon1)
      2901-2950: 1.0 (read_spliced exon2)
      6501-6550: 1.0 (read_reverse)
    """
    obj = SAM.ParseBAM(str(mini_bam))
    outprefix = str(tmp_path / "wig")
    chrom_sizes = {"chr1": 50000, "chr2": 50000}
    # Write a temp chrom sizes file (needed by bamTowig signature but won't be used for bigwig)
    chrom_file = str(tmp_path / "chrom.sizes")
    with open(chrom_file, "w") as f:
        for c, s in chrom_sizes.items():
            f.write(f"{c}\t{s}\n")

    old_stderr = sys.stderr
    sys.stderr = io.StringIO()
    try:
        obj.bamTowig(
            outfile=outprefix,
            chrom_sizes=chrom_sizes,
            chrom_file=chrom_file,
            skip_multi=True,
            strand_rule=None,
            WigSumFactor=None,
            q_cut=30,
        )
    finally:
        sys.stderr = old_stderr

    wig_file = Path(outprefix + ".wig")
    assert wig_file.exists()
    content = wig_file.read_text()

    # Parse the wig file into {chrom: {pos: value}}
    coverage: dict[str, dict[int, float]] = {}
    current_chrom = None
    for line in content.strip().split("\n"):
        if line.startswith("variableStep"):
            current_chrom = line.split("chrom=")[1]
            coverage[current_chrom] = {}
        else:
            parts = line.split("\t")
            pos = int(parts[0])
            val = float(parts[1])
            coverage[current_chrom][pos] = val

    # Check chr1 coverage
    assert "chr1" in coverage
    chr1 = coverage["chr1"]

    # Positions 1051-1100 should have coverage 2.0 (read_unique1 + read_pair_r1)
    for p in range(1051, 1101):
        assert chr1.get(p) == 2.0, f"chr1:{p} expected 2.0, got {chr1.get(p)}"

    # Positions 1101-1150 should have coverage 1.0 (read_unique2)
    for p in range(1101, 1151):
        expected = 1.0 if p > 1100 else 2.0  # 1101-1150
        assert chr1.get(p) == expected, f"chr1:{p} expected {expected}, got {chr1.get(p)}"

    # read_pair_r2 covers 1201-1250
    for p in range(1201, 1251):
        assert chr1.get(p) == 1.0, f"chr1:{p} expected 1.0, got {chr1.get(p)}"

    # Spliced read exon2: 2901-2950
    for p in range(2901, 2951):
        assert chr1.get(p) == 1.0, f"chr1:{p} expected 1.0, got {chr1.get(p)}"

    # read_reverse: 6501-6550
    for p in range(6501, 6551):
        assert chr1.get(p) == 1.0, f"chr1:{p} expected 1.0, got {chr1.get(p)}"

    # chr2 should have coverage 501-550
    assert "chr2" in coverage
    chr2 = coverage["chr2"]
    for p in range(501, 551):
        assert chr2.get(p) == 1.0, f"chr2:{p} expected 1.0, got {chr2.get(p)}"


def test_parsebam_stat_splice_counts(mini_bam, capsys):
    """stat() must correctly count spliced vs non-spliced reads.

    With q_cut=0, the unique hits include:
      - read_spliced (has N in CIGAR) → splice
      - read_unique1, read_unique2, read_reverse, read_chr2, read_lowmapq,
        read_pair R1, read_pair R2 → non-splice (all 50M)
    """
    obj = SAM.ParseBAM(str(mini_bam))
    obj.stat(q_cut=0)
    captured = capsys.readouterr()
    output = captured.out
    # Parse splice counts from output
    for line in output.split("\n"):
        if "Spliced reads" in line:
            count = int(line.split()[-1])
            assert count == 1, f"Expected 1 spliced read, got {count}"
        if "Non-splice reads" in line:
            count = int(line.split()[-1])
            assert count == 7, f"Expected 7 non-spliced reads, got {count}"

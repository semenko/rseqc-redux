"""Tests for helper functions in scripts/ that can be unit-tested without BAM files."""

import math

from bx.intervals import Intersecter, Interval

from scripts.FPKM_count import build_range
from scripts.geneBody_coverage import pearson_moment_coefficient, valid_name
from scripts.junction_annotation import generate_bed12, generate_interact
from scripts.read_distribution import cal_size, foundone
from scripts.read_hexamer import file_exist
from scripts.RNA_fragment_size import overlap_length2
from scripts.RPKM_saturation import square_error
from scripts.tin import shannon_entropy, tin_score, uniqify

# --- scripts.RPKM_saturation: square_error ---


def test_square_error_basic():
    lst = [2.0, 4.0, 6.0, 8.0, 10.0]
    result = square_error(lst)
    assert result is not None
    assert len(result) == 5
    # Last element should be 0 (true_rpkm == lst[-1])
    assert result[-1] == 0.0
    # Each SE = abs(i - true_rpkm) / true_rpkm
    assert abs(result[0] - abs(2.0 - 10.0) / 10.0) < 1e-10


def test_square_error_zero_true_rpkm():
    """If last element is 0, return None."""
    lst = [1.0, 2.0, 0.0]
    assert square_error(lst) is None


def test_square_error_constant():
    """If all values are same, range is 0 → None."""
    lst = [5.0, 5.0, 5.0]
    assert square_error(lst) is None


def test_square_error_two_elements():
    lst = [5.0, 10.0]
    result = square_error(lst)
    assert result is not None
    assert result[0] == 0.5  # abs(5-10)/10
    assert result[1] == 0.0  # abs(10-10)/10


# --- scripts.geneBody_coverage: valid_name, pearson_moment_coefficient ---


def test_valid_name_simple():
    assert valid_name("sample1") == "sample1"


def test_valid_name_with_spaces():
    result = valid_name("my sample")
    assert result == "my_sample"


def test_valid_name_starts_with_digit():
    result = valid_name("1sample")
    assert result.startswith("V")


def test_valid_name_special_chars():
    result = valid_name("sample@#$.bam")
    assert "@" not in result
    assert "#" not in result
    assert "$" not in result
    # Period should be preserved (valid R char)
    assert ".bam" in result


def test_valid_name_multiple_spaces():
    result = valid_name("my  sample  name")
    assert "_" in result


def test_pearson_moment_coefficient_symmetric():
    """Symmetric distribution should have skewness near 0."""
    lst = [1.0, 2.0, 3.0, 4.0, 5.0]
    result = pearson_moment_coefficient(lst)
    assert isinstance(result, float)


def test_pearson_moment_coefficient_right_skewed():
    """Right-skewed distribution should have positive skewness."""
    lst = [1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 10.0, 10.0, 100.0]
    result = pearson_moment_coefficient(lst)
    assert isinstance(result, float)


# --- scripts.read_distribution: cal_size, foundone ---


def test_cal_size_basic():
    entries = [["chr1", 100, 200], ["chr1", 300, 500]]
    assert cal_size(entries) == 300  # 100 + 200


def test_cal_size_empty():
    assert cal_size([]) == 0


def test_cal_size_single():
    assert cal_size([["chr1", 0, 1000]]) == 1000


def test_foundone_found():
    ranges = {"CHR1": Intersecter()}
    ranges["CHR1"].add_interval(Interval(100, 200))
    assert foundone("CHR1", ranges, 150, 151) > 0


def test_foundone_not_found():
    ranges = {"CHR1": Intersecter()}
    ranges["CHR1"].add_interval(Interval(100, 200))
    assert foundone("CHR1", ranges, 300, 301) == 0


def test_foundone_missing_chrom():
    ranges = {"CHR1": Intersecter()}
    ranges["CHR1"].add_interval(Interval(100, 200))
    assert foundone("CHR2", ranges, 150, 151) == 0


def test_foundone_empty_ranges():
    assert foundone("CHR1", {}, 100, 200) == 0


# --- scripts.RNA_fragment_size: overlap_length2 ---


def test_overlap_length2_full_overlap():
    lst1 = [[100, 200]]
    lst2 = [[100, 200]]
    result = overlap_length2(lst1, lst2)
    assert result == 101  # range(100, 201) = 101 positions


def test_overlap_length2_partial_overlap():
    lst1 = [[100, 200]]
    lst2 = [[150, 250]]
    result = overlap_length2(lst1, lst2)
    assert result == 51  # range(150, 201) = 51 positions


def test_overlap_length2_no_overlap():
    lst1 = [[100, 200]]
    lst2 = [[300, 400]]
    result = overlap_length2(lst1, lst2)
    assert result == 0


def test_overlap_length2_multiple_exons():
    lst1 = [[100, 200], [300, 400]]
    lst2 = [[150, 350]]
    result = overlap_length2(lst1, lst2)
    # exon1 overlaps [150,200]: 51 positions
    # exon2 overlaps [300,350]: 51 positions
    assert result == 51 + 51


# --- scripts.tin: uniqify, shannon_entropy, tin_score ---


def test_uniqify_basic():
    assert uniqify([1, 2, 2, 3, 3, 4]) == [1, 2, 3, 4]


def test_uniqify_no_duplicates():
    assert uniqify([1, 2, 3]) == [1, 2, 3]


def test_uniqify_all_same():
    assert uniqify([5, 5, 5]) == [5]


def test_uniqify_empty():
    assert uniqify([]) == []


def test_uniqify_preserves_order():
    assert uniqify([3, 1, 2, 1, 3]) == [3, 1, 2]


def test_shannon_entropy_uniform():
    """Uniform distribution should have max entropy."""
    arg = [1.0, 1.0, 1.0, 1.0]
    result = shannon_entropy(arg)
    assert abs(result - math.log(4)) < 1e-10


def test_shannon_entropy_single():
    """Single element: P=1, log(1)=0, entropy=0."""
    assert shannon_entropy([5.0]) == 0


def test_shannon_entropy_two_equal():
    arg = [3.0, 3.0]
    result = shannon_entropy(arg)
    assert abs(result - math.log(2)) < 1e-10


def test_shannon_entropy_skewed():
    """Skewed distribution should have lower entropy than uniform."""
    uniform = shannon_entropy([1.0, 1.0, 1.0, 1.0])
    skewed = shannon_entropy([10.0, 1.0, 1.0, 1.0])
    assert skewed < uniform


def test_tin_score_empty():
    assert tin_score([], 100) == 0


def test_tin_score_uniform():
    """Uniform coverage → TIN close to 100."""
    cvg = [10.0] * 100
    result = tin_score(cvg, 100)
    assert abs(result - 100.0) < 1e-10


def test_tin_score_with_zeros():
    """Zeros in coverage are filtered before entropy calculation."""
    cvg = [10.0, 0.0, 10.0, 0.0, 10.0]
    result = tin_score(cvg, 5)
    # 3 equal values → entropy = ln(3), TIN = 100*exp(ln(3))/5 = 100*3/5 = 60
    assert abs(result - 60.0) < 1e-10


def test_tin_score_single_position():
    """Single non-zero position → entropy=0 → TIN = 100*1/length."""
    cvg = [10.0]
    result = tin_score(cvg, 100)
    assert abs(result - 1.0) < 1e-10


# --- scripts.FPKM_count: build_range ---


def test_build_range_basic(tmp_path):
    """build_range should parse BED12 and build interval trees."""
    bed = tmp_path / "test.bed"
    bed.write_text("chr1\t1000\t5000\tgene1\t0\t+\t1000\t5000\t0\t2\t500,500,\t0,3500,\n")
    result = build_range(str(bed))
    assert "CHR1" in result  # build_range uppercases chrom
    # Should have intervals in the tree
    assert len(result["CHR1"].find(1000, 1500)) > 0  # first exon
    assert len(result["CHR1"].find(4500, 5000)) > 0  # second exon


def test_build_range_skips_comments(tmp_path):
    bed = tmp_path / "test.bed"
    bed.write_text("# comment\ntrack name=test\nchr1\t1000\t5000\tgene1\t0\t+\t1000\t5000\t0\t1\t4000,\t0,\n")
    result = build_range(str(bed))
    assert "CHR1" in result


def test_build_range_skips_bad_lines(tmp_path):
    bed = tmp_path / "test.bed"
    bed.write_text("bad line\nchr1\t1000\t5000\tgene1\t0\t+\t1000\t5000\t0\t1\t4000,\t0,\n")
    result = build_range(str(bed))
    assert "CHR1" in result


def test_build_range_empty_file(tmp_path):
    bed = tmp_path / "test.bed"
    bed.write_text("")
    result = build_range(str(bed))
    assert result == {}


# --- scripts.read_hexamer: file_exist ---


def test_file_exist_true(tmp_path):
    f = tmp_path / "test.txt"
    f.write_text("hello")
    assert file_exist(str(f)) is True


def test_file_exist_false():
    assert file_exist("/nonexistent/file.txt") is False


# --- scripts.junction_annotation: generate_bed12, generate_interact ---

JUNCTION_XLS_CONTENT = (
    "chrom\tintron_st(0-based)\tintron_end(1-based)\tread_count\tannotation\n"
    "chr1\t1000\t2000\t50\tannotated\n"
    "chr1\t3000\t4000\t25\tpartial_novel\n"
    "chr1\t5000\t6000\t10\tcomplete_novel\n"
)


def test_generate_bed12_basic(tmp_path):
    xls = tmp_path / "test.junction.xls"
    xls.write_text(JUNCTION_XLS_CONTENT)
    generate_bed12(str(xls))
    bed = tmp_path / "test.junction.bed"
    assert bed.exists()
    lines = bed.read_text().strip().split("\n")
    assert len(lines) == 3
    # Check first line is annotated (red)
    fields = lines[0].split("\t")
    assert fields[0] == "chr1"
    assert fields[8] == "205,0,0"  # annotated = red
    assert fields[9] == "2"  # blockCount


def test_generate_bed12_colors(tmp_path):
    xls = tmp_path / "test.junction.xls"
    xls.write_text(JUNCTION_XLS_CONTENT)
    generate_bed12(str(xls))
    bed = tmp_path / "test.junction.bed"
    lines = bed.read_text().strip().split("\n")
    colors = [line.split("\t")[8] for line in lines]
    assert colors == ["205,0,0", "0,205,0", "0,0,205"]


def test_generate_bed12_empty(tmp_path):
    xls = tmp_path / "test.junction.xls"
    xls.write_text("chrom\tintron_st(0-based)\tintron_end(1-based)\tread_count\tannotation\n")
    generate_bed12(str(xls))
    bed = tmp_path / "test.junction.bed"
    assert bed.exists()
    assert bed.read_text().strip() == ""


def test_generate_interact_basic(tmp_path):
    xls = tmp_path / "test.junction.xls"
    xls.write_text(JUNCTION_XLS_CONTENT)
    generate_interact(str(xls), "sample.bam")
    interact = tmp_path / "test.junction.Interact.bed"
    assert interact.exists()
    lines = interact.read_text().strip().split("\n")
    assert lines[0].startswith("track type=interact")
    assert len(lines) == 4  # 1 header + 3 junctions


def test_generate_interact_colors(tmp_path):
    xls = tmp_path / "test.junction.xls"
    xls.write_text(JUNCTION_XLS_CONTENT)
    generate_interact(str(xls), "sample.bam")
    interact = tmp_path / "test.junction.Interact.bed"
    lines = interact.read_text().strip().split("\n")[1:]  # skip header
    colors = [line.split("\t")[7] for line in lines]
    assert colors == ["205,0,0", "0,205,0", "0,0,205"]


def test_generate_interact_bam_in_header(tmp_path):
    xls = tmp_path / "test.junction.xls"
    xls.write_text(JUNCTION_XLS_CONTENT)
    generate_interact(str(xls), "my_sample.bam")
    interact = tmp_path / "test.junction.Interact.bed"
    header = interact.read_text().split("\n")[0]
    assert "my_sample.bam" in header

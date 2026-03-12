"""Tests for rseqc.fasta."""

import importlib

from rseqc.fasta import Fasta


def test_import():
    """Verify that rseqc.fasta can be imported."""
    mod = importlib.import_module("rseqc.fasta")
    assert mod is not None


def test_fasta_from_file(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\nACGTACGT\n>chr2\nTTTTAAAA\n")
    obj = Fasta(str(fa))
    assert "chr1" in obj.seqs
    assert "chr2" in obj.seqs
    assert obj.seqs["chr1"] == "ACGTACGT"
    assert obj.seqs["chr2"] == "TTTTAAAA"


def test_fasta_get_names(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGT\n>seq2\nTGCA\n")
    obj = Fasta(str(fa))
    assert obj.getNames() == ["seq1", "seq2"]


def test_fasta_get_seq(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGT\n")
    obj = Fasta(str(fa))
    assert obj.getSeq("seq1") == "ACGT"


def test_fasta_get_all_seqs(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGT\n>seq2\nTGCA\n")
    obj = Fasta(str(fa))
    seqs = obj.getSeq()
    assert len(seqs) == 2


def test_fasta_add_seq():
    obj = Fasta()
    obj.addSeq("test", "acgt")
    assert obj.seqs["test"] == "ACGT"


def test_fasta_rev_comp(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGT\n")
    obj = Fasta(str(fa))
    result = obj.revComp("seq1")
    assert result == "ACGT"  # ACGT is its own reverse complement


def test_fasta_rev_comp_asymmetric(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nAAAA\n")
    obj = Fasta(str(fa))
    result = obj.revComp("seq1")
    assert result == "TTTT"


def test_fasta_fetch_seq(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\nACGTACGTACGT\n")
    obj = Fasta(str(fa))
    result = obj.fetchSeq(chr="chr1", st=2, end=6)
    assert result == "GTAC"


def test_fasta_multiline(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\nACGT\nTGCA\n")
    obj = Fasta(str(fa))
    assert obj.seqs["chr1"] == "ACGTTGCA"


def test_fasta_get_seq_len(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\nACGTACGT\n")
    obj = Fasta(str(fa))
    lengths = obj.getSeqLen()
    assert lengths["chr1"] == 8


def test_fasta_empty_file(tmp_path):
    """Bug #2 regression: empty FASTA should not raise UnboundLocalError."""
    fa = tmp_path / "empty.fa"
    fa.write_text("")
    obj = Fasta(str(fa))
    assert obj.seqs == {}
    assert obj.IDs == []


def test_fasta_no_header(tmp_path):
    """Bug #2 regression: FASTA with no header lines should handle gracefully."""
    fa = tmp_path / "noheader.fa"
    fa.write_text("ACGTACGT\n")
    obj = Fasta(str(fa))
    assert obj.seqs == {}
    assert obj.IDs == []


def test_fasta_add_seq_duplicate(capsys):
    """Adding a duplicate ID prints warning and does not overwrite."""
    obj = Fasta()
    obj.addSeq("test", "ACGT")
    obj.addSeq("test", "TTTT")
    assert obj.seqs["test"] == "ACGT"  # not overwritten
    captured = capsys.readouterr()
    assert "already exists" in captured.err


def test_fasta_get_seq_len_single(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\nACGTACGT\n>chr2\nTTTT\n")
    obj = Fasta(str(fa))
    result = obj.getSeqLen("chr1")
    assert result["chr1"] == 8
    assert "chr2" not in result


def test_fasta_get_seq_len_missing(capsys):
    obj = Fasta()
    obj.addSeq("chr1", "ACGT")
    result = obj.getSeqLen("missing")
    assert "missing" not in result or result.get("missing") is None
    captured = capsys.readouterr()
    assert "Not found" in captured.err


def test_fasta_count_base_default(tmp_path, capsys):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nAACCGGTT\n")
    obj = Fasta(str(fa))
    obj.countBase()
    captured = capsys.readouterr()
    assert "seq1" in captured.out
    assert "8" in captured.out  # total length


def test_fasta_count_base_pattern(tmp_path, capsys):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGTACGT\n")
    obj = Fasta(str(fa))
    obj.countBase(pattern="ACG")
    captured = capsys.readouterr()
    assert "2" in captured.out  # ACG appears twice


def test_fasta_cal_entropy():
    obj = Fasta()
    obj.addSeq("test", "ACGTACGTACGTACGT")
    results = list(obj.cal_entropy(length=1))
    assert len(results) == 1
    name, entropy = results[0]
    assert name == "test"
    assert entropy > 0  # equal base distribution should give max entropy


def test_fasta_cal_entropy_homopolymer():
    obj = Fasta()
    obj.addSeq("test", "AAAAAAAAAAAA")
    results = list(obj.cal_entropy(length=1))
    name, entropy = results[0]
    # Only A present (prop=1, info=0), C/G/T have 0 count → skipped → entropy=0
    assert entropy == 0.0


def test_fasta_rev_comp_all_seqs(tmp_path):
    """revComp with no seqID prints and returns first seq's rev comp."""
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nAAAA\n")
    obj = Fasta(str(fa))
    result = obj.revComp()
    assert result == "TTTT"


def test_fasta_get_uniq_seqs(capsys):
    obj = Fasta()
    obj.addSeq("a", "ACGT")
    obj.addSeq("b", "TGCA")
    obj.addSeq("c", "ACGT")  # duplicate of a
    obj.getUniqSeqs()
    captured = capsys.readouterr()
    # Should print unique seqs with count
    assert "ACGT" in captured.out
    assert "TGCA" in captured.out
    assert "_2" in captured.out  # ACGT appears twice


def test_fasta_print_seqs(tmp_path, capsys):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGT\n")
    obj = Fasta(str(fa))
    obj.printSeqs(n=2)
    captured = capsys.readouterr()
    assert ">seq1" in captured.out
    assert "AC" in captured.out
    assert "GT" in captured.out


def test_fasta_find_pattern(tmp_path):
    obj = Fasta()
    obj.addSeq("chr1", "AACGTACGTAA")
    outfile = str(tmp_path / "out.bed")
    obj.findPattern("CGT", outfile, rev=False)
    with open(outfile) as f:
        lines = f.readlines()
    assert len(lines) == 2  # CGT at positions 2 and 6
    assert "chr1" in lines[0]


def test_fasta_find_pattern_with_rev(tmp_path):
    obj = Fasta()
    obj.addSeq("chr1", "AACGTACGTAA")
    outfile = str(tmp_path / "out.bed")
    obj.findPattern("CGT", outfile, rev=True)
    with open(outfile) as f:
        lines = f.readlines()
    # Forward CGT matches + reverse complement ACG matches
    assert len(lines) >= 2


def test_fasta_find_pattern_specific_seq(tmp_path):
    obj = Fasta()
    obj.addSeq("chr1", "AACGTAA")
    obj.addSeq("chr2", "TTTTTT")
    outfile = str(tmp_path / "out.bed")
    obj.findPattern("CGT", outfile, seqID="chr1", rev=False)
    with open(outfile) as f:
        content = f.read()
    assert "chr1" in content


def test_fasta_fetch_seq_missing_chrom():
    obj = Fasta()
    obj.addSeq("chr1", "ACGTACGT")
    result = obj.fetchSeq(chr="chrZ", st=0, end=5)
    assert result == ""


def test_fasta_fetch_seq_from_file(tmp_path):
    obj = Fasta()
    obj.addSeq("chr1", "ACGTACGTACGT")
    bedfile = tmp_path / "input.bed"
    bedfile.write_text("chr1\t2\t6\n")
    outfile = str(tmp_path / "out.fa")
    result = obj.fetchSeq(infile=str(bedfile), outfile=outfile)
    assert result is None  # returns None when using infile
    with open(outfile) as f:
        content = f.read()
    assert "GTAC" in content


def test_fasta_fetch_seq_from_file_with_strand(tmp_path):
    obj = Fasta()
    obj.addSeq("chr1", "ACGTACGTACGT")
    bedfile = tmp_path / "input.bed"
    # BED6 with minus strand
    bedfile.write_text("chr1\t0\t4\tname\t0\t-\n")
    outfile = str(tmp_path / "out.fa")
    obj.fetchSeq(infile=str(bedfile), outfile=outfile)
    with open(outfile) as f:
        content = f.read()
    assert "strand=-" in content


def test_fasta_lowercase_to_upper(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nacgt\n")
    obj = Fasta(str(fa))
    assert obj.seqs["seq1"] == "ACGT"

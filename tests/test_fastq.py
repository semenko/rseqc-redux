"""Tests for rseqc.fastq."""

import bz2
import gzip
import importlib

import pandas as pd

from rseqc import fastq


def test_import():
    """Verify that rseqc.fastq can be imported."""
    mod = importlib.import_module("rseqc.fastq")
    assert mod is not None


# --- _open_file ---


def test_open_file_plain_text(tmp_path):
    f = tmp_path / "test.txt"
    f.write_text("line1\nline2\nline3\n")
    lines = list(fastq._open_file(str(f)))
    assert lines == ["line1", "line2", "line3"]


def test_open_file_plain_text_empty(tmp_path):
    f = tmp_path / "empty.txt"
    f.write_text("")
    lines = list(fastq._open_file(str(f)))
    assert lines == []


def test_open_file_strips_carriage_return(tmp_path):
    f = tmp_path / "crlf.txt"
    f.write_bytes(b"line1\r\nline2\r\n")
    lines = list(fastq._open_file(str(f)))
    assert lines == ["line1", "line2"]


def test_open_file_gzip(tmp_path):
    f = tmp_path / "test.gz"
    with gzip.open(str(f), "wb") as gz:
        gz.write(b"gzline1\ngzline2\n")
    lines = list(fastq._open_file(str(f)))
    assert lines == ["gzline1", "gzline2"]


def test_open_file_bz2(tmp_path):
    f = tmp_path / "test.bz2"
    with bz2.open(str(f), "wb") as bz:
        bz.write(b"line1\nline2\n")
    lines = list(fastq._open_file(str(f)))
    assert lines == ["line1", "line2"]


# --- fasta_iter ---


def test_fasta_iter(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGTACGT\n>seq2\nTTTTAAAA\n")
    seqs = list(fastq.fasta_iter(str(fa)))
    assert seqs == ["ACGTACGT", "TTTTAAAA"]


def test_fasta_iter_empty_lines(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\n\nACGT\n\n")
    seqs = list(fastq.fasta_iter(str(fa)))
    assert seqs == ["ACGT"]


# --- fastq_iter ---


def test_fastq_iter_seq(tmp_path):
    fq = tmp_path / "test.fq"
    fq.write_text("@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTTTTAAAA\n+\nIIIIIIII\n")
    seqs = list(fastq.fastq_iter(str(fq), mode="seq"))
    assert seqs == ["ACGTACGT", "TTTTAAAA"]


def test_fastq_iter_qual(tmp_path):
    fq = tmp_path / "test.fq"
    fq.write_text("@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTTTTAAAA\n+\nJJJJJJJJ\n")
    quals = list(fastq.fastq_iter(str(fq), mode="qual"))
    assert quals == ["IIIIIIII", "JJJJJJJJ"]


# --- seq2countMat ---


def test_seq2countmat_basic(tmp_path):
    fq = tmp_path / "test.fq"
    fq.write_text("@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n")
    s_obj = fastq.fastq_iter(str(fq), mode="seq")
    mat = fastq.seq2countMat(s_obj, limit=None)
    assert len(mat) == 4  # 4 positions
    assert mat[0]["A"] == 2
    assert mat[1]["C"] == 2


# --- write_matrix_csv ---


class TestWriteMatrixCsv:
    """Verify write_matrix_csv matches original pandas-based output."""

    def test_seq_count_matrix_matches_pandas(self, tmp_path):
        """sc_seqLogo default call: write_matrix_csv(mat, outfile, index_label="Index").

        The old pandas code did NOT sort columns — column order depended on dict
        insertion order, making it non-deterministic. The new code sorts columns
        alphabetically, which is strictly better. We verify that row labels, column
        labels, and all cell values match (ignoring column order).
        """
        mat = {
            0: {"A": 10, "C": 5, "G": 3, "T": 2},
            1: {"A": 1, "C": 8, "T": 11},
            2: {"A": 7, "G": 6, "N": 1, "T": 6},
        }

        # Reproduce old pandas path: from_dict → fillna → .T → to_csv
        old_df = pd.DataFrame.from_dict(mat).fillna(0).astype(int)
        old_df = old_df.T
        # Sort columns to match new code's deterministic behavior
        old_df = old_df[sorted(old_df.columns)]
        pandas_file = str(tmp_path / "pandas.csv")
        old_df.to_csv(pandas_file, index=True, index_label="Index")

        new_file = str(tmp_path / "new.csv")
        fastq.write_matrix_csv(mat, new_file, index_label="Index")

        with open(pandas_file) as f:
            pandas_content = f.read()
        with open(new_file) as f:
            new_content = f.read()

        assert new_content == pandas_content

    def test_qual_count_matrix_matches_pandas(self, tmp_path):
        """sc_seqQual count call: transpose=True, sort_index_descending=True."""
        mat = {
            0: {30: 5, 35: 10, 40: 3},
            1: {30: 2, 40: 8},
            2: {35: 7, 40: 6},
        }

        # Old pandas path for sc_seqQual count CSV:
        # qual2countMat returned pd.DataFrame.from_dict(mat).fillna(0).T
        #   → rows=positions, cols=q_scores
        # sc_seqQual then did:
        #   qual_mat = qual_mat.T → rows=q_scores, cols=positions
        #   qual_mat.sort_index(ascending=False)
        #   qual_mat.to_csv(outfile, index=True, index_label="Index")
        old_df = pd.DataFrame.from_dict(mat).fillna(0).astype(int)
        old_df.index.name = "pos"
        old_df = old_df.T  # rows=positions, cols=q_scores
        old_df = old_df.T  # back to rows=q_scores, cols=positions
        old_df.sort_index(inplace=True, ascending=False)
        pandas_file = str(tmp_path / "pandas.csv")
        old_df.to_csv(pandas_file, index=True, index_label="Index")

        new_file = str(tmp_path / "new.csv")
        fastq.write_matrix_csv(
            mat,
            new_file,
            index_label="Index",
            transpose=True,
            sort_index_descending=True,
        )

        with open(pandas_file) as f:
            pandas_content = f.read()
        with open(new_file) as f:
            new_content = f.read()

        assert new_content == pandas_content

    def test_qual_percent_matrix_matches_pandas(self, tmp_path):
        """sc_seqQual percent call: transpose=True, sort_index_descending=True, normalize=True."""
        mat = {
            0: {30: 5, 35: 10, 40: 3},
            1: {30: 2, 40: 8},
            2: {35: 7, 40: 6},
        }

        # Old pandas path for sc_seqQual percent CSV:
        # Same as count, but then: qual_mat_per = qual_mat / qual_mat.sum()
        old_df = pd.DataFrame.from_dict(mat).fillna(0)
        old_df.index.name = "pos"
        old_df = old_df.T.T  # round-trip through both transposes
        old_df.sort_index(inplace=True, ascending=False)
        old_df_per = old_df / old_df.sum()
        pandas_file = str(tmp_path / "pandas.csv")
        old_df_per.to_csv(pandas_file, index=True, index_label="Index")

        new_file = str(tmp_path / "new.csv")
        fastq.write_matrix_csv(
            mat,
            new_file,
            index_label="Index",
            transpose=True,
            sort_index_descending=True,
            normalize=True,
        )

        with open(pandas_file) as f:
            pandas_lines = f.read().strip().split("\n")
        with open(new_file) as f:
            new_lines = f.read().strip().split("\n")

        # Headers must match exactly
        assert new_lines[0] == pandas_lines[0]
        # Same number of rows
        assert len(new_lines) == len(pandas_lines)
        # Values must match (compare as floats to handle formatting differences)
        for pandas_line, new_line in zip(pandas_lines[1:], new_lines[1:]):
            p_parts = pandas_line.split(",")
            n_parts = new_line.split(",")
            assert p_parts[0] == n_parts[0]  # row label
            for p_val, n_val in zip(p_parts[1:], n_parts[1:]):
                assert abs(float(p_val) - float(n_val)) < 1e-12

    def test_non_transpose_basic(self, tmp_path):
        """Non-transpose mode: rows=outer keys, cols=inner keys."""
        mat = {0: {"A": 3, "C": 1}, 1: {"A": 2, "T": 5}}
        outfile = str(tmp_path / "basic.csv")
        fastq.write_matrix_csv(mat, outfile, index_label="Index")

        with open(outfile) as f:
            lines = f.read().strip().split("\n")
        assert lines[0] == "Index,A,C,T"
        assert lines[1] == "0,3,1,0"
        assert lines[2] == "1,2,0,5"

    def test_empty_mat(self, tmp_path):
        """Empty dict should produce header-only CSV."""
        outfile = str(tmp_path / "empty.csv")
        fastq.write_matrix_csv({}, outfile, index_label="Index")

        with open(outfile) as f:
            content = f.read().strip()
        assert content == "Index"

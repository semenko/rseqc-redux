"""
Created on Sun Aug  2 15:43:45 2020

@author: m102324
"""

from __future__ import annotations

import bz2
import collections
import csv
import gzip
import logging
import sys
from collections.abc import Generator, Iterable
from typing import Any

import logomaker
import matplotlib.pyplot as plt


def _open_file(path: str) -> Generator[str, None, None]:
    """Open a plain, gzip, or bz2 file and yield stripped lines."""
    fh: Any
    if path.endswith((".gz", ".Z", ".z")):
        fh = gzip.open(path, "rb")
    elif path.endswith((".bz", ".bz2", ".bzip2")):
        fh = bz2.open(path, "rb")
    else:
        fh = open(path, "rb")
    try:
        for line in fh:
            yield line.decode("utf8").strip().replace("\r", "")
    finally:
        fh.close()


def fasta_iter(infile: str) -> Generator[str, None, None]:
    """
    Generate seq or qual string.

    Parameters
    ----------
    infile : str
            Input fasta file.

    Yields
    ------
    str
            String of nucleotides or quality scores.
    """
    logging.info(f'Reading FASTA file "{infile}" ...')
    for line in _open_file(infile):
        line = line.strip()
        if len(line) == 0:
            continue
        if line.startswith(">"):
            continue
        yield line


def fastq_iter(infile: str, mode: str = "seq") -> Generator[str, None, None]:
    """
    Generate seq or qual string.

    Parameters
    ----------
    infile : str
            Input fastq file.
    mode : str
            Must be 'seq' or 'qual'.

    Yields
    ------
    str
            String of nucleotides or quality scores.
    """
    logging.info(f'Reading FASTQ file "{infile}" ...')
    count = 0
    s = ""
    q = ""
    for line in _open_file(infile):
        line = line.strip()
        if len(line) == 0:
            continue
        count += 1
        if count == 2:
            s = line
        if count == 4:
            q = line
            if mode == "seq":
                yield s
            elif mode == "qual":
                yield q
            count = 0


def qual2countMat(q_obj: Iterable[str], limit: int | None, step_size: int = 100000) -> dict[int, dict[int, int]]:
    """
    Generate quality-score count matrix as a nested dict.

    Parameters
    ----------
    q_obj : iterable
            Generator returned by fastq_iter
    limit : int
            Only read this many sequences.
    step_size : int
            Output progress report every step_size.

    Return
    ------
    dict mapping position -> {quality_score: count}
    """
    dat: dict[int, dict[int, int]] = collections.defaultdict(dict)
    count = 0
    for qstr in q_obj:
        count += 1
        for i, q in enumerate(qstr):
            q_score = ord(q) - 33
            if q_score not in dat[i]:
                dat[i][q_score] = 1
            else:
                dat[i][q_score] += 1
        if count % step_size == 0:
            print(f"{count} quality sequences finished\r", end=" ", file=sys.stderr)
        if limit is not None:
            if count >= limit:
                break

    logging.info(f"{count} quality sequences finished")
    return dict(dat)


def seq2countMat(
    s_obj: Iterable[str], limit: int | None, step_size: int = 100000, exclude_N: bool = False
) -> dict[int, dict[str, int]]:
    """
    Generate nucleotide count matrix as a nested dict.

    Parameters
    ----------
    s_obj : iterable
            Generator returned by fastq_iter or fasta_iter.
    limit : int
            Only read this many sequences.
    step_size : int
            Output progress report when step_size sequences have been processed.
    exclude_N : bool
            If True, sequences containing "N" will be skipped.

    Return
    ------
    dict mapping position -> {base: count}
    """
    mat: dict[int, dict[str, int]] = collections.defaultdict(dict)
    count = 0
    for s in s_obj:
        count += 1
        if (exclude_N is True) and "N" in s:
            continue
        for indx, base in enumerate(s):
            if base not in mat[indx]:
                mat[indx][base] = 1
            else:
                mat[indx][base] += 1
        if count % step_size == 0:
            print(f"{count} sequences finished\r", end=" ", file=sys.stderr)
        if limit is not None:
            if count >= limit:
                break
    logging.info(f"{count} sequences finished")
    return dict(mat)


def write_matrix_csv(
    mat: dict[int, dict[Any, int]],
    outfile: str,
    index_label: str = "Index",
    transpose: bool = False,
    sort_index_descending: bool = False,
    normalize: bool = False,
) -> None:
    """Write a nested dict matrix to CSV.

    Parameters
    ----------
    mat : dict
        Nested dict {outer_key: {inner_key: value}}.
    outfile : str
        Output CSV file path.
    index_label : str
        Label for the index column.
    transpose : bool
        If True, inner keys become rows and outer keys become columns.
    sort_index_descending : bool
        If True, sort row keys in descending order.
    normalize : bool
        If True, normalize each column to sum to 1.
    """
    if transpose:
        # Rows = inner keys, Columns = outer keys (sorted)
        col_keys = sorted(mat.keys())
        row_keys_set: set[Any] = set()
        for d in mat.values():
            row_keys_set.update(d.keys())
        row_keys = sorted(row_keys_set, reverse=sort_index_descending)

        # Build column sums for normalization
        col_sums: dict[int, float] = {}
        if normalize:
            for ck in col_keys:
                col_sums[ck] = sum(mat[ck].get(rk, 0) for rk in row_keys)

        with open(outfile, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([index_label] + [str(c) for c in col_keys])
            for rk in row_keys:
                if normalize:
                    row = [mat[ck].get(rk, 0) / col_sums[ck] if col_sums[ck] else 0 for ck in col_keys]
                else:
                    row = [mat[ck].get(rk, 0) for ck in col_keys]
                writer.writerow([rk] + row)
    else:
        # Rows = outer keys, Columns = inner keys (sorted)
        row_keys_2 = sorted(mat.keys(), reverse=sort_index_descending)
        col_keys_set: set[Any] = set()
        for d in mat.values():
            col_keys_set.update(d.keys())
        col_keys_2 = sorted(col_keys_set)

        with open(outfile, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([index_label] + [str(c) for c in col_keys_2])
            for rk in row_keys_2:
                row = [mat[rk].get(ck, 0) for ck in col_keys_2]
                writer.writerow([rk] + row)


def make_logo(
    mat: dict[int, dict[str, int]],
    outfile: str,
    exclude_N: bool = False,
    font_name: str = "sans",
    stack_order: str = "big_on_top",
    flip_below: bool = True,
    shade_below: float = 0.0,
    fade_below: float = 0.0,
    highlight_start: int = 0,
    highlight_end: int = 15,
    oformat: str = "pdf",
) -> None:
    """
    Call logomaker (https://github.com/jbkinney/logomaker) to make nucleotide logo.

    Parameters
    ----------
    mat : dict
            Nested dict {position: {base: count}}.
    exclude_N : bool
            If True, sequences containing "N" will be skipped.
    font_name : str
            The character font to use when rendering the logo.
            run logomaker.list_font_names() to get a list of
            valid font names.
    stack_order : str
            Must be 'big_on_top', 'small_on_top', or 'fixed'.
            'big_on_top' : nucleotide with the highest frequency will be on top;
            'small_on_top' : nucleotide with the lowest frequency will be on top;
            'fixed' : nucleotides from top to bottom are in the same order as characters appear in the data frame.
    flip_below : bool
            If True, characters below the x-axis (which correspond to negative
            values in the matrix) will be flipped upside down
    shade_below : float
            The amount of shading to use for characters drawn below the x-axis.
            Must be in [0.0, 1.0]
    fade_below : float
            The amount of fading to use for characters drawn below the x-axis.
            Must be in [0.0, 1.0]
    highlight_start : int
            Highlight logo from this position. Must be within [0, len(logo)-1].
    highlight_end : int
            Highlight logo to this position. Must be within [0, len(logo)-1].
    """
    import pandas as pd

    # Convert dict to DataFrame for logomaker (rows=positions, columns=bases)
    df = pd.DataFrame.from_dict(mat, orient="index").fillna(0)
    df.sort_index(inplace=True)

    logging.info("Making logo ...")
    logging.debug(f"Font name is: {font_name}")
    logging.debug(f"Stack order is: {stack_order}")
    logging.debug(f"Flip below flag: {flip_below}")
    logging.debug(f"Shade below score: {shade_below:f}")
    logging.debug(f"Fade below score: {fade_below:f}")
    if exclude_N:
        logging.info("'N' will be excluded.")
        color = {"A": "green", "C": "blue", "G": "orange", "T": "red"}
    else:
        logging.info("'N' will be kept.")
        color = {"A": "green", "C": "blue", "G": "orange", "T": "red", "N": "grey"}
    if (highlight_start is not None) and (highlight_end is not None):
        if highlight_start < 0 or highlight_end < 0 or highlight_start > highlight_end:
            logging.error("Incorrect highlight positions")
            sys.exit(1)
    logging.info(f'Mean-centered logo saved to "{outfile}.logo_mean_centered.{oformat}".')
    logo = logomaker.Logo(
        df,
        center_values=True,
        color_scheme=color,
        font_name=font_name,
        stack_order=stack_order,
        flip_below=flip_below,
        shade_below=shade_below,
        fade_below=fade_below,
    )
    if isinstance(highlight_start, int) and isinstance(highlight_end, int):
        logging.info(f"Highlight logo from {highlight_start} to {highlight_end}")
        logo.highlight_position_range(pmin=highlight_start, pmax=highlight_end)
    plt.savefig(outfile + f".logo.mean_centered.{oformat.lower()}")

    logging.info(f'Logo saved to "{outfile}.logo.{oformat}".')
    logo = logomaker.Logo(
        df,
        center_values=False,
        color_scheme=color,
        font_name=font_name,
        stack_order=stack_order,
        flip_below=flip_below,
        shade_below=shade_below,
        fade_below=fade_below,
    )
    if isinstance(highlight_start, int) and isinstance(highlight_end, int):
        logging.info(f"Highlight logo from {highlight_start} to {highlight_end}")
        logo.highlight_position_range(pmin=highlight_start, pmax=highlight_end)
    plt.savefig(outfile + f".logo.{oformat.lower()}")

"""Shared utility functions for CLI scripts."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from collections.abc import Generator
from time import strftime
from typing import Any, NamedTuple

from bx.intervals import Intersecter, Interval

import rseqc


def create_parser(description: str | None = None) -> argparse.ArgumentParser:
    """Create an ArgumentParser with standard --version flag."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--version", action="version", version=rseqc.__version__)
    return parser


def add_input_bam_arg(
    parser: argparse.ArgumentParser,
    *,
    help: str = "Alignment file in BAM or SAM format.",
    dest: str = "input_file",
    long: str = "--input-file",
) -> None:
    """Add the standard -i input BAM/SAM argument."""
    parser.add_argument("-i", long, dest=dest, help=help)


def add_mapq_arg(
    parser: argparse.ArgumentParser,
    *,
    help: str = (
        'Minimum mapping quality (phred scaled) for an alignment to be called "uniquely mapped". default=%(default)s'
    ),
) -> None:
    """Add the standard -q/--mapq argument (default=30)."""
    parser.add_argument("-q", "--mapq", type=int, dest="map_qual", default=30, help=help)


def add_output_prefix_arg(
    parser: argparse.ArgumentParser,
    *,
    help: str = "Prefix of output files(s).",
) -> None:
    """Add the standard -o/--out-prefix argument."""
    parser.add_argument("-o", "--out-prefix", dest="output_prefix", help=help)


def add_refgene_arg(
    parser: argparse.ArgumentParser,
    *,
    help: str = "Reference gene model in bed fomat. [required]",
    dest: str = "ref_gene_model",
) -> None:
    """Add the standard -r/--refgene argument."""
    parser.add_argument("-r", "--refgene", dest=dest, help=help)


def validate_files_exist(*paths: str) -> None:
    """Exit with error if any of the given file paths do not exist."""
    for path in paths:
        if not os.path.exists(path):
            print("\n\n" + path + " does NOT exists" + "\n", file=sys.stderr)
            sys.exit(1)


def printlog(mesg: str) -> None:
    """Print progress message to stderr with timestamp."""
    mesg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
    print(mesg, file=sys.stderr)


def build_bitsets(entries: list[list[Any]]) -> dict[str, Intersecter]:
    """Build interval tree from list of [chrom, start, end] entries."""
    ranges: dict[str, Intersecter] = {}
    for entry in entries:
        chrom = entry[0].upper()
        st = int(entry[1])
        end = int(entry[2])
        if chrom not in ranges:
            ranges[chrom] = Intersecter()
        ranges[chrom].add_interval(Interval(st, end))
    return ranges


def load_chromsize(file: str) -> dict[str, int]:
    """Read chrom.size file (tab/space-separated: chrom, size)."""
    chromSize: dict[str, int] = {}
    with open(file, "r") as _fh:
        for line in _fh:
            if line.startswith("#"):
                continue
            if not line.strip():
                continue
            fields = line.strip().split()
            chromSize[fields[0]] = int(fields[1])
    return chromSize


def run_rscript(script_path: str) -> None:
    """Run an R script via Rscript, silently handling missing Rscript."""
    try:
        subprocess.run(["Rscript", script_path], check=False)
    except OSError:
        print("Cannot generate pdf file from " + script_path, file=sys.stderr)


class BED12Record(NamedTuple):
    """Parsed BED12 record."""

    chrom: str
    tx_start: int
    tx_end: int
    gene_name: str
    strand: str
    exon_starts: list[int]
    exon_ends: list[int]
    fields: list[str]


def iter_bed12(bedfile: str) -> Generator[BED12Record, None, None]:
    """Iterate over BED12 records, skipping comment/track/browser headers.

    Malformed lines are skipped with a warning to stderr.
    """
    with open(bedfile, "r") as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            try:
                fields = line.split()
                chrom = fields[0]
                tx_start = int(fields[1])
                tx_end = int(fields[2])
                gene_name = fields[3]
                strand = fields[5]
                exon_starts = [int(x) for x in fields[11].rstrip(",\n").split(",")]
                exon_starts = [x + tx_start for x in exon_starts]
                exon_ends = [int(x) for x in fields[10].rstrip(",\n").split(",")]
                exon_ends = [x + y for x, y in zip(exon_starts, exon_ends)]
            except (IndexError, ValueError):
                print(
                    "[NOTE:input bed must be 12-column] skipped this line: " + line,
                    end=" ",
                    file=sys.stderr,
                )
                continue
            yield BED12Record(
                chrom=chrom,
                tx_start=tx_start,
                tx_end=tx_end,
                gene_name=gene_name,
                strand=strand,
                exon_starts=exon_starts,
                exon_ends=exon_ends,
                fields=fields,
            )

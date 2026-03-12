"""Shared utility functions for CLI scripts."""

from __future__ import annotations

import subprocess
import sys
from time import strftime
from typing import Any

from bx.intervals import Intersecter, Interval


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

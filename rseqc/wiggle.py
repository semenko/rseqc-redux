# Liguo Wang
# 04/13/2011

from __future__ import annotations

import math
import re
import sys

import bx.wiggle
from bx.binned_array import BinnedArray
from bx_extras.fpconst import isNaN


class ParseWig:
    """provie methods to manipulate wiggle format file. For wiggle format see:
    http://genome.ucsc.edu/goldenPath/help/wiggle.html"""

    def __init__(self, wigFile: str):
        """read wig file, creat wig obj"""
        self.scores: dict[str, BinnedArray] = {}
        self.num_re = re.compile(r"[\d\.\-\+]+")
        with open(wigFile) as fh:
            for i, (chrom, pos, val) in enumerate(bx.wiggle.Reader(fh)):
                chrom = chrom.upper()
                if chrom not in self.scores:
                    self.scores[chrom] = BinnedArray()
                self.scores[chrom][pos] = val
                if i % 100000 == 0:
                    print("%i datapoints loaded \r" % i)
            print("total " + str(i) + " points loaded")

    def fetch_all_scores(self, chr: str, st: int, end: int) -> list[float]:
        """fetch all wiggle scores defined by st and end.  NOTE:
        1)st and end are 0-based, half-open. (st,end]
        2)points without score are indicated as "nan"
        """
        chr = chr.upper()
        return [self.scores[chr][i] for i in range(st, end)]

    def fetch_max_scores(self, chr: str, st: int, end: int) -> float:
        """fetch maximum score defined by chr, st, end
        1)st and end are 0-based, half-open. (st,end]
        """

        chr = chr.upper()
        return float(max([self.scores[chr][i] for i in range(st, end)]))

    def fetch_min_scores(self, chr: str, st: int, end: int) -> float:
        """fetch minimum score defined by chr, st, end
        1)st and end are 0-based, half-open. (st,end]
        """

        chr = chr.upper()
        return float(min([self.scores[chr][i] for i in range(st, end)]))

    def fetch_avg_scores(self, chr: str, st: int, end: int) -> float:
        """fetch average score defined by chr, st, end
        1)st and end are 0-based, half-open. (st,end]
        """

        chr = chr.upper()
        total = 0.0
        for i in range(st, end):
            val = self.scores[chr][i]
            if not math.isnan(val):
                total += val
        return total / (end - st)

    def fetch_sum_scores(self, chr: str, st: int, end: int) -> float:
        """fetch sum score defined by chr, st, end
        1)st and end are 0-based, half-open. (st,end]
        """

        chr = chr.upper()
        total = 0.0
        for i in range(st, end):
            val = self.scores[chr][i]
            if not math.isnan(val):
                total += val
        return total


class ParseWig2:
    """provie methods to manipulate wiggle format file. For wiggle format see:
    http://genome.ucsc.edu/goldenPath/help/wiggle.html. The same coordinate could occur more than
    one time in wig file, and the scores will be sumed up. Slower than ParseWig"""

    def __init__(self, wigFile: str):
        """read wig file, creat wig obj"""
        self.scores: dict[str, BinnedArray] = {}
        self.num_re = re.compile(r"[\d\.\-\+]+")
        with open(wigFile) as fh:
            for i, (chrom, pos, val) in enumerate(bx.wiggle.Reader(fh)):
                chrom = chrom.upper()
                if chrom not in self.scores:
                    self.scores[chrom] = BinnedArray()
                tmp = self.scores[chrom][pos]
                if isNaN(tmp):
                    self.scores[chrom][pos] = val
                else:
                    self.scores[chrom][pos] += val
                if i % 100000 == 0:
                    print("%i datapoints loaded \r" % i)
            print("total " + str(i) + " points loaded")

    def fetch_all_scores_by_range(self, chr: str, st: int, end: int) -> list[float]:
        '''fetch all wiggle scores defined by st and end.  NOTE:
        1)st and end are 0-based, half-open. (st,end]
        2)points without score are indicated as "nan"'''
        chr = chr.upper()
        return [self.scores[chr][i] for i in range(st, end)]

    def fetch_all_scores_by_positions(self, chr: str, lst: list[int]) -> list[float]:
        '''fetch all wiggle scores defined by st and end.  NOTE:
        2)points without score are indicated as "nan"'''
        chr = chr.upper()
        return [self.scores[chr][i] for i in lst]

    def fetch_max_scores_by_range(self, chr: str, st: int, end: int) -> float:
        """fetch maximum score defined by chr, st, end
        1)st and end are 0-based, half-open. (st,end]
        """
        chr = chr.upper()
        return float(max([self.scores[chr][i] for i in range(st, end)]))

    def fetch_max_scores_by_positions(self, chr: str, lst: list[int]) -> float:
        """fetch maximum score defined by chr, st, end"""

        chr = chr.upper()
        return float(max([self.scores[chr][i] for i in lst]))

    def fetch_min_scores_by_range(self, chr: str, st: int, end: int) -> float:
        """fetch minimum score defined by chr, st, end
        1)st and end are 0-based, half-open. (st,end]
        """
        chr = chr.upper()
        return float(min([self.scores[chr][i] for i in range(st, end)]))

    def fetch_min_scores_by_positions(self, chr: str, lst: list[int]) -> float:
        """fetch minimum score defined by chr, st, end"""
        chr = chr.upper()
        return float(min([self.scores[chr][i] for i in lst]))

    def fetch_avg_scores_by_range(self, chr: str, st: int, end: int) -> float:
        """fetch average score defined by chr, st, end
        1)st and end are 0-based, half-open. (st,end]
        """
        chr = chr.upper()
        total = 0.0
        for i in range(st, end):
            val = self.scores[chr][i]
            if not math.isnan(val):
                total += val
        return total / (end - st)

    def fetch_avg_scores_by_positions(self, chr: str, lst: list[int]) -> float:
        """fetch average score defined by chr, st, end"""
        chr = chr.upper()
        total = 0.0
        count = 0
        for i in lst:
            val = self.scores[chr][i]
            if not math.isnan(val):
                total += val
                count += 1
        return total / count

    def fetch_sum_scores_by_range(self, chr: str, st: int, end: int) -> float:
        """fetch sum score defined by chr, st, end"""
        chr = chr.upper()
        total = 0.0
        for i in range(st, end):
            val = self.scores[chr][i]
            if not math.isnan(val):
                total += val
        return total

    def fetch_sum_scores_by_positions(self, chr: str, lst: list[int]) -> float:
        """fetch sum score defined by chr, st, end"""
        chr = chr.upper()
        total = 0.0
        for i in lst:
            val = self.scores[chr][i]
            if not math.isnan(val):
                total += val
        return total

    def distriub_wig(self, bed: str, till_count: int = 100) -> None:
        """calculate coverage over bed file (only consider exon regions). The mRNA sequences in input
        bed file will be cut into 100 tills of equal size"""

        print("Reading " + bed + " ...", file=sys.stderr)
        with open(bed, "r") as _fh:
            for line in _fh:
                try:
                    if line.startswith(("#", "track", "browser")):
                        continue
                    fields = line.rstrip("\r\n").split()
                    txStart = int(fields[1])
                    exon_start = list(map(int, fields[11].rstrip(",").split(",")))
                    exon_start = [x + txStart for x in exon_start]
                    exon_end = list(map(int, fields[10].rstrip(",").split(",")))
                    exon_end = [x + y for x, y in zip(exon_start, exon_end)]
                except (IndexError, ValueError):
                    print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=" ", file=sys.stderr)
                    continue

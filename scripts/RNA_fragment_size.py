#!/usr/bin/env python
"""
calculate fragment size for each gene/transcript. For each transcript/gene, it Will report:
1) # of fragment that was used.
2) mean of fragment size
3) median of fragment size
4) stdev of fragment size
"""

from __future__ import annotations

import os
import sys
from collections.abc import Generator

import pysam
from numpy import mean, median, std

from rseqc.cli_common import add_mapq_arg, add_refgene_arg, create_parser
from rseqc.SAM import _pysam_iter


def overlap_length2(lst1: list[list[int]], lst2: list[list[int]]) -> int:
    overlap_len = 0
    for x in lst1:
        for y in lst2:
            overlap_len += max(0, min(x[-1], y[-1]) - max(x[0], y[0]) + 1)
    return overlap_len


def fragment_size(
    bedfile: str, samfile: pysam.AlignmentFile, qcut: int = 30, ncut: int = 5
) -> Generator[str, None, None]:
    """calculate the fragment size for each gene"""
    with open(bedfile, "r") as _fh:
        for line in _fh:
            exon_range = []
            if line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            chrom = fields[0]
            tx_start = int(fields[1])
            tx_end = int(fields[2])
            geneName = fields[3]
            exon_starts = [int(x) for x in fields[11].rstrip(",\n").split(",")]
            exon_starts = [x + tx_start for x in exon_starts]
            exon_ends = [int(x) for x in fields[10].rstrip(",\n").split(",")]
            exon_ends = [x + y for x, y in zip(exon_starts, exon_ends)]
            geneID = "\t".join([str(i) for i in (chrom, tx_start, tx_end, geneName)])

            for st, end in zip(exon_starts, exon_ends):
                exon_range.append([st + 1, end + 1])

            try:
                alignedReads = samfile.fetch(chrom, tx_start, tx_end)
            except (KeyError, ValueError):
                yield "\t".join([str(i) for i in (geneID, 0, 0, 0)])
                continue

            frag_sizes = []
            for aligned_read in _pysam_iter(alignedReads):
                if not aligned_read.is_paired:  # skip single sequencing
                    continue
                if aligned_read.is_read2:
                    continue
                if aligned_read.mate_is_unmapped:
                    continue
                if aligned_read.is_qcfail:
                    continue  # skip low quanlity
                if aligned_read.is_duplicate:
                    continue  # skip duplicate read
                if aligned_read.is_secondary:
                    continue  # skip non primary hit
                if aligned_read.mapq < qcut:
                    continue

                read_st = aligned_read.pos
                mate_st = aligned_read.pnext
                if read_st > mate_st:
                    (read_st, mate_st) = (mate_st, read_st)
                if read_st < tx_start or mate_st > tx_end:
                    continue
                read_len = aligned_read.qlen
                map_range = [[read_st + 1, mate_st]]
                frag_len = overlap_length2(exon_range, map_range) + read_len
                frag_sizes.append(frag_len)
            if len(frag_sizes) < ncut:
                yield "\t".join([str(i) for i in (geneID, len(frag_sizes), 0, 0, 0)])
            else:
                yield "\t".join(
                    [str(i) for i in (geneID, len(frag_sizes), mean(frag_sizes), median(frag_sizes), std(frag_sizes))]
                )


def main() -> None:
    parser = create_parser(__doc__)
    parser.add_argument("-i", "--input", dest="input_file", help="Input BAM file")
    add_refgene_arg(
        parser,
        dest="refgene_bed",
        help="Reference gene model in BED format. Must be strandard 12-column BED file. [required]",
    )
    add_mapq_arg(parser)
    parser.add_argument(
        "-n",
        "--frag-num",
        type=int,
        dest="fragment_num",
        default=3,
        help="Minimum number of fragment. default=%(default)s",
    )

    args = parser.parse_args()

    if not (args.input_file and args.refgene_bed):
        parser.print_help()
        sys.exit(1)
    if not os.path.exists(args.input_file + ".bai"):
        print("cannot find index file of input BAM file", file=sys.stderr)
        print(args.input_file + ".bai" + " does not exists", file=sys.stderr)
        sys.exit(1)

    for file in (args.input_file, args.refgene_bed):
        if not os.path.exists(file):
            print(file + " does NOT exists" + "\n", file=sys.stderr)
            sys.exit(1)

    print(
        "\t".join(
            [
                str(i)
                for i in ("chrom", "tx_start", "tx_end", "symbol", "frag_count", "frag_mean", "frag_median", "frag_std")
            ]
        )
    )
    for tmp in fragment_size(args.refgene_bed, pysam.AlignmentFile(args.input_file), args.map_qual, args.fragment_num):
        print(tmp)


if __name__ == "__main__":
    main()

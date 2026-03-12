#!/usr/bin/env python
"""Manipulate two bigwig files"""

import argparse
import sys

import numpy
import pyBigWig

from rseqc import BED, twoList


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="5.0.2")

    parser.add_argument("-i", "--bwfile1", dest="BigWig_File1", help="One BigWig file.")
    parser.add_argument(
        "-j",
        "--bwfile2",
        dest="BigWig_File2",
        help="Another BigWig file. Both BigWig files should use the same reference genome.",
    )
    parser.add_argument(
        "-a",
        "--action",
        dest="action",
        help=(
            "After pairwise align two bigwig files, perform the"
            ' follow actions (Only select one keyword): "Add" = add'
            ' signals. "Average" = average signals. "Division" ='
            " divide bigwig2 from bigwig1. Add 1 to both bigwig."
            ' "Max" = pick the signal that is larger. "Min" = pick'
            ' the signal that is smaller. "Product" = multiply'
            ' signals. "Subtract" = subtract signals in 2nd bigwig'
            " file from the corresponding ones in the 1st bigwig"
            ' file. "geometricMean" = take the geometric mean of'
            " signals."
        ),
    )
    parser.add_argument("-o", "--output", dest="output_wig", help="Output wig file")
    parser.add_argument(
        "-c",
        "--chunk",
        type=int,
        dest="chunk_size",
        default=100000,
        help=(
            "Chromosome chunk size. Each chromosome will be cut into"
            " small chunks of this size. Decrease chunk size will save"
            " more RAM. default=%(default)s (bp)"
        ),
    )
    args = parser.parse_args()

    if not (args.BigWig_File1 and args.BigWig_File2 and args.output_wig):
        parser.print_help()
        sys.exit(1)
    with open(args.output_wig, "w") as OUT:
        bw1 = pyBigWig.open(args.BigWig_File1)
        bw2 = pyBigWig.open(args.BigWig_File2)

        print("Get chromosome sizes from BigWig header ...", file=sys.stderr)
        chrom_sizes = {}
        for chr, size in bw1.chroms().items():
            chrom_sizes[chr] = size
        for chr, size in bw2.chroms().items():
            chrom_sizes[chr] = size

        for chr_name, chr_size in chrom_sizes.items():  # iterate each chrom
            print("Processing " + chr_name + " ...", file=sys.stderr)
            OUT.write("variableStep chrom=" + chr_name + "\n")
            for interval in BED.tillingBed(chrName=chr_name, chrSize=chr_size, stepSize=args.chunk_size):
                if (bw1.stats(chr_name, interval[1], interval[2])[0] is None) and (
                    bw2.stats(chr_name, interval[1], interval[2])[0] is None
                ):
                    continue
                coord = interval[1]
                try:
                    bw_signal1 = bw1.values(chr_name, interval[1], interval[2])
                except (RuntimeError, KeyError):
                    bw_signal1 = numpy.array()
                try:
                    bw_signal2 = bw2.values(chr_name, interval[1], interval[2])
                except (RuntimeError, KeyError):
                    bw_signal2 = numpy.array()
                if bw_signal1 is None and bw_signal2 is None:
                    continue
                if numpy.isnan(numpy.nansum(bw_signal1)) and numpy.isnan(numpy.nansum(bw_signal2)):
                    continue
                if len(bw_signal1) == 0 and len(bw_signal2) == 0:
                    continue
                bw_signal1 = numpy.nan_to_num(bw_signal1)
                bw_signal2 = numpy.nan_to_num(bw_signal2)

                call_back = getattr(twoList, args.action)
                for v in call_back(bw_signal1, bw_signal2):
                    coord += 1
                    if v != 0:
                        print("%d\t%.2f" % (coord, v), file=OUT)


if __name__ == "__main__":
    main()

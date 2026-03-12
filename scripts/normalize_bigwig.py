#!/usr/bin/env python
"""Normalize bigwig signal to fixed wigsum (equivelent to total reads). Output wiggle file"""

import argparse
import collections
import sys
from itertools import groupby
from operator import itemgetter

import numpy
import pyBigWig

from rseqc import BED


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="5.0.2")

    parser.add_argument("-i", "--bwfile", dest="BigWig_File", help="Input BigWig file. [required]")
    parser.add_argument("-o", "--output", dest="output_wig", help="Output wig file. [required]")
    parser.add_argument(
        "-t",
        "--wigsum",
        type=int,
        dest="total_wigsum",
        default=100000000,
        help="Specified wigsum. 100000000 equals to coverage of 1 million 100nt reads. default=%(default)s  [optional]",
    )
    parser.add_argument(
        "-r",
        "--refgene",
        dest="refgene_bed",
        help="Reference gene model in bed format. [optional]",
    )
    parser.add_argument(
        "-c",
        "--chunk",
        type=int,
        dest="chunk_size",
        default=500000,
        help=(
            "Chromosome chunk size. Each chromosome will be cut into"
            " small chunks of this size. Decrease chunk size will save"
            " more RAM. default=%(default)s (bp) [optional]"
        ),
    )
    parser.add_argument(
        "-f",
        "--format",
        dest="out_format",
        default="bgr",
        help='Output format. either "wig" or "bgr". "bgr" save disk space but make program slower. default=%(default)s',
    )
    args = parser.parse_args()

    if not (args.BigWig_File and args.output_wig):
        parser.print_help()
        sys.exit(1)

    with open(args.output_wig, "w") as OUT:
        bw = pyBigWig.open(args.BigWig_File)

        if bw.isBigWig():
            pass
        else:
            print("%s is not a bigwig file!" % args.BigWig_File, file=sys.stderr)
            sys.exit(1)

        print("Get chromosome sizes from BigWig header ...", file=sys.stderr)
        chrom_sizes = {}
        for chr, size in bw.chroms().items():
            chrom_sizes[chr] = size

        exons = []
        WIG_SUM = 0.0
        if args.refgene_bed:
            print("Extract exons from " + args.refgene_bed, file=sys.stderr)
            obj = BED.ParseBED(args.refgene_bed)
            exons = obj.getExon()
            print("Merge overlapping exons ...", file=sys.stderr)
            exons = BED.unionBed3(exons)
            print("Calculate wigsum covered by " + args.refgene_bed + " only", file=sys.stderr)
            for chrom, st, end in exons:
                if bw.stats(chrom, st, end)[0] is None:
                    continue
                bw_signal = bw.values(chrom, st, end)
                tmp = numpy.nansum(
                    bw_signal
                )  # nan will be ignored. but if all items are 'nan', the result summay is 'nan' NOT 0
                if numpy.isnan(tmp):
                    continue
                WIG_SUM += tmp
            print("Total wigsum is %.2f\n" % WIG_SUM, file=sys.stderr)
        else:
            print("Calculate wigsum from " + args.BigWig_File, file=sys.stderr)
            for chr_name, chr_size in list(chrom_sizes.items()):  # iterate each chrom
                if bw.stats(chr_name, 0, chr_size)[0] is None:
                    print("Skip " + chr_name + "!", file=sys.stderr)
                    continue

                print("Processing " + chr_name + " ...", file=sys.stderr)
                for interval in BED.tillingBed(chrName=chr_name, chrSize=chr_size, stepSize=args.chunk_size):
                    if bw.stats(interval[0], interval[1], interval[2])[0] is None:
                        continue
                    bw_signal = bw.values(interval[0], interval[1], interval[2])
                    tmp = numpy.nansum(bw_signal)
                    if numpy.isnan(tmp):
                        continue
                    WIG_SUM += tmp
            print("\nTotal wigsum is %.2f\n" % WIG_SUM, file=sys.stderr)

        try:
            weight = args.total_wigsum / WIG_SUM
        except ZeroDivisionError:
            print("Error, WIG_SUM cannot be 0", file=sys.stderr)
            sys.exit(1)

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        print("Normalizing bigwig file ...", file=sys.stderr)
        for chr_name, chr_size in list(chrom_sizes.items()):  # iterate each chrom
            if bw.stats(chr_name, 0, chr_size)[0] is None:
                print("Skip " + chr_name + "!", file=sys.stderr)
                continue

            if args.out_format.upper() == "WIG":
                print("Writing " + chr_name + " ...", file=sys.stderr)
                OUT.write("variableStep chrom=" + chr_name + "\n")
                for interval in BED.tillingBed(chrName=chr_name, chrSize=chr_size, stepSize=args.chunk_size):
                    coord = interval[1]
                    bw_signal = bw.values(chr_name, interval[1], interval[2])
                    tmp = numpy.nansum(bw_signal)
                    if numpy.isnan(tmp):
                        continue
                    bw_signal = numpy.nan_to_num(bw_signal) * weight
                    for v in bw_signal:
                        coord += 1
                        if v != 0:
                            print("%d\t%.2f" % (coord, v), file=OUT)
            elif args.out_format.upper() == "BGR":
                print("Writing " + chr_name + " ...", file=sys.stderr)
                # OUT.write('variableStep chrom='+chr_name+'\n')
                for interval in BED.tillingBed(chrName=chr_name, chrSize=chr_size, stepSize=args.chunk_size):
                    v2p = collections.defaultdict(list)  # value to position
                    range2p = {}  # coorindate range to value, bedgraph. #[start]=[len,value]
                    coord = interval[1]
                    bw_signal = bw.values(chr_name, interval[1], interval[2])
                    tmp = numpy.nansum(bw_signal)
                    if numpy.isnan(tmp):
                        continue
                    bw_signal = numpy.nan_to_num(bw_signal) * weight
                    for v in bw_signal:
                        coord += 1
                        if v != 0:
                            v2p[v].append(coord)
                    for v in v2p:
                        for k, g in groupby(enumerate(v2p[v]), lambda i_x: i_x[0] - i_x[1]):
                            for group in [list(map(itemgetter(1), g))]:
                                range2p[group[0] - 1] = [len(group), v]
                    for i in sorted(range2p):
                        print(
                            chr_name + "\t" + str(i) + "\t" + str(i + range2p[i][0]) + "\t" + str(range2p[i][1]),
                            file=OUT,
                        )
            else:
                print("unknown output format", file=sys.stderr)
                sys.exit(1)


if __name__ == "__main__":
    main()

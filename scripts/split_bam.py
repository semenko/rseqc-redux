#!/usr/bin/env python
"""-------------------------------------------------------------------------------------------------
Split bam file according to input gene list (bed).
-------------------------------------------------------------------------------------------------"""

import os
import sys

if sys.version_info[0] != 3:
    print(
        "\nYou are using python"
        + str(sys.version_info[0])
        + "."
        + str(sys.version_info[1])
        + " This verion of RSeQC needs python3!\n",
        file=sys.stderr,
    )
    sys.exit()

import argparse

import pysam
from bx.intervals import Intersecter, Interval

from qcmodule import BED


def searchit(exon_range, exon_list):
    """return 1 if find, return 0 if cannot find"""
    for chrom, st, end in exon_list:
        if chrom.upper() not in exon_range:
            return 0
        elif len(exon_range[chrom].find(st, end)) >= 1:
            return 1
    return 0


def build_bitsets(list):
    """build intevalTree from list"""
    ranges = {}
    for entry in list:
        chrom = entry[0].upper()
        st = int(entry[1])
        end = int(entry[2])
        if chrom not in ranges:
            ranges[chrom] = Intersecter()
        ranges[chrom].add_interval(Interval(st, end))
    return ranges


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="5.0.2")
    parser.add_argument(
        "-i",
        "--input-file",
        dest="input_file",
        help="Alignment file in BAM or SAM format. BAM file should be sorted and indexed.",
    )
    parser.add_argument(
        "-r",
        "--genelist",
        dest="gene_list",
        help="Gene list in bed foramt. All reads hits to exon regions (defined by this gene list) will be saved into one BAM file, the remaining reads will saved into another BAM file.",
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        dest="output_prefix",
        help='Prefix of output BAM files. "prefix.in.bam" file contains reads mapped to the gene list specified by "-r", "prefix.ex.bam" contains reads that cannot mapped to gene list. "prefix.junk.bam" contains qcfailed or unmapped reads.',
    )
    args = parser.parse_args()

    if not (args.input_file and args.gene_list):
        parser.print_help()
        sys.exit(0)
    if not os.path.exists(args.gene_list):
        print("\n\n" + args.gene_list + " does NOT exists" + "\n", file=sys.stderr)
        # parser.print_help()
        sys.exit(0)
    if not os.path.exists(args.input_file):
        print("\n\n" + args.input_file + " does NOT exists" + "\n", file=sys.stderr)
        sys.exit(0)

    # build bitset for gene list
    print("reading " + args.gene_list + " ... ", end=" ", file=sys.stderr)
    obj = BED.ParseBED(args.gene_list)
    exons = obj.getExon()
    exon_ranges = build_bitsets(exons)
    print("Done", file=sys.stderr)

    samfile = pysam.Samfile(args.input_file, "rb")
    out1 = pysam.Samfile(
        args.output_prefix + ".in.bam", "wb", template=samfile
    )  # bam file containing reads hit to exon region
    out2 = pysam.Samfile(
        args.output_prefix + ".ex.bam", "wb", template=samfile
    )  # bam file containing reads not hit to exon region
    out3 = pysam.Samfile(
        args.output_prefix + ".junk.bam", "wb", template=samfile
    )  # bam file containing reads not hit to exon region

    total_alignment = 0
    in_alignment = 0
    ex_alignment = 0
    bad_alignment = 0
    print("spliting " + args.input_file + " ...", end=" ", file=sys.stderr)
    try:
        while 1:
            aligned_read = next(samfile)
            total_alignment += 1

            if aligned_read.is_qcfail:
                bad_alignment += 1
                out3.write(aligned_read)
                continue
            if aligned_read.is_unmapped:
                bad_alignment += 1
                out3.write(aligned_read)
                continue

            chrom = samfile.getrname(aligned_read.tid)
            chrom = chrom.upper()
            read_start = aligned_read.pos
            mate_start = aligned_read.mpos

            # read_exons = bam_cigar.fetch_exon(chrom, aligned_read.pos, aligned_read.cigar)
            if aligned_read.mate_is_unmapped:  # only one end mapped
                if chrom not in exon_ranges:
                    out2.write(aligned_read)
                    ex_alignment += 1
                    continue
                else:
                    if len(exon_ranges[chrom].find(read_start, read_start + 1)) >= 1:
                        out1.write(aligned_read)
                        in_alignment += 1
                        continue
                    elif len(exon_ranges[chrom].find(read_start, read_start + 1)) == 0:
                        out2.write(aligned_read)
                        ex_alignment += 1
                        continue
            else:  # both end mapped
                if chrom not in exon_ranges:
                    out2.write(aligned_read)
                    ex_alignment += 1
                    continue
                else:
                    if (len(exon_ranges[chrom].find(read_start, read_start + 1)) >= 1) or (
                        len(exon_ranges[chrom].find(mate_start, mate_start + 1)) >= 1
                    ):
                        out1.write(aligned_read)
                        in_alignment += 1
                    else:
                        out2.write(aligned_read)
                        ex_alignment += 1

    except StopIteration:
        print("Done", file=sys.stderr)

    print("%-55s%d" % ("Total records:", total_alignment))
    print("%-55s%d" % (args.output_prefix + ".in.bam (Alignments consumed by input gene list):", in_alignment))
    print("%-55s%d" % (args.output_prefix + ".ex.bam (Alignments not consumed by input gene list):", ex_alignment))
    print("%-55s%d" % (args.output_prefix + ".junk.bam (qcfailed, unmapped reads):", bad_alignment))


if __name__ == "__main__":
    main()

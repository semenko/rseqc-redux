#!/usr/bin/env python
"""=================================================================================================
Infer RNA-seq experiment design from SAM/BAM file. This module will determine if the RNA-seq
experiment is:
1) pair-end or single-end
2) if experiment is strand-specific, how reads were stranded.
 * For pair-end RNA-seq, there are two different ways to strand reads:
  i) 1++,1--,2+-,2-+
     read1 mapped to '+' strand indicates parental gene on '+' strand
     read1 mapped to '-' strand indicates parental gene on '-' strand
     read2 mapped to '+' strand indicates parental gene on '-' strand
     read2 mapped to '-' strand indicates parental gene on '+' strand
  ii) 1+-,1-+,2++,2--
     read1 mapped to '+' strand indicates parental gene on '-' strand
     read1 mapped to '-' strand indicates parental gene on '+' strand
     read2 mapped to '+' strand indicates parental gene on '+' strand
     read2 mapped to '-' strand indicates parental gene on '-' strand
 * For single-end RNA-seq, there are two different ways to strand reads:
  i) ++,--
     read mapped to '+' strand indicates parental gene on '+' strand
     read mapped to '-' strand indicates parental gene on '-' strand
  ii) +-,-+
     read mapped to '+' strand indicates parental gene on '-' strand
     read mapped to '-' strand indicates parental gene on '+' strand

 NOTE:
        You don't need to know the RNA sequencing protocol before mapping your reads to the reference
        genome. Mapping your RNA-seq reads as if they were non-strand specific, this script can
        "guess" how RNA-seq reads were stranded.
================================================================================================="""

import argparse
import os
import sys
from time import strftime

from rseqc import SAM


def printlog(mesg):
    """print progress into stderr and log file"""
    mesg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
    LOG = open("class.log", "a")
    print(mesg, file=sys.stderr)
    print(mesg, file=LOG)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="5.0.2")
    parser.add_argument(
        "-i",
        "--input-file",
        dest="input_file",
        help="Input alignment file in SAM or BAM format",
    )
    parser.add_argument("-r", "--refgene", dest="refgene_bed", help="Reference gene model in bed fomat.")
    parser.add_argument(
        "-s",
        "--sample-size",
        type=int,
        dest="sample_size",
        default=200000,
        help="Number of reads sampled from SAM/BAM file. default=%(default)s",
    )
    parser.add_argument(
        "-q",
        "--mapq",
        type=int,
        dest="map_qual",
        default=30,
        help=(
            "Minimum mapping quality (phred scaled) for an alignment"
            ' to be considered as "uniquely mapped".'
            " default=%(default)s"
        ),
    )

    args = parser.parse_args()

    if not (args.input_file and args.refgene_bed):
        parser.print_help()
        print("\n\n" + __doc__, file=sys.stderr)
        sys.exit(0)
    for f in (args.input_file, args.refgene_bed):
        if not os.path.exists(f):
            print("\n\n" + f + " does NOT exists." + "\n", file=sys.stderr)
            sys.exit(0)
    if args.sample_size < 1000:
        print("Warn: Sample Size too small to give a accurate estimation", file=sys.stderr)
    obj = SAM.ParseBAM(args.input_file)
    (protocol, sp1, sp2, other) = obj.configure_experiment(
        refbed=args.refgene_bed, sample_size=args.sample_size, q_cut=args.map_qual
    )
    if other < 0:
        other = 0.0
    if protocol == "PairEnd":
        print("\n\nThis is PairEnd Data")
        print("Fraction of reads failed to determine: %.4f" % other)
        print('Fraction of reads explained by "1++,1--,2+-,2-+": %.4f' % sp1)
        print('Fraction of reads explained by "1+-,1-+,2++,2--": %.4f' % sp2)

    elif protocol == "SingleEnd":
        print("\n\nThis is SingleEnd Data")
        print("Fraction of reads failed to determine: %.4f" % other)
        print('Fraction of reads explained by "++,--": %.4f' % sp1)
        print('Fraction of reads explained by "+-,-+": %.4f' % sp2)

    else:
        print("Unknown Data type")
    # print mesg


if __name__ == "__main__":
    main()

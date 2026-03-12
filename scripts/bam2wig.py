#!/usr/bin/env python
"""
Convert BAM file into wig file. BAM file must be sorted and indexed using SAMtools.
Note: SAM format file is not supported.
"""

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
from time import strftime

from qcmodule import SAM


def printlog(mesg):
    """print progress into stderr and log file"""
    mesg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
    LOG = open("class.log", "a")
    print(mesg, file=sys.stderr)
    print(mesg, file=LOG)


def load_chromsize(file):
    """read chrom.size file"""
    chromSize = {}
    for line in open(file, "r"):
        if line.startswith("#"):
            continue
        if not line.strip():
            continue
        fields = line.strip().split()
        chromSize[fields[0]] = int(fields[1])
    return chromSize


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="5.0.2")
    parser.add_argument(
        "-i",
        "--input-file",
        dest="input_file",
        help="Alignment file in BAM format. BAM file must be sorted and indexed using samTools. .bam and .bai files should be placed in the same directory.",
    )
    parser.add_argument(
        "-s",
        "--chromSize",
        dest="chromSize",
        help='Chromosome size file. Tab or space separated text file with 2 columns: first column is chromosome name/ID, second column is chromosome size. Chromosome name (such as "chr1") should be consistent between this file and the BAM file.',
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        dest="output_prefix",
        help='Prefix of output wiggle files(s). One wiggle file will be generated for non strand-specific data, two wiggle files ("Prefix_Forward.wig" and "Prefix_Reverse.wig") will be generated for strand-specific RNA-seq data.',
    )
    parser.add_argument(
        "-t",
        "--wigsum",
        type=int,
        dest="total_wigsum",
        help="Specified wigsum. Eg: 1,000,000,000 equals to coverage of 10 million 100nt reads. Ignore this option to disable normalization",
    )
    parser.add_argument(
        "-u", "--skip-multi-hits", action="store_true", dest="skip_multi", help="Skip non-unique hit reads."
    )
    parser.add_argument(
        "-d",
        "--strand",
        dest="strand_rule",
        default=None,
        help="How read(s) were stranded during sequencing. For example: --strand='1++,1--,2+-,2-+' means that this is a pair-end, strand-specific RNA-seq data, and the strand rule is: read1 mapped to '+' => parental gene on '+'; read1 mapped to '-' => parental gene on '-'; read2 mapped to '+' => parental gene on '-'; read2 mapped to '-' => parental gene on '+'.  If you are not sure about the strand rule, run 'infer_experiment.py' default=%(default)s (Not a strand specific RNA-seq data).",
    )
    parser.add_argument(
        "-q",
        "--mapq",
        type=int,
        dest="map_qual",
        default=30,
        help='Minimum mapping quality to determine "uniquely mapped". default=%(default)s',
    )

    args = parser.parse_args()

    if not (args.output_prefix and args.input_file and args.chromSize and args.output_prefix):
        parser.print_help()
        sys.exit(0)
    for file in (args.input_file, args.chromSize):
        if not os.path.exists(file):
            print("\n\n" + file + " does NOT exists" + "\n", file=sys.stderr)
            sys.exit(0)
    if not os.path.exists(args.input_file + ".bai"):
        print("index file " + args.input_file + ".bai" + " does not exists", file=sys.stderr)
        sys.exit(0)

    if args.skip_multi:
        print("Skip multi-hits:True")
    else:
        print("Skip multi-hits:False")

    chromSizes = load_chromsize(args.chromSize)

    norm_factor = None
    if args.total_wigsum:
        obj = SAM.ParseBAM(args.input_file)
        wig_sum = obj.calWigSum(chrom_sizes=chromSizes, skip_multi=args.skip_multi)
        print("\n\ntotal wigsum is:" + str(wig_sum) + "\n", file=sys.stderr)
        try:
            norm_factor = args.total_wigsum / wig_sum
        except Exception:
            norm_factor = None

    obj = SAM.ParseBAM(args.input_file)
    obj.bamTowig(
        outfile=args.output_prefix,
        chrom_sizes=chromSizes,
        chrom_file=args.chromSize,
        q_cut=args.map_qual,
        skip_multi=args.skip_multi,
        strand_rule=args.strand_rule,
        WigSumFactor=norm_factor,
    )


if __name__ == "__main__":
    main()

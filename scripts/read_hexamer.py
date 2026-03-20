#!/usr/bin/env python
"""Calculate hexamer (6-mer) frequency from read sequences."""

import os
import sys

from rseqc import FrameKmer
from rseqc.cli_common import create_parser


def main() -> None:
    parser = create_parser(__doc__)
    parser.add_argument(
        "-i",
        "--input",
        dest="input_read",
        help=(
            "Read sequence in fasta or fastq format. Multiple"
            " fasta/fastq files should be separated by ','."
            " For example: read.fq,read2.fa,read3,fa"
        ),
    )
    parser.add_argument(
        "-r",
        "--refgenome",
        dest="ref_genome",
        help="Reference genome sequence in fasta format. Optional",
    )
    parser.add_argument(
        "-g",
        "--refgene",
        dest="ref_gene",
        help="Reference mRNA sequence in fasta format. Optional",
    )
    args = parser.parse_args()

    if not args.input_read:
        parser.print_help()
        sys.exit(1)

    read_table = {}
    read_file_names = []  # base name
    read_file_sum = {}  # sum of hexamer

    for read_file in args.input_read.split(","):
        if not os.path.exists(read_file):
            print(read_file, " does NOT exist!", file=sys.stderr)
            continue
        print("Calculate hexamer of " + read_file + " file ...", end=" ", file=sys.stderr)
        read_table[os.path.basename(read_file)] = FrameKmer.kmer_freq_file(
            fastafile=read_file, word_size=6, step_size=1, frame=0
        )
        read_file_names.append(os.path.basename(read_file))
        read_file_sum[os.path.basename(read_file)] = float(sum(read_table[os.path.basename(read_file)].values()))
        print("Done", file=sys.stderr)

    if args.ref_genome and os.path.exists(args.ref_genome):
        print("Calculate hexamer of " + args.ref_genome + " file ...", end=" ", file=sys.stderr)
        read_table[os.path.basename(args.ref_genome)] = FrameKmer.kmer_freq_file(
            fastafile=args.ref_genome, word_size=6, step_size=1, frame=0
        )
        read_file_names.append(os.path.basename(args.ref_genome))
        read_file_sum[os.path.basename(args.ref_genome)] = float(
            sum(read_table[os.path.basename(args.ref_genome)].values())
        )
        print("Done.", file=sys.stderr)

    if args.ref_gene and os.path.exists(args.ref_gene):
        print("Calculate hexamer of " + args.ref_gene + " file ...", end=" ", file=sys.stderr)
        read_table[os.path.basename(args.ref_gene)] = FrameKmer.kmer_freq_file(
            fastafile=args.ref_gene, word_size=6, step_size=1, frame=0
        )
        read_file_names.append(os.path.basename(args.ref_gene))
        read_file_sum[os.path.basename(args.ref_gene)] = float(
            sum(read_table[os.path.basename(args.ref_gene)].values())
        )
        print("Done.", file=sys.stderr)

    print("Hexamer" + "\t" + "\t".join(read_file_names))

    for kmer in FrameKmer.all_possible_kmer(6):
        if "N" in kmer:
            continue
        print(kmer + "\t", end=" ")
        try:
            print("\t".join([str(read_table[name][kmer] / (read_file_sum[name])) for name in read_file_names]))
        except ZeroDivisionError:
            print("\t".join([str(read_table[name][kmer] / (read_file_sum[name] + 1)) for name in read_file_names]))


if __name__ == "__main__":
    main()

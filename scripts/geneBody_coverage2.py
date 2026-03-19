#!/usr/bin/env python
"""
Calculate the RNA-seq reads coverage over gene body.
This module uses bigwig file as input.
"""

import sys
from collections import defaultdict

from numpy import nan_to_num
from pyBigWig import open as openBigWig

from rseqc import mystat
from rseqc.cli_common import (
    add_output_prefix_arg,
    add_refgene_arg,
    create_parser,
    iter_bed12,
    run_rscript,
    validate_files_exist,
)


def coverageGeneBody_bigwig(bigFile: str, refbed: str, outfile: str, gtype: str = "png") -> None:
    """Calculate reads coverage over gene body, from 5'to 3'. each gene will be equally divided
    into 100 regsions. bigFile is bigwig format file"""

    if refbed is None:
        print("You must specify a bed file representing gene model", file=sys.stderr)
        sys.exit(1)

    bw = openBigWig(bigFile, "r")
    # Get chromosomes present in the bigwig file
    chroms = bw.chroms().keys()

    print("calculating coverage over gene body ...", file=sys.stderr)
    coverage = defaultdict(int)
    flag = 0
    gene_count = 0
    for record in iter_bed12(refbed):
        chrom = record.chrom
        strand = record.strand

        # Skip chromosomes present in the bed file but not present in the bigwig file.
        # This could happen with PATCHES or Unplaced chromosomes.
        if chrom not in chroms:
            continue

        # Count gene if it was properly read
        gene_count += 1

        gene_all_base = []
        mRNA_len = 0
        flag = 0
        for st, end in zip(record.exon_starts, record.exon_ends):
            gene_all_base.extend(range(st + 1, end + 1))  # 0-based coordinates on genome
            mRNA_len = len(gene_all_base)
            if mRNA_len < 100:
                flag = 1
                break
        if flag == 1:
            continue

        # Sort coordinates according to the strand
        if strand == "-":
            gene_all_base.sort(reverse=True)
        else:
            gene_all_base.sort(reverse=False)

        # Get 100 points from each gene's coordinates
        percentile_base = []
        percentile_base = mystat.percentile_list(gene_all_base)

        for i in range(0, len(percentile_base)):
            sig = bw.values(chrom, percentile_base[i] - 1, percentile_base[i])
            coverage[i] += nan_to_num(sig[0])

        print(" \t%d genes finished\r" % gene_count, end=" ", file=sys.stderr)

    # Close bigwig file
    bw.close()
    print("\n", file=sys.stderr)

    x_coord = []
    y_coord = []
    with open(outfile + ".geneBodyCoverage.txt", "w") as handle:
        handle.write("percentile\tcount\n")
        for i in coverage:
            x_coord.append(str(i))
            y_coord.append(str(coverage[i]))
            handle.write("%i\t%i\n" % (i, coverage[i]))

    with open(outfile + ".geneBodyCoverage_plot.r", "w") as handle:
        handle.write("%s('%s')\n" % (gtype, outfile + ".geneBodyCoverage." + gtype))
        handle.write("x=1:100\n")
        handle.write("y=c(%s)\n" % ",".join(y_coord))
        handle.write(
            "plot(x, y/%s, xlab=\"percentile of gene body (5'->3')\", ylab='average wigsum', type='s')\n" % gene_count
        )
        handle.write("dev.off()\n")


def main() -> None:
    parser = create_parser(__doc__)
    parser.add_argument(
        "-i",
        "--input-file",
        dest="input_file",
        help="Coverage signal file in bigwig format",
    )
    add_refgene_arg(parser)
    add_output_prefix_arg(parser, help="Prefix of output files(s). [required]")
    parser.add_argument(
        "-t",
        "--graph-type",
        dest="graph_type",
        default="pdf",
        help='Graphic file type in "pdf", "jpeg", "bmp", "bmp", "tiff" or "png".default=%(default)s [optional]',
    )
    args = parser.parse_args()

    gt = args.graph_type.lower()
    if gt not in ("pdf", "png", "bmp", "jpeg", "tiff"):
        print("graphic file type must be 'pdf' or 'png'", end="\n", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if not (args.output_prefix and args.input_file and args.ref_gene_model):
        parser.print_help()
        sys.exit(1)

    validate_files_exist(args.ref_gene_model, args.input_file)

    coverageGeneBody_bigwig(args.input_file, args.ref_gene_model, args.output_prefix, gtype=args.graph_type)
    run_rscript(args.output_prefix + ".geneBodyCoverage_plot.r")


if __name__ == "__main__":
    main()

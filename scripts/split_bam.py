#!/usr/bin/env python
"""Split BAM file according to an input gene list (BED format)."""

import sys

import pysam

from rseqc import BED
from rseqc.cli_common import (
    _pysam_iter,
    add_input_bam_arg,
    add_output_prefix_arg,
    build_bitsets,
    create_parser,
    validate_files_exist,
)


def main() -> None:
    parser = create_parser(__doc__)
    add_input_bam_arg(parser, help="Alignment file in BAM or SAM format. BAM file should be sorted and indexed.")
    parser.add_argument(
        "-r",
        "--genelist",
        dest="gene_list",
        help=(
            "Gene list in bed format. All reads hits to exon regions"
            " (defined by this gene list) will be saved into one BAM"
            " file, the remaining reads will saved into another BAM"
            " file."
        ),
    )
    add_output_prefix_arg(
        parser,
        help=(
            'Prefix of output BAM files. "prefix.in.bam" file'
            " contains reads mapped to the gene list specified by"
            ' "-r", "prefix.ex.bam" contains reads that cannot'
            ' mapped to gene list. "prefix.junk.bam" contains'
            " qcfailed or unmapped reads."
        ),
    )
    args = parser.parse_args()

    if not (args.input_file and args.gene_list):
        parser.print_help()
        sys.exit(1)
    validate_files_exist(args.gene_list, args.input_file)

    # build bitset for gene list
    print("reading " + args.gene_list + " ... ", end=" ", file=sys.stderr)
    obj = BED.ParseBED(args.gene_list)
    exons = obj.getExon()
    exon_ranges = build_bitsets(exons)
    print("Done", file=sys.stderr)

    with (
        pysam.AlignmentFile(args.input_file, "rb") as samfile,
        pysam.AlignmentFile(args.output_prefix + ".in.bam", "wb", template=samfile) as out_in,
        pysam.AlignmentFile(args.output_prefix + ".ex.bam", "wb", template=samfile) as out_ex,
        pysam.AlignmentFile(args.output_prefix + ".junk.bam", "wb", template=samfile) as out_junk,
    ):
        total_alignment = 0
        in_alignment = 0
        ex_alignment = 0
        bad_alignment = 0
        print("spliting " + args.input_file + " ...", end=" ", file=sys.stderr)
        for aligned_read in _pysam_iter(samfile):
            total_alignment += 1

            if aligned_read.is_qcfail or aligned_read.is_unmapped:
                bad_alignment += 1
                out_junk.write(aligned_read)
                continue

            chrom = aligned_read.reference_name.upper()

            bitset = exon_ranges.get(chrom)
            if bitset is None:
                out_ex.write(aligned_read)
                ex_alignment += 1
                continue

            read_hit = len(bitset.find(aligned_read.pos, aligned_read.pos + 1)) >= 1
            mate_hit = (
                not aligned_read.mate_is_unmapped and len(bitset.find(aligned_read.mpos, aligned_read.mpos + 1)) >= 1
            )

            if read_hit or mate_hit:
                out_in.write(aligned_read)
                in_alignment += 1
            else:
                out_ex.write(aligned_read)
                ex_alignment += 1

    print("Done", file=sys.stderr)

    print(f"{'Total records:':<55}{total_alignment}")
    print(f"{args.output_prefix + '.in.bam (Alignments consumed by input gene list):':<55}{in_alignment}")
    print(f"{args.output_prefix + '.ex.bam (Alignments not consumed by input gene list):':<55}{ex_alignment}")
    print(f"{args.output_prefix + '.junk.bam (qcfailed, unmapped reads):':<55}{bad_alignment}")


if __name__ == "__main__":
    main()

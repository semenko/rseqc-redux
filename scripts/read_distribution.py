#!/usr/bin/env python
"""
Check reads distribution over exon, intron, UTR, intergenic ... etc
The following reads will be skipped:
        qc_failed
        PCR duplicate
        Unmapped
        Non-primary (or secondary)
"""

import sys

from bx.intervals import Intersecter, Interval

from rseqc import BED, SAM, bam_cigar
from rseqc.cli_common import add_input_bam_arg, add_refgene_arg, build_bitsets, create_parser, validate_files_exist
from rseqc.SAM import _pysam_iter


def cal_size(list: list[list]) -> int:
    """calcualte bed list total size"""
    size = 0
    for entry in list:
        size += entry[2] - entry[1]
    return size


def foundone(chrom: str, ranges: dict, st: int, end: int) -> int:
    found = 0
    if chrom in ranges:
        found = len(ranges[chrom].find(st, end))
    return found


def build_unified_tree(
    labeled_lists: list[tuple[str, list[list]]],
) -> dict[str, Intersecter]:
    """Build a single interval tree per chrom with labeled intervals.

    Each entry in labeled_lists is (label, bed3_list) where bed3_list is
    a list of [chrom, start, end] entries.
    """
    tree: dict[str, Intersecter] = {}
    for label, entries in labeled_lists:
        for entry in entries:
            chrom = entry[0].upper()
            st = int(entry[1])
            end = int(entry[2])
            if chrom not in tree:
                tree[chrom] = Intersecter()
            tree[chrom].add_interval(Interval(st, end, value=label))
    return tree


def process_gene_model(gene_model: str) -> tuple:
    print("processing " + gene_model + " ...", end=" ", file=sys.stderr)
    obj = BED.ParseBED(gene_model)
    utr_3 = obj.getUTR(utr=3)
    utr_5 = obj.getUTR(utr=5)
    cds_exon = obj.getCDSExon()
    intron = obj.getIntron()

    intron = BED.unionBed3(intron)
    cds_exon = BED.unionBed3(cds_exon)
    utr_5 = BED.unionBed3(utr_5)
    utr_3 = BED.unionBed3(utr_3)

    utr_5 = BED.subtractBed3(utr_5, cds_exon)
    utr_3 = BED.subtractBed3(utr_3, cds_exon)
    intron = BED.subtractBed3(intron, cds_exon)
    intron = BED.subtractBed3(intron, utr_5)
    intron = BED.subtractBed3(intron, utr_3)

    intergenic_up_1kb = obj.getIntergenic(direction="up", size=1000)
    intergenic_down_1kb = obj.getIntergenic(direction="down", size=1000)
    intergenic_up_5kb = obj.getIntergenic(direction="up", size=5000)
    intergenic_down_5kb = obj.getIntergenic(direction="down", size=5000)
    intergenic_up_10kb = obj.getIntergenic(direction="up", size=10000)
    intergenic_down_10kb = obj.getIntergenic(direction="down", size=10000)

    # merge integenic region
    intergenic_up_1kb = BED.unionBed3(intergenic_up_1kb)
    intergenic_up_5kb = BED.unionBed3(intergenic_up_5kb)
    intergenic_up_10kb = BED.unionBed3(intergenic_up_10kb)
    intergenic_down_1kb = BED.unionBed3(intergenic_down_1kb)
    intergenic_down_5kb = BED.unionBed3(intergenic_down_5kb)
    intergenic_down_10kb = BED.unionBed3(intergenic_down_10kb)

    # purify intergenic region
    intergenic_up_1kb = BED.subtractBed3(intergenic_up_1kb, cds_exon)
    intergenic_up_1kb = BED.subtractBed3(intergenic_up_1kb, utr_5)
    intergenic_up_1kb = BED.subtractBed3(intergenic_up_1kb, utr_3)
    intergenic_up_1kb = BED.subtractBed3(intergenic_up_1kb, intron)
    intergenic_down_1kb = BED.subtractBed3(intergenic_down_1kb, cds_exon)
    intergenic_down_1kb = BED.subtractBed3(intergenic_down_1kb, utr_5)
    intergenic_down_1kb = BED.subtractBed3(intergenic_down_1kb, utr_3)
    intergenic_down_1kb = BED.subtractBed3(intergenic_down_1kb, intron)

    # purify intergenic region
    intergenic_up_5kb = BED.subtractBed3(intergenic_up_5kb, cds_exon)
    intergenic_up_5kb = BED.subtractBed3(intergenic_up_5kb, utr_5)
    intergenic_up_5kb = BED.subtractBed3(intergenic_up_5kb, utr_3)
    intergenic_up_5kb = BED.subtractBed3(intergenic_up_5kb, intron)
    intergenic_down_5kb = BED.subtractBed3(intergenic_down_5kb, cds_exon)
    intergenic_down_5kb = BED.subtractBed3(intergenic_down_5kb, utr_5)
    intergenic_down_5kb = BED.subtractBed3(intergenic_down_5kb, utr_3)
    intergenic_down_5kb = BED.subtractBed3(intergenic_down_5kb, intron)

    # purify intergenic region
    intergenic_up_10kb = BED.subtractBed3(intergenic_up_10kb, cds_exon)
    intergenic_up_10kb = BED.subtractBed3(intergenic_up_10kb, utr_5)
    intergenic_up_10kb = BED.subtractBed3(intergenic_up_10kb, utr_3)
    intergenic_up_10kb = BED.subtractBed3(intergenic_up_10kb, intron)
    intergenic_down_10kb = BED.subtractBed3(intergenic_down_10kb, cds_exon)
    intergenic_down_10kb = BED.subtractBed3(intergenic_down_10kb, utr_5)
    intergenic_down_10kb = BED.subtractBed3(intergenic_down_10kb, utr_3)
    intergenic_down_10kb = BED.subtractBed3(intergenic_down_10kb, intron)

    # build intervalTree
    cds_exon_ranges = build_bitsets(cds_exon)
    utr_5_ranges = build_bitsets(utr_5)
    utr_3_ranges = build_bitsets(utr_3)
    intron_ranges = build_bitsets(intron)
    interg_ranges_up_1kb_ranges = build_bitsets(intergenic_up_1kb)
    interg_ranges_up_5kb_ranges = build_bitsets(intergenic_up_5kb)
    interg_ranges_up_10kb_ranges = build_bitsets(intergenic_up_10kb)
    interg_ranges_down_1kb_ranges = build_bitsets(intergenic_down_1kb)
    interg_ranges_down_5kb_ranges = build_bitsets(intergenic_down_5kb)
    interg_ranges_down_10kb_ranges = build_bitsets(intergenic_down_10kb)

    # build unified labeled interval tree for single-lookup classification
    unified = build_unified_tree(
        [
            ("cds_exon", cds_exon),
            ("utr_5", utr_5),
            ("utr_3", utr_3),
            ("intron", intron),
            ("intergenic_up_1kb", intergenic_up_1kb),
            ("intergenic_up_5kb", intergenic_up_5kb),
            ("intergenic_up_10kb", intergenic_up_10kb),
            ("intergenic_down_1kb", intergenic_down_1kb),
            ("intergenic_down_5kb", intergenic_down_5kb),
            ("intergenic_down_10kb", intergenic_down_10kb),
        ]
    )

    exon_size = cal_size(cds_exon)
    intron_size = cal_size(intron)
    utr3_size = cal_size(utr_3)
    utr5_size = cal_size(utr_5)
    int_up1k_size = cal_size(intergenic_up_1kb)
    int_up5k_size = cal_size(intergenic_up_5kb)
    int_up10k_size = cal_size(intergenic_up_10kb)
    int_down1k_size = cal_size(intergenic_down_1kb)
    int_down5k_size = cal_size(intergenic_down_5kb)
    int_down10k_size = cal_size(intergenic_down_10kb)

    print("Done", file=sys.stderr)
    return (
        cds_exon_ranges,
        intron_ranges,
        utr_5_ranges,
        utr_3_ranges,
        interg_ranges_up_1kb_ranges,
        interg_ranges_up_5kb_ranges,
        interg_ranges_up_10kb_ranges,
        interg_ranges_down_1kb_ranges,
        interg_ranges_down_5kb_ranges,
        interg_ranges_down_10kb_ranges,
        exon_size,
        intron_size,
        utr5_size,
        utr3_size,
        int_up1k_size,
        int_up5k_size,
        int_up10k_size,
        int_down1k_size,
        int_down5k_size,
        int_down10k_size,
        unified,
    )


def main() -> None:
    parser = create_parser(__doc__)
    add_input_bam_arg(parser)
    add_refgene_arg(parser, help="Reference gene model in bed format.")
    args = parser.parse_args()

    if not (args.input_file and args.ref_gene_model):
        parser.print_help()
        sys.exit(1)
    validate_files_exist(args.ref_gene_model, args.input_file)

    # build bitset
    (
        cds_exon_r,
        intron_r,
        utr_5_r,
        utr_3_r,
        intergenic_up_1kb_r,
        intergenic_up_5kb_r,
        intergenic_up_10kb_r,
        intergenic_down_1kb_r,
        intergenic_down_5kb_r,
        intergenic_down_10kb_r,
        cds_exon_base,
        intron_base,
        utr_5_base,
        utr_3_base,
        intergenic_up1kb_base,
        intergenic_up5kb_base,
        intergenic_up10kb_base,
        intergenic_down1kb_base,
        intergenic_down5kb_base,
        intergenic_down10kb_base,
        unified,
    ) = process_gene_model(args.ref_gene_model)

    intron_read = 0
    cds_exon_read = 0
    utr_5_read = 0
    utr_3_read = 0

    intergenic_up1kb_read = 0
    intergenic_down1kb_read = 0
    intergenic_up5kb_read = 0
    intergenic_down5kb_read = 0
    intergenic_up10kb_read = 0
    intergenic_down10kb_read = 0

    totalReads = 0
    totalFrags = 0
    unAssignFrags = 0
    obj = SAM.ParseBAM(args.input_file)

    R_qc_fail = 0
    R_duplicate = 0
    R_nonprimary = 0
    R_unmap = 0

    print("processing " + args.input_file + " ...", end=" ", file=sys.stderr)
    for aligned_read in _pysam_iter(obj.samfile):
        if aligned_read.is_qcfail:  # skip QC fail read
            R_qc_fail += 1
            continue
        if aligned_read.is_duplicate:  # skip duplicate read
            R_duplicate += 1
            continue
        if aligned_read.is_secondary:  # skip non primary hit
            R_nonprimary += 1
            continue
        if aligned_read.is_unmapped:  # skip unmap read
            R_unmap += 1
            continue
        totalReads += 1
        chrom = obj.samfile.getrname(aligned_read.tid)
        chrom = chrom.upper()
        exons = bam_cigar.fetch_exon(aligned_read.pos, aligned_read.cigar)
        totalFrags += len(exons)

        for exn in exons:
            mid = int(exn[0]) + int((int(exn[1]) - int(exn[0])) / 2)

            # Single unified tree lookup instead of up to 11 separate lookups
            if chrom in unified:
                labels = {h.value for h in unified[chrom].find(mid, mid)}
            else:
                labels = set()

            if "cds_exon" in labels:
                cds_exon_read += 1
                continue
            elif "utr_5" in labels and "utr_3" not in labels:
                utr_5_read += 1
                continue
            elif "utr_3" in labels and "utr_5" not in labels:
                utr_3_read += 1
                continue
            elif "utr_3" in labels and "utr_5" in labels:
                unAssignFrags += 1
                continue
            elif "intron" in labels:
                intron_read += 1
                continue
            elif "intergenic_up_10kb" in labels and "intergenic_down_10kb" in labels:
                unAssignFrags += 1
                continue
            elif "intergenic_up_1kb" in labels:
                intergenic_up1kb_read += 1
                intergenic_up5kb_read += 1
                intergenic_up10kb_read += 1
            elif "intergenic_up_5kb" in labels:
                intergenic_up5kb_read += 1
                intergenic_up10kb_read += 1
            elif "intergenic_up_10kb" in labels:
                intergenic_up10kb_read += 1

            elif "intergenic_down_1kb" in labels:
                intergenic_down1kb_read += 1
                intergenic_down5kb_read += 1
                intergenic_down10kb_read += 1
            elif "intergenic_down_5kb" in labels:
                intergenic_down5kb_read += 1
                intergenic_down10kb_read += 1
            elif "intergenic_down_10kb" in labels:
                intergenic_down10kb_read += 1
            else:
                unAssignFrags += 1
    print("Finished\n", file=sys.stderr)

    print("%-30s%d" % ("Total Reads", totalReads))
    print("%-30s%d" % ("Total Tags", totalFrags))
    print("%-30s%d" % ("Total Assigned Tags", totalFrags - unAssignFrags))

    print("=====================================================================")
    print("%-20s%-20s%-20s%-20s" % ("Group", "Total_bases", "Tag_count", "Tags/Kb"))
    print(
        "%-20s%-20d%-20d%-18.2f"
        % ("CDS_Exons", cds_exon_base, cds_exon_read, cds_exon_read * 1000.0 / (cds_exon_base + 1))
    )
    print("%-20s%-20d%-20d%-18.2f" % ("5'UTR_Exons", utr_5_base, utr_5_read, utr_5_read * 1000.0 / (utr_5_base + 1)))
    print("%-20s%-20d%-20d%-18.2f" % ("3'UTR_Exons", utr_3_base, utr_3_read, utr_3_read * 1000.0 / (utr_3_base + 1)))
    print("%-20s%-20d%-20d%-18.2f" % ("Introns", intron_base, intron_read, intron_read * 1000.0 / (intron_base + 1)))

    print(
        "%-20s%-20d%-20d%-18.2f"
        % (
            "TSS_up_1kb",
            intergenic_up1kb_base,
            intergenic_up1kb_read,
            intergenic_up1kb_read * 1000.0 / (intergenic_up1kb_base + 1),
        )
    )
    print(
        "%-20s%-20d%-20d%-18.2f"
        % (
            "TSS_up_5kb",
            intergenic_up5kb_base,
            intergenic_up5kb_read,
            intergenic_up5kb_read * 1000.0 / (intergenic_up5kb_base + 1),
        )
    )
    print(
        "%-20s%-20d%-20d%-18.2f"
        % (
            "TSS_up_10kb",
            intergenic_up10kb_base,
            intergenic_up10kb_read,
            intergenic_up10kb_read * 1000.0 / (intergenic_up10kb_base + 1),
        )
    )
    print(
        "%-20s%-20d%-20d%-18.2f"
        % (
            "TES_down_1kb",
            intergenic_down1kb_base,
            intergenic_down1kb_read,
            intergenic_down1kb_read * 1000.0 / (intergenic_down1kb_base + 1),
        )
    )
    print(
        "%-20s%-20d%-20d%-18.2f"
        % (
            "TES_down_5kb",
            intergenic_down5kb_base,
            intergenic_down5kb_read,
            intergenic_down5kb_read * 1000.0 / (intergenic_down5kb_base + 1),
        )
    )
    print(
        "%-20s%-20d%-20d%-18.2f"
        % (
            "TES_down_10kb",
            intergenic_down10kb_base,
            intergenic_down10kb_read,
            intergenic_down10kb_read * 1000.0 / (intergenic_down10kb_base + 1),
        )
    )
    print("=====================================================================")


if __name__ == "__main__":
    main()

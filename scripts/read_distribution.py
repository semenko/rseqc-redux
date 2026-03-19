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
from typing import NamedTuple

from bx.intervals import Intersecter, Interval

from rseqc import BED, SAM, bam_cigar
from rseqc.cli_common import add_input_bam_arg, add_refgene_arg, create_parser, validate_files_exist
from rseqc.SAM import _pysam_iter


def cal_size(list: list[list]) -> int:
    """calcualte bed list total size"""
    size = 0
    for entry in list:
        size += entry[2] - entry[1]
    return size


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


class GeneModelResult(NamedTuple):
    """Result from process_gene_model()."""

    exon_size: int
    intron_size: int
    utr5_size: int
    utr3_size: int
    int_up1k_size: int
    int_up5k_size: int
    int_up10k_size: int
    int_down1k_size: int
    int_down5k_size: int
    int_down10k_size: int
    unified: dict[str, Intersecter]


def process_gene_model(gene_model: str) -> GeneModelResult:
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

    # purify intergenic regions by subtracting genic regions
    genic_regions = [cds_exon, utr_5, utr_3, intron]
    intergenic_list = [
        intergenic_up_1kb,
        intergenic_down_1kb,
        intergenic_up_5kb,
        intergenic_down_5kb,
        intergenic_up_10kb,
        intergenic_down_10kb,
    ]
    for i, ig in enumerate(intergenic_list):
        for genic in genic_regions:
            ig = BED.subtractBed3(ig, genic)
        intergenic_list[i] = ig
    (
        intergenic_up_1kb,
        intergenic_down_1kb,
        intergenic_up_5kb,
        intergenic_down_5kb,
        intergenic_up_10kb,
        intergenic_down_10kb,
    ) = intergenic_list

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
    return GeneModelResult(
        exon_size=exon_size,
        intron_size=intron_size,
        utr5_size=utr5_size,
        utr3_size=utr3_size,
        int_up1k_size=int_up1k_size,
        int_up5k_size=int_up5k_size,
        int_up10k_size=int_up10k_size,
        int_down1k_size=int_down1k_size,
        int_down5k_size=int_down5k_size,
        int_down10k_size=int_down10k_size,
        unified=unified,
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

    # build gene model
    model = process_gene_model(args.ref_gene_model)
    unified = model.unified
    cds_exon_base = model.exon_size
    intron_base = model.intron_size
    utr_5_base = model.utr5_size
    utr_3_base = model.utr3_size
    intergenic_up1kb_base = model.int_up1k_size
    intergenic_up5kb_base = model.int_up5k_size
    intergenic_up10kb_base = model.int_up10k_size
    intergenic_down1kb_base = model.int_down1k_size
    intergenic_down5kb_base = model.int_down5k_size
    intergenic_down10kb_base = model.int_down10k_size

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

    print(f"{'Total Reads':<30}{totalReads}")
    print(f"{'Total Tags':<30}{totalFrags}")
    print(f"{'Total Assigned Tags':<30}{totalFrags - unAssignFrags}")

    print("=====================================================================")
    print(f"{'Group':<20}{'Total_bases':<20}{'Tag_count':<20}{'Tags/Kb':<20}")
    cds_tags_kb = cds_exon_read * 1000.0 / (cds_exon_base + 1)
    print(f"{'CDS_Exons':<20}{cds_exon_base:<20}{cds_exon_read:<20}{cds_tags_kb:<18.2f}")
    utr5_label = "5'UTR_Exons"
    utr3_label = "3'UTR_Exons"
    print(f"{utr5_label:<20}{utr_5_base:<20}{utr_5_read:<20}{utr_5_read * 1000.0 / (utr_5_base + 1):<18.2f}")
    print(f"{utr3_label:<20}{utr_3_base:<20}{utr_3_read:<20}{utr_3_read * 1000.0 / (utr_3_base + 1):<18.2f}")
    print(f"{'Introns':<20}{intron_base:<20}{intron_read:<20}{intron_read * 1000.0 / (intron_base + 1):<18.2f}")

    print(
        f"{'TSS_up_1kb':<20}{intergenic_up1kb_base:<20}{intergenic_up1kb_read:<20}"
        f"{intergenic_up1kb_read * 1000.0 / (intergenic_up1kb_base + 1):<18.2f}"
    )
    print(
        f"{'TSS_up_5kb':<20}{intergenic_up5kb_base:<20}{intergenic_up5kb_read:<20}"
        f"{intergenic_up5kb_read * 1000.0 / (intergenic_up5kb_base + 1):<18.2f}"
    )
    print(
        f"{'TSS_up_10kb':<20}{intergenic_up10kb_base:<20}{intergenic_up10kb_read:<20}"
        f"{intergenic_up10kb_read * 1000.0 / (intergenic_up10kb_base + 1):<18.2f}"
    )
    print(
        f"{'TES_down_1kb':<20}{intergenic_down1kb_base:<20}{intergenic_down1kb_read:<20}"
        f"{intergenic_down1kb_read * 1000.0 / (intergenic_down1kb_base + 1):<18.2f}"
    )
    print(
        f"{'TES_down_5kb':<20}{intergenic_down5kb_base:<20}{intergenic_down5kb_read:<20}"
        f"{intergenic_down5kb_read * 1000.0 / (intergenic_down5kb_base + 1):<18.2f}"
    )
    print(
        f"{'TES_down_10kb':<20}{intergenic_down10kb_base:<20}{intergenic_down10kb_read:<20}"
        f"{intergenic_down10kb_read * 1000.0 / (intergenic_down10kb_base + 1):<18.2f}"
    )
    print("=====================================================================")


if __name__ == "__main__":
    main()

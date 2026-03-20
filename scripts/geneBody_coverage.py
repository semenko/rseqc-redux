#!/usr/bin/env python
"""
Calculate the RNA-seq reads coverage over gene body.

Note:
1) Only input sorted and indexed BAM file(s). SAM format is not supported.
2) Genes/transcripts with mRNA length < 100 will be skipped (Number specified to "-l" cannot be < 100).
"""

import collections
import operator
import sys
from os.path import basename

import numpy as np
import pysam
from numpy import mean, std

from rseqc import mystat
from rseqc.cli_common import (
    add_output_prefix_arg,
    add_refgene_arg,
    create_parser,
    get_bam_files,
    iter_bed12,
    printlog,
    run_rscript,
    validate_files_exist,
)


def valid_name(s: str) -> str:
    """make sure the string 's' is valid name for R variable"""
    symbols = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_."
    digit = "0123456789"
    rid = "_".join(i for i in s.split())  # replace space(s) with '_'
    if rid[0] in digit:
        rid = "V" + rid
    tmp = ""
    for i in rid:
        if i in symbols:
            tmp = tmp + i
        else:
            tmp = tmp + "_"
    return tmp


def pearson_moment_coefficient(lst: list[float]) -> float:
    """measure skewness"""
    mid_value = lst[len(lst) // 2]
    sigma = std(lst, ddof=1)
    return mean([((i - mid_value) / sigma) ** 3 for i in lst])


def genebody_percentile(refbed: str, mRNA_len_cut: int = 100) -> dict:
    """
    return percentile points of gene body
    mRNA length < mRNA_len_cut will be skipped
    """
    if refbed is None:
        print("You must specify a bed file representing gene model\n", file=sys.stderr)
        sys.exit(1)

    g_percentiles = {}
    transcript_count = 0
    for record in iter_bed12(refbed):
        geneID = "_".join(
            [str(j) for j in (record.chrom, record.tx_start, record.tx_end, record.gene_name, record.strand)]
        )
        transcript_count += 1
        gene_all_base = []
        for st, end in zip(record.exon_starts, record.exon_ends):
            gene_all_base.extend(list(range(st + 1, end + 1)))  # 1-based coordinates on genome
        if len(gene_all_base) < mRNA_len_cut:
            continue
        g_percentiles[geneID] = (
            record.chrom,
            record.strand,
            mystat.percentile_list(gene_all_base),
        )  # get 100 points from each gene's coordinates
    printlog("Total " + str(transcript_count) + " transcripts loaded", logfile="log.txt")
    return g_percentiles


_CHUNK_SIZE = 5_000_000  # max bases per batched count_coverage() call


def genebody_coverage(bam: str, position_list: dict) -> dict:
    """
    position_list is dict returned from genebody_percentile
    position is 1-based genome coordinate

    Batches count_coverage() calls per chromosome to avoid redundant BAM iteration.
    Instead of one count_coverage() per gene (each of which internally iterates all
    overlapping reads), nearby genes share a single call.
    """
    with pysam.AlignmentFile(bam, "rb") as samfile:
        aggreagated_cvg = collections.defaultdict(int)

        # Group genes by chromosome for batched count_coverage
        by_chrom: dict[str, list[tuple[str, list[int]]]] = collections.defaultdict(list)
        for chrom, strand, positions in position_list.values():
            by_chrom[chrom].append((strand, positions))

        gene_finished = 0
        for chrom, genes in by_chrom.items():
            # Sort genes by start position so we can batch nearby ones
            genes.sort(key=lambda g: g[1][0])

            batch_start_idx = 0
            while batch_start_idx < len(genes):
                # Start a new chunk from this gene
                range_start = genes[batch_start_idx][1][0] - 1
                if range_start < 0:
                    range_start = 0
                range_end = genes[batch_start_idx][1][-1]

                # Extend chunk to include nearby genes within _CHUNK_SIZE
                batch_end_idx = batch_start_idx + 1
                while batch_end_idx < len(genes):
                    candidate_end = genes[batch_end_idx][1][-1]
                    if candidate_end - range_start > _CHUNK_SIZE:
                        break
                    range_end = max(range_end, candidate_end)
                    batch_end_idx += 1

                try:
                    a_arr, c_arr, g_arr, t_arr = samfile.count_coverage(
                        chrom, range_start, range_end, quality_threshold=0, read_callback="all"
                    )
                except (KeyError, ValueError):
                    batch_start_idx = batch_end_idx
                    continue

                total = np.array(a_arr, dtype=np.int32)
                total += np.array(c_arr, dtype=np.int32)
                total += np.array(g_arr, dtype=np.int32)
                total += np.array(t_arr, dtype=np.int32)
                del a_arr, c_arr, g_arr, t_arr

                for strand, positions in genes[batch_start_idx:batch_end_idx]:
                    tmp = []
                    for pos in positions:
                        idx = pos - 1 - range_start
                        if 0 <= idx < len(total):
                            tmp.append(int(total[idx]))
                        else:
                            tmp.append(0)

                    if strand == "-":
                        tmp = tmp[::-1]
                    for i in range(len(tmp)):
                        aggreagated_cvg[i] += tmp[i]
                    gene_finished += 1

                    if gene_finished % 100 == 0:
                        print(f"\t{gene_finished} transcripts finished\r", end=" ", file=sys.stderr)

                batch_start_idx = batch_end_idx

        return aggreagated_cvg


def Rcode_write(dataset: list, file_prefix: str, format: str = "pdf", colNum: int = 100) -> None:
    """generate R script for visualization"""
    with open(file_prefix + ".r", "w") as ROUT:
        names = []
        datas = []
        for name, data, tmp in dataset:
            names.append(name)
            datas.append(data)
            print(name + " <- c(" + ",".join([str(i) for i in data]) + ")", file=ROUT)

        tick_pos = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        tick_lab = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

        fmt_lower = format.lower()

        # do not generate heatmap if only 1 sample
        if len(names) >= 3:
            print(
                "data_matrix" + " <- matrix(c(" + ",".join(names) + "), byrow=T, " + "ncol=" + str(colNum) + ")",
                file=ROUT,
            )
            print("rowLabel <- c(" + ",".join(['"' + i + '"' for i in names]) + ")", file=ROUT)
            print("\n", file=ROUT)
            print(f'{fmt_lower}("{file_prefix}.heatMap.{fmt_lower}")', file=ROUT)
            print("rc <- cm.colors(ncol(data_matrix))", file=ROUT)
            tick_pos_csv = ",".join([str(i) for i in tick_pos])
            tick_lab_csv = ",".join(['"' + str(i) + '"' for i in tick_lab])
            heatmap_cmd = (
                "heatmap(data_matrix"
                ', scale=c("none"),keep.dendro=F, labRow = rowLabel '
                ",Colv = NA,Rowv = NA,labCol=NA,col=cm.colors(256),"
                "margins = c(6, 8),ColSideColors = rc,"
                "cexRow=1,cexCol=1,"
                "xlab=\"Gene body percentile (5'->3')\", "
                f"add.expr=x_axis_expr <- axis(side=1,at=c({tick_pos_csv}),"
                f"labels=c({tick_lab_csv})))"
            )
            print(heatmap_cmd, file=ROUT)
            print("dev.off()", file=ROUT)

        print("\n", file=ROUT)

        print(f'{fmt_lower}("{file_prefix}.curves.{fmt_lower}")', file=ROUT)
        print(f"x=1:{colNum}", file=ROUT)
        print(
            f'icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))({len(names)})',
            file=ROUT,
        )

        if names:
            plot_cmd = (
                f"plot(x,{names[0]},type='l',xlab=\"Gene body percentile (5'->3')\","
                ' ylab="Coverage",lwd=0.8,col=icolor[1])'
            )

            if len(names) == 1:
                print(plot_cmd, file=ROUT)

            elif len(names) >= 2 and len(names) <= 6:
                print(plot_cmd, file=ROUT)
                for i in range(1, len(names)):
                    print(f"lines(x,{names[i]},type='l',col=icolor[{i + 1}])", file=ROUT)
                legend_names = ",".join(["'" + str(n) + "'" for n in names])
                print(
                    f"legend(0,1,fill=icolor[{1}:{len(names)}], legend=c({legend_names}))",
                    file=ROUT,
                )

            elif len(names) > 6:
                print("layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,2), 4, 4, byrow = TRUE))", file=ROUT)
                print(plot_cmd, file=ROUT)
                for i in range(1, len(names)):
                    print(f"lines(x,{names[i]},type='l',col=icolor[{i + 1}])", file=ROUT)
                print("par(mar=c(1,0,2,1))", file=ROUT)
                print("plot.new()", file=ROUT)
                legend_names = ",".join(["'" + str(n) + "'" for n in names])
                print(
                    f"legend(0,1,fill=icolor[{1}:{len(names)}],legend=c({legend_names}))",
                    file=ROUT,
                )

        print("dev.off()", file=ROUT)


def main() -> None:
    parser = create_parser(__doc__)
    parser.add_argument(
        "-i",
        "--input",
        dest="input_files",
        help=(
            'Input file(s) in BAM format. "-i" takes these input:'
            " 1) a single BAM file."
            ' 2) "," separated BAM files.'
            " 3) directory containing one or more bam files."
            " 4) plain text file containing the path of one or more"
            " bam file (Each row is a BAM file path). All BAM files"
            " should be sorted and indexed using samtools."
        ),
    )
    add_refgene_arg(parser)
    parser.add_argument(
        "-l",
        "--minimum_length",
        type=int,
        default=100,
        dest="min_mRNA_length",
        help='Minimum mRNA length (bp). mRNA smaller than "min_mRNA_length" will be skipped. default=%(default)s',
    )
    parser.add_argument(
        "-f",
        "--format",
        dest="output_format",
        default="pdf",
        help="Output file format, 'pdf', 'png' or 'jpeg'. default=%(default)s",
    )
    add_output_prefix_arg(parser, help="Prefix of output files(s). [required]")
    args = parser.parse_args()

    if not (args.output_prefix and args.input_files and args.ref_gene_model):
        parser.print_help()
        sys.exit(1)

    validate_files_exist(args.ref_gene_model)
    if args.min_mRNA_length < 100:
        print('The number specified to "-l" cannot be smaller than 100.' + "\n", file=sys.stderr)
        sys.exit(1)

    with open(args.output_prefix + ".geneBodyCoverage.txt", "w") as OUT1:
        print("Percentile\t" + "\t".join([str(i) for i in range(1, 101)]), file=OUT1)

        printlog("Read BED file (reference gene model) ...", logfile="log.txt")
        gene_percentiles = genebody_percentile(refbed=args.ref_gene_model, mRNA_len_cut=args.min_mRNA_length)

        printlog("Get BAM file(s) ...", logfile="log.txt")
        bamfiles = get_bam_files(args.input_files)
        for f in bamfiles:
            print("\t" + f, file=sys.stderr)

        file_container = []
        for bamfile in bamfiles:
            printlog("Processing " + basename(bamfile) + " ...", logfile="log.txt")
            cvg = genebody_coverage(bamfile, gene_percentiles)
            if len(cvg) == 0:
                print("\nCannot get coverage signal from " + basename(bamfile) + " ! Skip", file=sys.stderr)
                continue
            tmp = valid_name(basename(bamfile).replace(".bam", ""))  # scrutinize R identifer
            if file_container.count(tmp) == 0:
                print(tmp + "\t" + "\t".join([str(cvg[k]) for k in sorted(cvg)]), file=OUT1)
            else:
                print(
                    tmp + "." + str(file_container.count(tmp)) + "\t" + "\t".join([str(cvg[k]) for k in sorted(cvg)]),
                    file=OUT1,
                )
            file_container.append(tmp)

    dataset = []
    with open(args.output_prefix + ".geneBodyCoverage.txt", "r") as _fh:
        for line in _fh:
            line = line.strip()
            if line.startswith("Percentile"):
                continue
            f = line.split()
            name = f[0]
            dat = [float(i) for i in f[1:]]
            skewness = pearson_moment_coefficient(dat)
            dataset.append((name, [(i - min(dat)) / (max(dat) - min(dat)) for i in dat], skewness))
    dataset.sort(key=operator.itemgetter(2), reverse=True)

    print("\n\n", file=sys.stderr)
    print("\tSample\tSkewness", file=sys.stderr)
    for a, b, c in dataset:
        print("\t" + a + "\t" + str(c), file=sys.stderr)
    Rcode_write(dataset, args.output_prefix + ".geneBodyCoverage", format=args.output_format)

    printlog("Running R script ...", logfile="log.txt")
    run_rscript(args.output_prefix + ".geneBodyCoverage.r")


if __name__ == "__main__":
    main()

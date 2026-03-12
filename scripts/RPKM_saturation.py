#!/usr/bin/env python
"""-------------------------------------------------------------------------------------------------
For each gene, check whether the RPKM value was saturated or not. Saturation analysis is based on
re-sampling. For example, sample 5%, 10%, ... , 95%, 100% from total mapped reads, then
calculate RPKM value for each step. Strand specific sequencing protocol is supported.
-------------------------------------------------------------------------------------------------"""

import argparse
import collections
import operator
import os
import subprocess
import sys
from time import strftime

import numpy as np

from rseqc import SAM


def printlog(mesg):
    """print progress into stderr and log file"""
    mesg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
    print(mesg, file=sys.stderr)
    with open("class.log", "a") as LOG:
        print(mesg, file=LOG)


def normalize(lst):
    """normalize all numbers between 0 and 1"""
    norm_lst = []
    if max(lst) == min(lst):
        return norm_lst
    if max(lst) - min(lst) == 0:
        return norm_lst
    for i in lst:
        norm_lst.append((i - min(lst)) / (max(lst) - min(lst)))
    return norm_lst


def square_error(lst):
    """transform list into normalized squared error (squared error divided by range)"""
    SE = []
    true_rpkm = lst[-1]
    rang = max(lst) - min(lst)
    if true_rpkm == 0:
        return None
    if rang == 0:
        return None
    for i in lst:
        SE.append(abs(i - true_rpkm) / true_rpkm)
    return SE


def show_saturation(infile, outfile, rpkm_cut=0.01):

    RPKM_values = collections.defaultdict(list)
    RPKM_mean = {}
    gene_count = 0
    Quan = {"Q1": [0, 0.25], "Q2": [0.25, 0.5], "Q3": [0.5, 0.75], "Q4": [0.75, 1]}
    with open(infile) as _fh:
        for line in _fh:
            line = line.strip()
            fields = line.split()
            if fields[0].startswith("#"):
                head = [i.replace("%", "") for i in fields[6:]]
                continue
            mykey = "\t".join(fields[0:6])
            myvalue = [float(i) for i in fields[6:]]
            if max(myvalue) == 0:
                continue
            if max(myvalue) - min(myvalue) == 0:
                continue
            if np.mean(myvalue) < rpkm_cut:
                continue

            RPKM_values[mykey] = square_error(myvalue)
            RPKM_mean[mykey] = np.mean(myvalue)
            gene_count += 1
            if len(head) == 0:
                print("No head line found, exit.", file=sys.stderr)
                sys.exit(1)

    with open(outfile, "w") as ROUT:
        print("pdf('%s')" % (outfile.replace(".r", ".pdf")), file=ROUT)
        print("par(mfrow=c(2,2))", file=ROUT)
        for quantile in sorted(Quan):
            line_count = 0
            norm_RPKM = collections.defaultdict(list)
            for k, v in sorted(iter(RPKM_mean.items()), key=operator.itemgetter(1)):
                line_count += 1
                if (line_count > gene_count * Quan[quantile][0]) and (line_count <= gene_count * Quan[quantile][1]):
                    for i, j in enumerate(RPKM_values[k]):
                        norm_RPKM[head[i]].append(str(j))
            print("name=c(%s)" % (",".join(head[:-1])), file=ROUT)
            for i in head[:-1]:
                print("S%s=c(%s)" % (i, ",".join(norm_RPKM[i])), file=ROUT)
            print(
                "boxplot(%s,names=name,outline=F,ylab='Percent Relative Error',main='%s',xlab='Resampling percentage')"
                % (",".join(["100*S" + i for i in head[:-1]]), quantile),
                file=ROUT,
            )
        print("dev.off()", file=ROUT)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="5.0.2")
    parser.add_argument(
        "-i",
        "--input-file",
        dest="input_file",
        help="Alignment file in BAM or SAM format. [required]",
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        dest="output_prefix",
        help="Prefix of output files(s). [required]",
    )
    parser.add_argument(
        "-r",
        "--refgene",
        dest="refgene_bed",
        help="Reference gene model in bed fomat. [required]",
    )
    parser.add_argument(
        "-d",
        "--strand",
        dest="strand_rule",
        default=None,
        help=(
            "How read(s) were stranded during sequencing. For example:"
            " --strand='1++,1--,2+-,2-+' means that this is a pair-end,"
            " strand-specific RNA-seq, and the strand rule is: read1"
            " mapped to '+' => parental gene on '+'; read1 mapped to"
            " '-' => parental gene on '-'; read2 mapped to '+' =>"
            " parental gene on '-'; read2 mapped to '-' => parental"
            " gene on '+'. If you are not sure about the strand rule,"
            " run 'infer_experiment.py' default=%(default)s"
            " (Not a strand specific RNA-seq data)"
        ),
    )
    parser.add_argument(
        "-l",
        "--percentile-floor",
        type=int,
        dest="percentile_low_bound",
        default=5,
        help="Sampling starts from this percentile. A integer between 0 and 100. default=%(default)s",
    )
    parser.add_argument(
        "-u",
        "--percentile-ceiling",
        type=int,
        dest="percentile_up_bound",
        default=100,
        help="Sampling ends at this percentile. A integer between 0 and 100. default=%(default)s",
    )
    parser.add_argument(
        "-s",
        "--percentile-step",
        type=int,
        dest="percentile_step",
        default=5,
        help=(
            "Sampling frequency. Smaller value means more sampling"
            " times. A integer between 0 and 100."
            " default=%(default)s"
        ),
    )
    parser.add_argument(
        "-c",
        "--rpkm-cutoff",
        type=float,
        dest="rpkm_cutoff",
        default=0.01,
        help=(
            "Transcripts with RPKM smaller than this number will be ignored in visualization plot. default=%(default)s"
        ),
    )
    parser.add_argument(
        "-q",
        "--mapq",
        type=int,
        dest="map_qual",
        default=30,
        help=(
            "Minimum mapping quality (phred scaled) for an alignment"
            ' to be called "uniquely mapped". default=%(default)s'
        ),
    )

    args = parser.parse_args()

    if not (args.output_prefix and args.input_file):
        parser.print_help()
        sys.exit(1)
    if args.percentile_low_bound < 0 or args.percentile_low_bound > 100:
        print("percentile_low_bound must be larger than 0 and samller than 100", file=sys.stderr)
        sys.exit(1)
    if args.percentile_up_bound < 0 or args.percentile_up_bound > 100:
        print("percentile_up_bound must be larger than 0 and samller than 100", file=sys.stderr)
        sys.exit(1)
    if args.percentile_up_bound < args.percentile_low_bound:
        print("percentile_up_bound must be larger than percentile_low_bound", file=sys.stderr)
        sys.exit(1)
    if args.percentile_step < 0 or args.percentile_step > args.percentile_up_bound:
        print("percentile_step must be larger than 0 and samller than percentile_up_bound", file=sys.stderr)
        sys.exit(1)
    if os.path.exists(args.input_file):
        obj = SAM.ParseBAM(args.input_file)
        obj.saturation_RPKM(
            outfile=args.output_prefix,
            refbed=args.refgene_bed,
            sample_start=args.percentile_low_bound,
            sample_end=args.percentile_up_bound,
            sample_step=args.percentile_step,
            strand_rule=args.strand_rule,
            q_cut=args.map_qual,
        )
        show_saturation(
            infile=args.output_prefix + ".eRPKM.xls",
            outfile=args.output_prefix + ".saturation.r",
            rpkm_cut=args.rpkm_cutoff,
        )
        try:
            subprocess.call("Rscript " + args.output_prefix + ".saturation.r", shell=True)
        except Exception:
            pass
    else:
        print("\n\n" + args.input_file + " does NOT exists" + "\n", file=sys.stderr)
        # parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Manipulate two bigwig files"""

import sys
from typing import Any

import numpy
import pyBigWig
from numpy.typing import NDArray

from rseqc import BED
from rseqc.cli_common import create_parser


def _check_list(v1: NDArray[Any], v2: NDArray[Any]) -> None:
    """Check if the length of two arrays is the same."""
    if v1.size != v2.size:
        raise ValueError("the lenght of both arrays must be the same")


def Add(v1: NDArray[Any], v2: NDArray[Any]) -> Any:
    """Add two arrays."""
    _check_list(v1, v2)
    return v1.__add__(v2)


def Subtract(v1: NDArray[Any], v2: NDArray[Any]) -> Any:
    """Subtract v2 from v1."""
    _check_list(v1, v2)
    return v1.__sub__(v2)


def Product(v1: NDArray[Any], v2: NDArray[Any]) -> Any:
    """Return product of two arrays."""
    _check_list(v1, v2)
    return v1.__mul__(v2)


def Division(v1: NDArray[Any], v2: NDArray[Any]) -> Any:
    """Divide v1 by v2. Add 1 to both v1 and v2."""
    _check_list(v1, v2)
    return (v1 + 1) / (v2 + 1)


def Average(v1: NDArray[Any], v2: NDArray[Any]) -> Any:
    """Return arithmetic mean of two arrays."""
    _check_list(v1, v2)
    return v1.__add__(v2) / 2


def geometricMean(v1: NDArray[Any], v2: NDArray[Any]) -> Any:
    """Return geometric mean of two arrays."""
    _check_list(v1, v2)
    return (v1.__mul__(v2)) ** 0.5


def Max(v1: NDArray[Any], v2: NDArray[Any]) -> Any:
    """Pairwise comparison of two arrays. Return the max of each pair."""
    _check_list(v1, v2)
    return numpy.maximum(v1, v2)


def Min(v1: NDArray[Any], v2: NDArray[Any]) -> Any:
    """Pairwise comparison of two arrays. Return the min of each pair."""
    _check_list(v1, v2)
    return numpy.minimum(v1, v2)


_ACTIONS: dict[str, Any] = {
    "Add": Add,
    "Subtract": Subtract,
    "Product": Product,
    "Division": Division,
    "Average": Average,
    "geometricMean": geometricMean,
    "Max": Max,
    "Min": Min,
}


def main() -> None:
    parser = create_parser(__doc__)

    parser.add_argument("-i", "--bwfile1", dest="BigWig_File1", help="One BigWig file.")
    parser.add_argument(
        "-j",
        "--bwfile2",
        dest="BigWig_File2",
        help="Another BigWig file. Both BigWig files should use the same reference genome.",
    )
    parser.add_argument(
        "-a",
        "--action",
        dest="action",
        choices=list(_ACTIONS),
        help=(
            "After pairwise align two bigwig files, perform the"
            ' follow actions (Only select one keyword): "Add" = add'
            ' signals. "Average" = average signals. "Division" ='
            " divide bigwig2 from bigwig1. Add 1 to both bigwig."
            ' "Max" = pick the signal that is larger. "Min" = pick'
            ' the signal that is smaller. "Product" = multiply'
            ' signals. "Subtract" = subtract signals in 2nd bigwig'
            " file from the corresponding ones in the 1st bigwig"
            ' file. "geometricMean" = take the geometric mean of'
            " signals."
        ),
    )
    parser.add_argument("-o", "--output", dest="output_wig", help="Output wig file")
    parser.add_argument(
        "-c",
        "--chunk",
        type=int,
        dest="chunk_size",
        default=100000,
        help=(
            "Chromosome chunk size. Each chromosome will be cut into"
            " small chunks of this size. Decrease chunk size will save"
            " more RAM. default=%(default)s (bp)"
        ),
    )
    args = parser.parse_args()

    if not (args.BigWig_File1 and args.BigWig_File2 and args.output_wig):
        parser.print_help()
        sys.exit(1)
    with open(args.output_wig, "w") as OUT:
        bw1 = pyBigWig.open(args.BigWig_File1)
        bw2 = pyBigWig.open(args.BigWig_File2)
        try:
            print("Get chromosome sizes from BigWig header ...", file=sys.stderr)
            chrom_sizes = {}
            for chr, size in bw1.chroms().items():
                chrom_sizes[chr] = size
            for chr, size in bw2.chroms().items():
                chrom_sizes[chr] = size

            for chr_name, chr_size in chrom_sizes.items():  # iterate each chrom
                print("Processing " + chr_name + " ...", file=sys.stderr)
                OUT.write("variableStep chrom=" + chr_name + "\n")
                for interval in BED.tillingBed(chrName=chr_name, chrSize=chr_size, stepSize=args.chunk_size):
                    if (bw1.stats(chr_name, interval[1], interval[2])[0] is None) and (
                        bw2.stats(chr_name, interval[1], interval[2])[0] is None
                    ):
                        continue
                    coord = interval[1]
                    try:
                        bw_signal1 = bw1.values(chr_name, interval[1], interval[2])
                    except (RuntimeError, KeyError):
                        bw_signal1 = numpy.array()
                    try:
                        bw_signal2 = bw2.values(chr_name, interval[1], interval[2])
                    except (RuntimeError, KeyError):
                        bw_signal2 = numpy.array()
                    if bw_signal1 is None and bw_signal2 is None:
                        continue
                    if numpy.isnan(numpy.nansum(bw_signal1)) and numpy.isnan(numpy.nansum(bw_signal2)):
                        continue
                    if len(bw_signal1) == 0 and len(bw_signal2) == 0:
                        continue
                    bw_signal1 = numpy.nan_to_num(bw_signal1)
                    bw_signal2 = numpy.nan_to_num(bw_signal2)

                    call_back = _ACTIONS[args.action]
                    for v in call_back(bw_signal1, bw_signal2):
                        coord += 1
                        if v != 0:
                            print(f"{coord}\t{v:.2f}", file=OUT)
        finally:
            bw1.close()
            bw2.close()


if __name__ == "__main__":
    main()

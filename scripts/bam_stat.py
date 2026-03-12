#!/usr/bin/env python
"""
Summarizing mapping statistics of a BAM or SAM file.
"""

import sys

from rseqc import SAM
from rseqc.cli_common import add_input_bam_arg, add_mapq_arg, create_parser, validate_files_exist


def main() -> None:
    parser = create_parser(__doc__)
    add_input_bam_arg(parser)
    add_mapq_arg(
        parser,
        help='Minimum mapping quality (phred scaled) to determine "uniquely mapped" reads. default=%(default)s',
    )
    args = parser.parse_args()

    if not args.input_file:
        parser.print_help()
        sys.exit(1)
    validate_files_exist(args.input_file)

    obj = SAM.ParseBAM(args.input_file)
    obj.stat(q_cut=args.map_qual)


if __name__ == "__main__":
    main()

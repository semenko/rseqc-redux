"""manipulate CIGAR string represented as list of tuples (BAM file)
BAM OP  Description
0   M   alignment match
1   I   insertion to read. (relative to reference genome)
2   D   deletion from read. (relative to reference genome)
3   N   skipped region from the reference
4   S   soft clipping (clipped sequence present in SEQ)
5   H   hard clipping (clipped sequences NOT present in SEQ)
6   P   padding (silent deletion from padded reference)
7   =   sequence match
8   X   sequence mismatch

Example:
The tuple [ (0,3), (1,5), (0,2) ] refers to an alignment with 3 matches,
5 insertions and another 2 matches.

NOTE: only deal with Match, Gap, Soft Clip, Insertion, Deletion

"""


def fetch_exon(st: int, cigar: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """fetch exon regions defined by cigar. st must be zero based
    return list of tuple of (st, end)
    """
    chrom_st = st
    exon_bound = []
    for c, s in cigar:  # code and size
        if c == 0:  # match
            exon_bound.append((chrom_st, chrom_st + s))
            chrom_st += s
        elif c == 1:  # insertion to ref
            continue
        elif c == 2:  # deletion to ref
            chrom_st += s
        elif c == 3:  # gap or intron
            chrom_st += s
        elif c == 4:  # soft clipping — does NOT consume reference positions
            continue
        else:
            continue
    return exon_bound


def fetch_intron(st: int, cigar: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """fetch intron regions defined by cigar. st must be zero based
    return list of tuple of (st, end)
    """
    chrom_st = st
    intron_bound = []
    for c, s in cigar:  # code and size
        if c == 0:  # match
            chrom_st += s
        elif c == 1:  # insertion to ref
            continue
        elif c == 2:  # deletion to ref
            chrom_st += s
        elif c == 3:  # gap or intron
            intron_bound.append((chrom_st, chrom_st + s))
            chrom_st += s
        elif c == 4:  # soft clipping. We do NOT include soft clip as part of exon
            # chrom_st += s
            continue
        else:
            continue
    return intron_bound


def fetch_deletion_range(cigar: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """fetch deletion regions defined by cigar. st must be zero based
    return list of tuple of (st, end). 'st','end' is relative to the
    start of read.
    """
    del_bound = []
    st = 0
    for c, s in cigar:  # code and size
        if c == 0:  # match
            st += s
        elif c == 4:
            st += s  # soft clip
        elif c == 1:  # insertion to ref
            st += s
        elif c == 2:  # deletion to ref
            del_bound.append((st, s))  # only record the start position of deletion, and the deletion size
        elif c == 3:  # gap or intron
            continue
        else:
            continue
    return del_bound

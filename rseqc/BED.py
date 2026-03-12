# mypy: disable-error-code="assignment,var-annotated,operator,attr-defined,arg-type,misc"
import os
import sys
from collections.abc import Generator
from typing import Any

from bx.bitset_builders import binned_bitsets_from_list


class ParseBED:
    """manipulate BED (http://genome.ucsc.edu/FAQ/FAQformat.html) format file."""

    def __init__(self, bedFile: str):
        """This is constructor of ParseBED"""
        self.transtab = str.maketrans("ACGTNX", "TGCANX")
        self.f = open(bedFile, "r")
        self.fileName = os.path.basename(bedFile)
        self.ABS_fileName = bedFile

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    def close(self):
        if self.f and not self.f.closed:
            self.f.close()

    def getUTR(self, utr: int = 35) -> list[list[Any]]:
        """Extract UTR regions from input bed file (must be 12-column). output is 6-column bed format.
        When utr=35 [default], extract both 5' and 3' UTR. When utr=3, only extract 3' UTR. When utr=5,
        only extract 5' UTR"""

        ret_lst = []
        for line in self.f:
            if line.startswith("#"):
                continue
            if line.startswith("track"):
                continue
            if line.startswith("browser"):
                continue
            fields = line.rstrip("\r\n").split()
            chrom = fields[0]
            geneName = fields[3]
            strand = fields[5]
            txStart = int(fields[1])
            cdsStart = int(fields[6])
            cdsEnd = int(fields[7])
            exon_start = list(map(int, fields[11].rstrip(",").split(",")))
            exon_start = [x + txStart for x in exon_start]

            exon_end = list(map(int, fields[10].rstrip(",").split(",")))
            exon_end = [x + y for x, y in zip(exon_start, exon_end)]

            if strand == "+":
                if utr == 35 or utr == 5:
                    for st, end in zip(exon_start, exon_end):
                        if st < cdsStart:
                            utr_st = st
                            utr_end = min(end, cdsStart)
                            ret_lst.append([chrom, utr_st, utr_end, geneName, "0", strand])
                if utr == 35 or utr == 3:
                    for st, end in zip(exon_start, exon_end):
                        if end > cdsEnd:
                            utr_st = max(st, cdsEnd)
                            utr_end = end
                            ret_lst.append([chrom, utr_st, utr_end, geneName, "0", strand])
            if strand == "-":
                if utr == 35 or utr == 3:
                    for st, end in zip(exon_start, exon_end):
                        if st < cdsStart:
                            utr_st = st
                            utr_end = min(end, cdsStart)
                            ret_lst.append([chrom, utr_st, utr_end, geneName, "0", strand])
                if utr == 35 or utr == 5:
                    for st, end in zip(exon_start, exon_end):
                        if end > cdsEnd:
                            utr_st = max(st, cdsEnd)
                            utr_end = end
                            ret_lst.append([chrom, utr_st, utr_end, geneName, "0", strand])
        self.f.seek(0)
        return ret_lst

    def getExon(self) -> list[tuple[str, int, int]]:
        """Extract exon regions from input bed file (must be 12-column). output is 6-column Tab
        separated bed file, each row represents one exon"""

        ret_lst = []
        for f in self.f:
            f = f.strip().split()
            chrom = f[0]
            chrom_start = int(f[1])
            blockSizes = [int(i) for i in f[10].strip(",").split(",")]
            blockStarts = [chrom_start + int(i) for i in f[11].strip(",").split(",")]
            for base, offset in zip(blockStarts, blockSizes):
                ret_lst.append((chrom, base, base + offset))
        self.f.seek(0)
        return ret_lst

    def getTranscriptRanges(self) -> Generator[list[Any], None, None]:
        """Extract exon regions from input bed file (must be 12-column). Return ranges of
        transcript"""

        for line in self.f:
            try:
                if line.startswith("#"):
                    continue
                if line.startswith("track"):
                    continue
                if line.startswith("browser"):
                    continue
                fields = line.rstrip("\r\n").split()
                txStart = int(fields[1])
                txEnd = int(fields[2])
                chrom = fields[0]
                strand = fields[5]
                geneName = fields[3]
                yield ([chrom, txStart, txEnd, strand, geneName + ":" + chrom + ":" + str(txStart) + "-" + str(txEnd)])
            except (IndexError, ValueError):
                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=" ", file=sys.stderr)
                continue
        self.f.seek(0)

    def getCDSExon(self) -> list[list[Any]]:
        """Extract CDS exon regions from input bed file (must be 12-column)."""
        ret_lst = []
        for f in self.f:
            f = f.strip().split()
            chrom = f[0]
            chrom_start = int(f[1])
            cdsStart = int(f[6])
            cdsEnd = int(f[7])
            blockSizes = [int(i) for i in f[10].strip(",").split(",")]
            blockStarts = [chrom_start + int(i) for i in f[11].strip(",").split(",")]
            # grab cdsStart - cdsEnd
            for base, offset in zip(blockStarts, blockSizes):
                if (base + offset) < cdsStart:
                    continue
                if base > cdsEnd:
                    continue
                exon_start = max(base, cdsStart)
                exon_end = min(base + offset, cdsEnd)
                # cds_exons.append( (exon_start, exon_end) )
                ret_lst.append([chrom, exon_start, exon_end])
        self.f.seek(0)
        return ret_lst

    def getIntron(self) -> list[list[Any]]:
        """Extract Intron regions from input bed file (must be 12-column).  output is 6-column Tab
        separated bed file, each row represents one intron"""

        ret_lst = []
        for line in self.f:
            try:
                if line.startswith("#"):
                    continue
                if line.startswith("track"):
                    continue
                if line.startswith("browser"):
                    continue
                # Parse fields from gene tabls
                fields = line.split()
                chrom = fields[0]
                tx_start = int(fields[1])
                strand = fields[5].replace(" ", "_")
                if int(fields[9]) == 1:
                    continue

                exon_starts = list(map(int, fields[11].rstrip(",\n").split(",")))
                exon_starts = [x + tx_start for x in exon_starts]
                exon_ends = list(map(int, fields[10].rstrip(",\n").split(",")))
                exon_ends = [x + y for x, y in zip(exon_starts, exon_ends)]
                intron_start = exon_ends[:-1]
                intron_end = exon_starts[1:]

                if strand == "-":
                    for st, end in zip(intron_start, intron_end):
                        ret_lst.append([chrom, st, end])
                else:
                    for st, end in zip(intron_start, intron_end):
                        ret_lst.append([chrom, st, end])
            except (IndexError, ValueError):
                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=" ", file=sys.stderr)
                continue
        self.f.seek(0)
        return ret_lst

    def getIntergenic(self, direction: str = "up", size: int = 1000) -> list[list[Any]]:
        """get intergenic regions. direction=up or down or both."""

        ret_lst = []
        for line in self.f:
            if line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            chrom = fields[0]
            tx_start = int(fields[1])
            tx_end = int(fields[2])
            strand = fields[5]
            if direction == "up" or direction == "both":
                if strand == "-":
                    region_st = tx_end
                    region_end = tx_end + size
                else:
                    region_st = max(tx_start - size, 0)
                    region_end = tx_start
                ret_lst.append([chrom, region_st, region_end])
            if direction == "down" or direction == "both":
                if strand == "-":
                    region_st = max(0, tx_start - size)
                    region_end = tx_start
                else:
                    region_st = tx_end
                    region_end = tx_end + size
                ret_lst.append([chrom, region_st, region_end])
        self.f.seek(0)
        return ret_lst


def unionBed3(lst: list[list[Any]]) -> list[list[Any]]:
    """Take the union of 3 column bed files. return a new list"""
    bitsets = binned_bitsets_from_list(lst)
    ret_lst = []
    for chrom in bitsets:
        bits = bitsets[chrom]
        end = 0
        while True:
            start = bits.next_set(end)
            if start == bits.size:
                break
            end = bits.next_clear(start)
            ret_lst.append([chrom, start, end])
    bitsets = dict()
    return ret_lst


def intersectBed3(lst1: list[list[Any]], lst2: list[list[Any]]) -> list[list[Any]]:
    """Take the intersection of two bed files (3 column bed files)"""
    bits1 = binned_bitsets_from_list(lst1)
    bits2 = binned_bitsets_from_list(lst2)

    bitsets = dict()
    ret_lst = []
    for key in bits1:
        if key in bits2:
            bits1[key].iand(bits2[key])
            bitsets[key] = bits1[key]

    for chrom in bitsets:
        bits = bitsets[chrom]
        end = 0
        while True:
            start = bits.next_set(end)
            if start == bits.size:
                break
            end = bits.next_clear(start)
            ret_lst.append([chrom, start, end])
    bits1.clear()
    bits2.clear()
    bitsets.clear()
    return ret_lst


def subtractBed3(lst1: list[list[Any]], lst2: list[list[Any]]) -> list[list[Any]]:
    """subtrack lst2 from lst1"""
    bitsets1 = binned_bitsets_from_list(lst1)
    bitsets2 = binned_bitsets_from_list(lst2)

    ret_lst = []
    for chrom in bitsets1:
        bits1 = bitsets1[chrom]
        if chrom in bitsets2:
            bits2 = bitsets2[chrom]
            bits2.invert()
            bits1.iand(bits2)
        end = 0
        while True:
            start = bits1.next_set(end)
            if start == bits1.size:
                break
            end = bits1.next_clear(start)
            ret_lst.append([chrom, start, end])
    bitsets1 = dict()
    bitsets2 = dict()
    return ret_lst


def tillingBed(chrName: str, chrSize: int, stepSize: int = 10000) -> Generator[tuple[str, int, int], None, None]:
    """tilling whome genome into small sizes"""
    # tilling genome
    for start in range(0, chrSize, stepSize):
        end = start + stepSize
        if end < chrSize:
            yield (chrName, start, end)
        else:
            yield (chrName, start, chrSize)

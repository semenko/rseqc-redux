"""manipulate BAM/SAM file."""

from __future__ import annotations

import collections
import random
import re
import sys
from typing import Any

import numpy as np
import pyBigWig
import pysam
from bx.bitset import BinnedBitSet
from bx.bitset_builders import binned_bitsets_from_list
from bx.intervals import Intersecter, Interval

from rseqc import BED, bam_cigar
from rseqc.cli_common import _pysam_iter  # noqa: F401 — re-exported for backward compat

# Pre-compiled regex for parsing MD tags in mismatchProfile()
_MD_PAT = re.compile(r"(\d+)([A-Z]+)")

# Lookup table mapping ASCII byte values to NVC column indices (A=0, C=1, G=2, T=3, N=4, X=5)
_NVC_ASCII_MAP = np.full(256, 5, dtype=np.intp)
_NVC_ASCII_MAP[ord("A")] = 0
_NVC_ASCII_MAP[ord("C")] = 1
_NVC_ASCII_MAP[ord("G")] = 2
_NVC_ASCII_MAP[ord("T")] = 3
_NVC_ASCII_MAP[ord("N")] = 4


def _write_bigwig_chrom(
    bw: pyBigWig.pyBigWig,
    chr_name: str,
    coverage: np.ndarray,
    factor: float,
) -> None:
    """Write one chromosome's coverage data to an open BigWig file.

    ``coverage`` is a 1-based numpy array (index 0 unused).  Consecutive
    positions with the same value are merged into bedGraph-style intervals
    for efficient storage.  Coordinates are converted from 1-based to
    0-based half-open as required by pyBigWig.
    """
    nonzero_idx = np.nonzero(coverage)[0]
    if len(nonzero_idx) == 0:
        return

    values = coverage[nonzero_idx] * factor

    # Run-length encode: merge consecutive positions with the same value
    # Detect breakpoints where either position jumps or value changes
    pos_breaks = np.diff(nonzero_idx) != 1
    val_breaks = np.diff(values) != 0.0
    breaks = np.where(pos_breaks | val_breaks)[0] + 1

    # Build interval lists
    run_starts = np.empty(len(breaks) + 1, dtype=np.int64)
    run_starts[0] = 0
    run_starts[1:] = breaks
    run_ends_idx = np.empty(len(breaks) + 1, dtype=np.int64)
    run_ends_idx[:-1] = breaks
    run_ends_idx[-1] = len(nonzero_idx)

    # Convert 1-based positions to 0-based half-open
    starts = (nonzero_idx[run_starts] - 1).tolist()
    ends = nonzero_idx[run_ends_idx - 1].tolist()
    vals = values[run_starts].tolist()

    bw.addEntries([chr_name] * len(starts), starts, ends=ends, values=vals)


# BAM flag bits for QC filtering (SAM spec §1.4)
_QC_FAIL_FLAGS = 0x4 | 0x100 | 0x200 | 0x400  # unmapped | secondary | qcfail | duplicate


def _passes_qc(read: Any, q_cut: int) -> bool:
    """Return True if a read passes standard QC filters.

    Filters out: QC-failed, duplicate, secondary, unmapped, and low-MAPQ reads.
    Uses a single flag bitmask check instead of multiple property accesses.
    """
    if read.flag & _QC_FAIL_FLAGS:
        return False
    return bool(read.mapq >= q_cut)


def _parse_strand_rule(strand_rule: str | None) -> dict[str, str]:
    """Parse a strand rule string into a key→strand mapping.

    Returns an empty dict for non-strand-specific data.
    For paired-end (4 tokens like "1++,1--,2+-,2-+"): maps "1+" → "+", etc.
    For single-end (2 tokens like "++,--"): maps "+" → "+", etc.
    """
    if strand_rule is None:
        return {}
    parts = strand_rule.split(",")
    if len(parts) == 4:  # PairEnd, strand-specific
        return {p[0] + p[1]: p[2] for p in parts}
    if len(parts) == 2:  # SingleEnd, strand-specific
        return {p[0]: p[1] for p in parts}
    print(f"Unknown value of option: 'strand_rule' {strand_rule}", file=sys.stderr)
    sys.exit(1)


class ParseBAM:
    """This class provides fuctions to parsing/processing/transforming SAM or BAM files. The input
    file could be either SAM or BAM format file"""

    def __init__(self, inputFile: str):
        """constructor. input could be bam or sam"""
        try:
            self.samfile = pysam.AlignmentFile(inputFile, "rb")
            if len(self.samfile.header) == 0:  # type: ignore[arg-type]
                print("BAM/SAM file has no header section. Exit!", file=sys.stderr)
                sys.exit(1)
            self.bam_format = True
        except (OSError, ValueError):
            self.samfile = pysam.AlignmentFile(inputFile, "r")
            if len(self.samfile.header) == 0:  # type: ignore[arg-type]
                print("BAM/SAM file has no header section. Exit!", file=sys.stderr)
                sys.exit(1)
            self.bam_format = False

    def stat(self, q_cut: int = 30) -> None:
        """Calculate mapping statistics"""
        R_total = 0
        R_qc_fail = 0
        R_duplicate = 0
        R_nonprimary = 0
        R_unmap = 0

        R_multipleHit = 0
        R_uniqHit = 0  # all the following count should be sum to uniqHit

        R_read1 = 0
        R_read2 = 0
        R_reverse = 0
        R_forward = 0
        R_nonSplice = 0
        R_splice = 0
        R_properPair = 0
        R_pair_diff_chrom = 0

        if self.bam_format:
            print("Load BAM file ... ", end=" ", file=sys.stderr)
        else:
            print("Load SAM file ... ", end=" ", file=sys.stderr)

        for aligned_read in _pysam_iter(self.samfile):
            R_total += 1
            if aligned_read.is_qcfail:
                R_qc_fail += 1
                continue
            if aligned_read.is_duplicate:
                R_duplicate += 1
                continue
            if aligned_read.is_secondary:
                R_nonprimary += 1
                continue
            if aligned_read.is_unmapped:
                R_unmap += 1
                continue
            if aligned_read.mapq < q_cut:
                R_multipleHit += 1
                continue
            R_uniqHit += 1

            if aligned_read.is_read1:
                R_read1 += 1
            if aligned_read.is_read2:
                R_read2 += 1
            if aligned_read.is_reverse:
                R_reverse += 1
            else:
                R_forward += 1
            if any(c == 3 for c, _s in aligned_read.cigar):
                R_splice += 1
            else:
                R_nonSplice += 1
            if aligned_read.is_proper_pair:
                R_properPair += 1
                if aligned_read.tid != aligned_read.rnext:
                    R_pair_diff_chrom += 1
        print("Done", file=sys.stderr)

        print("\n#==================================================", file=sys.stdout)
        print("#All numbers are READ count", file=sys.stdout)
        print("#==================================================\n", file=sys.stdout)
        print(f"{'Total records:':<40}{R_total}", file=sys.stdout)
        print("\n", end="", file=sys.stdout)
        print(f"{'QC failed:':<40}{R_qc_fail}", file=sys.stdout)
        print(f"{'Optical/PCR duplicate:':<40}{R_duplicate}", file=sys.stdout)
        print(f"{'Non primary hits':<40}{R_nonprimary}", file=sys.stdout)
        print(f"{'Unmapped reads:':<40}{R_unmap}", file=sys.stdout)
        print(f"{'mapq < mapq_cut (non-unique):':<40}{R_multipleHit}", file=sys.stdout)
        print("\n", end="", file=sys.stdout)
        print(f"{'mapq >= mapq_cut (unique):':<40}{R_uniqHit}", file=sys.stdout)
        print(f"{'Read-1:':<40}{R_read1}", file=sys.stdout)
        print(f"{'Read-2:':<40}{R_read2}", file=sys.stdout)
        plus_label = "Reads map to '+':"
        minus_label = "Reads map to '-':"
        print(f"{plus_label:<40}{R_forward}", file=sys.stdout)
        print(f"{minus_label:<40}{R_reverse}", file=sys.stdout)
        print(f"{'Non-splice reads:':<40}{R_nonSplice}", file=sys.stdout)
        print(f"{'Splice reads:':<40}{R_splice}", file=sys.stdout)
        print(f"{'Reads mapped in proper pairs:':<40}{R_properPair}", file=sys.stdout)
        print(f"{'Proper-paired reads map to different chrom:':<40}{R_pair_diff_chrom}", file=sys.stdout)

    def configure_experiment(self, refbed: str, sample_size: int, q_cut: int = 30) -> list[str | float]:
        """Given a BAM/SAM file, this function will try to guess the RNA-seq experiment:
        1) single-end or pair-end
        2) strand_specific or not
        3) if it is strand-specific, what's the strand_ness of the protocol
        """

        # how many reads you want to sample
        count = 0
        p_strandness: dict[str, int] = collections.defaultdict(int)
        s_strandness: dict[str, int] = collections.defaultdict(int)
        # load reference gene model
        gene_ranges = {}
        print(f"Reading reference gene model {refbed} ...", end=" ", file=sys.stderr)
        with open(refbed, "r") as _fh:
            for line in _fh:
                try:
                    if line.startswith(("#", "track", "browser")):
                        continue
                    # Parse fields from gene tabls
                    fields = line.split()
                    chrom = fields[0]
                    txStart = int(fields[1])
                    txEnd = int(fields[2])
                    strand = fields[5]
                except (IndexError, ValueError):
                    print(f"[NOTE:input bed must be 12-column] skipped this line: {line}", file=sys.stderr)
                    continue
                if chrom not in gene_ranges:
                    gene_ranges[chrom] = Intersecter()
                gene_ranges[chrom].insert(txStart, txEnd, strand)
        print("Done", file=sys.stderr)

        # read SAM/BAM file
        print("Loading SAM/BAM file ... ", end=" ", file=sys.stderr)
        for aligned_read in _pysam_iter(self.samfile):
            if count >= sample_size:
                break
            if not _passes_qc(aligned_read, q_cut):
                continue

            chrom = self.samfile.getrname(aligned_read.tid)  # type: ignore[attr-defined]
            if aligned_read.is_paired:
                if aligned_read.is_read1:
                    read_id = "1"
                elif aligned_read.is_read2:
                    read_id = "2"
                if aligned_read.is_reverse:
                    map_strand = "-"
                else:
                    map_strand = "+"
                readStart = aligned_read.pos
                readEnd = readStart + aligned_read.qlen
                if chrom in gene_ranges:
                    hits = gene_ranges[chrom].find(readStart, readEnd)
                    if not hits:
                        continue
                    if len(hits) == 1:
                        strand_from_gene = hits[0]
                    else:
                        strand_from_gene = ":".join(sorted(set(hits)))
                    p_strandness[read_id + map_strand + strand_from_gene] += 1
                    count += 1
            else:
                if aligned_read.is_reverse:
                    map_strand = "-"
                else:
                    map_strand = "+"
                readStart = aligned_read.pos
                readEnd = readStart + aligned_read.qlen
                if chrom in gene_ranges:
                    hits = gene_ranges[chrom].find(readStart, readEnd)
                    if not hits:
                        continue
                    if len(hits) == 1:
                        strand_from_gene = hits[0]
                    else:
                        strand_from_gene = ":".join(sorted(set(hits)))
                    s_strandness[map_strand + strand_from_gene] += 1
                    count += 1
        print("Finished", file=sys.stderr)

        print(f"Total {count} usable reads were sampled", file=sys.stderr)
        protocol = "unknown"
        spec1 = 0.0
        spec2 = 0.0
        other = 0.0
        if len(p_strandness) > 0 and len(s_strandness) == 0:
            protocol = "PairEnd"
            # for k,v in p_strandness.items():
            spec1 = (p_strandness["1++"] + p_strandness["1--"] + p_strandness["2+-"] + p_strandness["2-+"]) / float(
                sum(p_strandness.values())
            )
            spec2 = (p_strandness["1+-"] + p_strandness["1-+"] + p_strandness["2++"] + p_strandness["2--"]) / float(
                sum(p_strandness.values())
            )
            other = 1 - spec1 - spec2

        elif len(s_strandness) > 0 and len(p_strandness) == 0:
            protocol = "SingleEnd"
            # for k,v in s_strandness.items():
            spec1 = (s_strandness["++"] + s_strandness["--"]) / float(sum(s_strandness.values()))
            spec2 = (s_strandness["+-"] + s_strandness["-+"]) / float(sum(s_strandness.values()))
            other = 1 - spec1 - spec2
        else:
            protocol = "Mixture"
            spec1 = 0
            spec2 = 0
            other = 0
        return [protocol, spec1, spec2, other]

    def bamTowig(
        self,
        outfile: str,
        chrom_sizes: dict[str, int],
        skip_multi: bool = True,
        strand_rule: str | None = None,
        WigSumFactor: float | None = None,
        q_cut: int = 30,
    ) -> None:
        """Convert BAM/SAM file to wig and BigWig files.

        ``chrom_sizes`` is a dict with chrom name as key and size as value.
        ``strand_rule`` should be determined from ``infer_experiment``, e.g.
        ``"1++,1--,2+-,2-+"``.  When ``WigSumFactor`` is provided, output
        values will be normalized (multiplied) by this factor.
        """

        strandRule = _parse_strand_rule(strand_rule)

        # Open WIG file handles
        if len(strandRule) == 0:
            FWO = open(f"{outfile}.wig", "w")
            RVO = None
        else:
            FWO = open(f"{outfile}.Forward.wig", "w")
            RVO = open(f"{outfile}.Reverse.wig", "w")

        # Open BigWig file handles
        header = list(chrom_sizes.items())
        bw_fwd = pyBigWig.open(f"{outfile}.bw" if len(strandRule) == 0 else f"{outfile}.Forward.bw", "w")
        bw_fwd.addHeader(header)
        bw_rev: pyBigWig.pyBigWig | None = None
        if len(strandRule) > 0:
            bw_rev = pyBigWig.open(f"{outfile}.Reverse.bw", "w")
            bw_rev.addHeader(header)

        try:
            read_id = ""

            for chr_name, chr_size in chrom_sizes.items():  # iterate each chrom
                try:
                    self.samfile.fetch(chr_name, 0, chr_size)
                except (KeyError, ValueError):
                    print(f"No alignments for {chr_name}. skipped", file=sys.stderr)
                    continue
                print(f"Processing {chr_name} ...", file=sys.stderr)
                if len(strandRule) == 0:
                    FWO.write(f"variableStep chrom={chr_name}\n")
                else:
                    FWO.write(f"variableStep chrom={chr_name}\n")
                    RVO.write(f"variableStep chrom={chr_name}\n")  # type: ignore[union-attr]
                # Use numpy arrays for coverage accumulation (1-based positions → index 0 unused)
                Fwig = np.zeros(chr_size + 1, dtype=np.float64)
                Rwig = np.zeros(chr_size + 1, dtype=np.float64) if len(strandRule) > 0 else None
                alignedReads = self.samfile.fetch(chr_name, 0, chr_size)
                for aligned_read in _pysam_iter(alignedReads):
                    if aligned_read.is_qcfail:
                        continue
                    if aligned_read.is_duplicate:
                        continue
                    if aligned_read.is_secondary:
                        continue
                    if aligned_read.is_unmapped:
                        continue
                    if skip_multi and aligned_read.mapq < q_cut:
                        continue
                    if aligned_read.is_paired:
                        if aligned_read.is_read1:
                            read_id = "1"
                        if aligned_read.is_read2:
                            read_id = "2"

                    if aligned_read.is_reverse:
                        map_strand = "-"
                    else:
                        map_strand = "+"

                    key = read_id + map_strand

                    for block in aligned_read.get_blocks():
                        start = block[0] + 1
                        end = block[1] + 1
                        if len(strandRule) == 0:
                            Fwig[start:end] += 1.0
                        else:
                            if strandRule[key] == "+":
                                Fwig[start:end] += 1.0
                            if strandRule[key] == "-":
                                Rwig[start:end] -= 1.0  # type: ignore[index]

                # Write non-zero positions (sparse output, matching original variableStep format)
                factor = WigSumFactor if WigSumFactor is not None else 1.0
                nonzero_f = np.nonzero(Fwig)[0]
                if len(nonzero_f) > 0:
                    positions = nonzero_f
                    values = Fwig[positions] * factor
                    lines = np.column_stack((positions, values))
                    np.savetxt(FWO, lines, fmt="%d\t%.2f")
                if len(strandRule) > 0:
                    nonzero_r = np.nonzero(Rwig)[0]  # type: ignore[arg-type]
                    if len(nonzero_r) > 0:
                        positions = nonzero_r
                        values = Rwig[positions] * factor  # type: ignore[index]
                        lines = np.column_stack((positions, values))
                        np.savetxt(RVO, lines, fmt="%d\t%.2f")  # type: ignore[arg-type]

                # Write BigWig data for this chromosome
                _write_bigwig_chrom(bw_fwd, chr_name, Fwig, factor)
                if bw_rev is not None and Rwig is not None:
                    _write_bigwig_chrom(bw_rev, chr_name, Rwig, factor)
        finally:
            FWO.close()
            if RVO is not None:
                RVO.close()
            bw_fwd.close()
            if bw_rev is not None:
                bw_rev.close()

    def calWigSum(self, chrom_sizes: dict[str, int], skip_multi: bool = True, q_cut: int = 30) -> float:
        """Calculate wigsum from BAM file.

        Uses the same mapq-based filtering as bamTowig for consistency.
        """

        print("Calcualte wigsum ... ", file=sys.stderr)
        wigsum = 0.0
        for chr_name, chr_size in chrom_sizes.items():  # iterate each chrom
            try:
                alignedReads = self.samfile.fetch(chr_name, 0, chr_size)
            except (KeyError, ValueError):
                print(f"No alignments for {chr_name}. skipped", file=sys.stderr)
                continue
            print(f"Processing {chr_name} ...", file=sys.stderr)

            for aligned_read in _pysam_iter(alignedReads):
                if aligned_read.is_qcfail:
                    continue
                if aligned_read.is_duplicate:
                    continue
                if aligned_read.is_secondary:
                    continue
                if aligned_read.is_unmapped:
                    continue
                if skip_multi and aligned_read.mapq < q_cut:
                    continue

                for block in aligned_read.get_blocks():
                    wigsum += block[1] - block[0]
        return wigsum

    def bam2fq(self, prefix: str, paired: bool = True) -> None:
        """Convert BAM/SAM into fastq files"""

        transtab = str.maketrans("ACGTNX", "TGCANX")

        if paired:
            with open(f"{prefix}.R1.fastq", "w") as OUT1, open(f"{prefix}.R2.fastq", "w") as OUT2:
                read1_count = 0
                read2_count = 0
                read_name = ""
                read_seq = ""
                read_qual = ""

                print("Convert BAM/SAM file into fastq format ... ", end=" ", file=sys.stderr)
                for aligned_read in _pysam_iter(self.samfile):
                    read_name = aligned_read.qname
                    read_seq = aligned_read.seq.upper()
                    read_qual = aligned_read.qual
                    if aligned_read.is_reverse:
                        read_seq = read_seq.translate(transtab)[::-1]
                        read_qual = read_qual[::-1]
                    if aligned_read.is_read1:
                        read1_count += 1
                        if not read_name.endswith("/1"):
                            print(f"@{read_name}/1", file=OUT1)
                        else:
                            print(f"@{read_name}", file=OUT1)
                        print(read_seq, file=OUT1)
                        print("+", file=OUT1)
                        print(read_qual, file=OUT1)
                    if aligned_read.is_read2:
                        read2_count += 1
                        if not read_name.endswith("/2"):
                            print(f"@{read_name}/2", file=OUT2)
                        else:
                            print(f"@{read_name}", file=OUT2)
                        print(read_seq, file=OUT2)
                        print("+", file=OUT2)
                        print(read_qual, file=OUT2)
                print("Done", file=sys.stderr)
            print(f"read_1 count: {read1_count}", file=sys.stderr)
            print(f"read_2 count: {read2_count}", file=sys.stderr)
        else:
            with open(f"{prefix}.fastq", "w") as OUT:
                read_count = 0
                read_name = ""
                read_seq = ""
                read_qual = ""

                print("Convert BAM/SAM file into fastq format ... ", end=" ", file=sys.stderr)
                for aligned_read in _pysam_iter(self.samfile):
                    read_name = aligned_read.qname
                    read_seq = aligned_read.seq.upper()
                    read_qual = aligned_read.qual
                    if aligned_read.is_reverse:
                        read_seq = read_seq.translate(transtab)[::-1]
                        read_qual = read_qual[::-1]
                    read_count += 1
                    print(f"@{read_name}", file=OUT)
                    print(read_seq, file=OUT)
                    print("+", file=OUT)
                    print(read_qual, file=OUT)
                print("Done", file=sys.stderr)
            print(f"read count: {read_count}", file=sys.stderr)

    def readsNVC(self, outfile: str, nx: bool = True, q_cut: int = 30) -> None:
        """for each read, calculate nucleotide frequency vs position"""
        outfile1 = f"{outfile}.NVC.xls"
        outfile2 = f"{outfile}.NVC_plot.r"
        with open(outfile1, "w") as FO, open(outfile2, "w") as RS:
            transtab = str.maketrans("ACGTNX", "TGCANX")
            read_len = 0
            base_freq: np.ndarray | None = None
            if self.bam_format:
                print("Read BAM file ... ", end=" ", file=sys.stderr)
            else:
                print("Read SAM file ... ", end=" ", file=sys.stderr)

            for aligned_read in _pysam_iter(self.samfile):
                if aligned_read.mapq < q_cut:
                    continue

                RNA_read = aligned_read.seq.upper()
                if aligned_read.is_reverse:
                    RNA_read = RNA_read.translate(transtab)[::-1]
                read_len = len(RNA_read)
                if base_freq is None:
                    base_freq = np.zeros((read_len, 6), dtype=np.int64)
                elif read_len > base_freq.shape[0]:
                    # Grow array for longer reads
                    new_freq = np.zeros((read_len, 6), dtype=np.int64)
                    new_freq[: base_freq.shape[0], :] = base_freq
                    base_freq = new_freq
                # Vectorized per-base counting via numpy lookup table
                seq_bytes = np.frombuffer(RNA_read.encode("ascii"), dtype=np.uint8)
                col_indices = _NVC_ASCII_MAP[seq_bytes]
                base_freq[np.arange(read_len), col_indices] += 1
            print("Done", file=sys.stderr)

            if base_freq is None:
                base_freq = np.zeros((0, 6), dtype=np.int64)

            print("generating data matrix ...", file=sys.stderr)
            print("Position\tA\tC\tG\tT\tN\tX", file=FO)
            a_count = []
            c_count = []
            g_count = []
            t_count = []
            n_count = []
            x_count = []
            for i in range(read_len):
                print(f"{i}\t", end=" ", file=FO)
                print(f"{base_freq[i, 0]}\t", end=" ", file=FO)
                a_count.append(str(base_freq[i, 0]))
                print(f"{base_freq[i, 1]}\t", end=" ", file=FO)
                c_count.append(str(base_freq[i, 1]))
                print(f"{base_freq[i, 2]}\t", end=" ", file=FO)
                g_count.append(str(base_freq[i, 2]))
                print(f"{base_freq[i, 3]}\t", end=" ", file=FO)
                t_count.append(str(base_freq[i, 3]))
                print(f"{base_freq[i, 4]}\t", end=" ", file=FO)
                n_count.append(str(base_freq[i, 4]))
                print(f"{base_freq[i, 5]}\t", file=FO)
                x_count.append(str(base_freq[i, 5]))

            # generating R scripts
            print("generating R script  ...", file=sys.stderr)
            pos_str = ",".join([str(i) for i in range(read_len)])
            print(f"position=c({pos_str})", file=RS)
            print(f"A_count=c({','.join(a_count)})", file=RS)
            print(f"C_count=c({','.join(c_count)})", file=RS)
            print(f"G_count=c({','.join(g_count)})", file=RS)
            print(f"T_count=c({','.join(t_count)})", file=RS)
            print(f"N_count=c({','.join(n_count)})", file=RS)
            print(f"X_count=c({','.join(x_count)})", file=RS)

            if nx:
                print("total= A_count + C_count + G_count + T_count + N_count + X_count", file=RS)
                print(
                    "ym=max(A_count/total,C_count/total,G_count/total,"
                    "T_count/total,N_count/total,X_count/total) + 0.05",
                    file=RS,
                )
                print(
                    "yn=min(A_count/total,C_count/total,G_count/total,T_count/total,N_count/total,X_count/total)",
                    file=RS,
                )

                print(f'pdf("{outfile}.NVC_plot.pdf")', file=RS)
                print(
                    'plot(position,A_count/total,type="o",pch=20,'
                    'ylim=c(yn,ym),col="dark green",'
                    'xlab="Position of Read",ylab="Nucleotide Frequency")',
                    file=RS,
                )
                print('lines(position,T_count/total,type="o",pch=20,col="red")', file=RS)
                print('lines(position,G_count/total,type="o",pch=20,col="blue")', file=RS)
                print('lines(position,C_count/total,type="o",pch=20,col="cyan")', file=RS)
                print('lines(position,N_count/total,type="o",pch=20,col="black")', file=RS)
                print('lines(position,X_count/total,type="o",pch=20,col="grey")', file=RS)
                print(
                    f"legend({read_len - 10},ym,"
                    'legend=c("A","T","G","C","N","X"),'
                    'col=c("dark green","red","blue","cyan","black","grey"),'
                    "lwd=2,pch=20,"
                    'text.col=c("dark green","red","blue","cyan","black","grey"))',
                    file=RS,
                )
                print("dev.off()", file=RS)
            else:
                print("total= A_count + C_count + G_count + T_count", file=RS)
                print("ym=max(A_count/total,C_count/total,G_count/total,T_count/total) + 0.05", file=RS)
                print("yn=min(A_count/total,C_count/total,G_count/total,T_count/total)", file=RS)

                print(f'pdf("{outfile}.NVC_plot.pdf")', file=RS)
                print(
                    'plot(position,A_count/total,type="o",pch=20,'
                    'ylim=c(yn,ym),col="dark green",'
                    'xlab="Position of Read",ylab="Nucleotide Frequency")',
                    file=RS,
                )
                print('lines(position,T_count/total,type="o",pch=20,col="red")', file=RS)
                print('lines(position,G_count/total,type="o",pch=20,col="blue")', file=RS)
                print('lines(position,C_count/total,type="o",pch=20,col="cyan")', file=RS)
                print(
                    f"legend({read_len - 10},ym,"
                    'legend=c("A","T","G","C"),'
                    'col=c("dark green","red","blue","cyan"),'
                    "lwd=2,pch=20,"
                    'text.col=c("dark green","red","blue","cyan"))',
                    file=RS,
                )
                print("dev.off()", file=RS)

    def readsQual_boxplot(self, outfile: str, shrink: int = 1000, q_cut: int = 30) -> None:
        """calculate phred quality score for each base in read (5->3)"""

        output = f"{outfile}.qual.r"
        with open(output, "w") as FO:
            if self.bam_format:
                print("Read BAM file ... ", end=" ", file=sys.stderr)
            else:
                print("Read SAM file ... ", end=" ", file=sys.stderr)

            quality: dict[int, dict[int, int]] = collections.defaultdict(dict)  # read_pos=>quality score=>count
            q_max = -1
            q_min = 10000
            q_list = []
            i_box = {}  # key is read postion,value is
            for aligned_read in _pysam_iter(self.samfile):
                if aligned_read.mapq < q_cut:
                    continue

                quals = aligned_read.query_qualities
                read_len = aligned_read.rlen
                if aligned_read.is_reverse:
                    quals = quals[::-1]

                for i, q in enumerate(quals):
                    if q > q_max:
                        q_max = q
                    if q < q_min:
                        q_min = q
                    try:
                        quality[i][q] += 1
                    except KeyError:
                        quality[i][q] = 1
            print("Done", file=sys.stderr)

            for p in range(0, read_len):
                val = []
                occurrence = []
                for q in range(q_min, q_max + 1):
                    if p in quality and q in quality[p]:
                        val.append(str(q))
                        occurrence.append(str(quality[p][q]))
                        q_list.append(str(quality[p][q]))
                    else:
                        q_list.append(str(0))
                i_box[p] = f"rep(c({','.join(val)}),times=c({','.join(occurrence)})/{shrink})"

            # generate R script for boxplot
            print(f"pdf('{outfile}.qual.boxplot.pdf')", file=FO)
            for i in sorted(i_box):
                print(f"p{i}<-{i_box[i]}", file=FO)
            boxplot_vars = ",".join([f"p{i}" for i in i_box])
            print(
                f'boxplot({boxplot_vars},xlab="Position of Read(5\'->3\')",ylab="Phred Quality Score",outline=F)',
                file=FO,
            )
            print("dev.off()", file=FO)

            # generate R script for heatmap
            print("\n", file=FO)
            print(f"pdf('{outfile}.qual.heatmap.pdf')", file=FO)
            print(f"qual=c({','.join(q_list)})", file=FO)
            print(f"mat=matrix(qual,ncol={read_len},byrow=F)", file=FO)
            print(
                "Lab.palette <- colorRampPalette("
                'c("blue", "orange", "red3","red2","red1","red"), '
                "space = \"rgb\",interpolate=c('spline'))",
                file=FO,
            )
            print(
                f"heatmap(mat,Rowv=NA,Colv=NA,"
                f'xlab="Position of Read",'
                f'ylab="Phred Quality Score",'
                f"labRow=seq(from={q_min},to={q_max}),"
                f'col = Lab.palette(256),scale="none" )',
                file=FO,
            )
            print("dev.off()", file=FO)

    def readGC(self, outfile: str, q_cut: int = 30) -> None:
        """GC content distribution of reads"""
        outfile1 = f"{outfile}.GC.xls"
        outfile2 = f"{outfile}.GC_plot.r"
        with open(outfile1, "w") as FO, open(outfile2, "w") as RS:
            gc_hist: dict[int, int] = collections.defaultdict(int)  # key is GC% * 100 (int), value is count

            if self.bam_format:
                print("Read BAM file ... ", end=" ", file=sys.stderr)
            else:
                print("Read SAM file ... ", end=" ", file=sys.stderr)

            for aligned_read in _pysam_iter(self.samfile):
                if aligned_read.is_unmapped:
                    continue
                if aligned_read.is_qcfail:
                    continue
                if aligned_read.mapq < q_cut:
                    continue
                RNA_read = aligned_read.seq.upper()
                gc_key = round((RNA_read.count("C") + RNA_read.count("G")) / len(RNA_read) * 10000)
                gc_hist[gc_key] += 1
            print("Done", file=sys.stderr)

            print("writing GC content ...", file=sys.stderr)
            print("GC%\tread_count", file=FO)
            for gc_key in gc_hist:
                gc_str = f"{gc_key / 100.0:4.2f}"
                print(f"{gc_str}\t{gc_hist[gc_key]}", file=FO)

            print("writing R script ...", file=sys.stderr)
            print(f'pdf("{outfile}.GC_plot.pdf")', file=RS)
            gc_strs = [f"{k / 100.0:4.2f}" for k in gc_hist]
            gc_vals = ",".join(str(v) for v in gc_hist.values())
            print(
                f"gc=rep(c({','.join(gc_strs)}),times=c({gc_vals}))",
                file=RS,
            )
            print(
                f"hist(gc,probability=T,breaks={100},"
                f'xlab="GC content (%)",ylab="Density of Reads",'
                f'border="blue",main="")',
                file=RS,
            )
            print("dev.off()", file=RS)

    def readDupRate(self, q_cut: int, outfile: str, up_bound: int = 500) -> None:
        """Calculate reads's duplicate rates"""
        outfile1 = f"{outfile}.seq.DupRate.xls"
        outfile2 = f"{outfile}.pos.DupRate.xls"
        outfile3 = f"{outfile}.DupRate_plot.r"
        with open(outfile1, "w") as SEQ, open(outfile2, "w") as POS, open(outfile3, "w") as RS:
            seqDup: dict[str, int] = collections.defaultdict(int)
            posDup: dict[str, int] = collections.defaultdict(int)

            seqDup_count: dict[int, int] = collections.defaultdict(int)
            posDup_count: dict[int, int] = collections.defaultdict(int)

            if self.bam_format:
                print("Load BAM file ... ", end=" ", file=sys.stderr)
            else:
                print("Load SAM file ... ", end=" ", file=sys.stderr)

            for aligned_read in _pysam_iter(self.samfile):
                exon_boundary = ""
                if aligned_read.is_unmapped:
                    continue  # skip unmapped read
                if aligned_read.is_qcfail:
                    continue  # skip low quality
                if aligned_read.mapq < q_cut:
                    continue
                RNA_read = aligned_read.seq.upper()
                seqDup[RNA_read] += 1  # key is read sequence

                chrom = self.samfile.getrname(aligned_read.tid)  # type: ignore[attr-defined]
                exon_blocks = aligned_read.get_blocks()
                exon_boundary = ":".join(f"{ex[0]}-{ex[1]}" for ex in exon_blocks)
                key = f"{chrom}:{aligned_read.pos}:{exon_boundary}"
                posDup[key] += 1
            print("Done", file=sys.stderr)

            print("report duplicte rate based on sequence ...", file=sys.stderr)
            print("Occurrence\tUniqReadNumber", file=SEQ)
            for i in seqDup.values():  # key is occurence, value is uniq reads number (based on seq)
                seqDup_count[i] += 1
            for k in sorted(seqDup_count.keys()):
                print(f"{k}\t{seqDup_count[k]}", file=SEQ)

            print("report duplicte rate based on mapping  ...", file=sys.stderr)
            print("Occurrence\tUniqReadNumber", file=POS)
            for i in posDup.values():  # key is occurence, value is uniq reads number (based on coord)
                posDup_count[i] += 1
            for k in sorted(posDup_count.keys()):
                print(f"{k}\t{posDup_count[k]}", file=POS)

            print("generate R script ...", file=sys.stderr)
            print(f"pdf('{outfile}.DupRate_plot.pdf')", file=RS)
            print("par(mar=c(5,4,4,5),las=0)", file=RS)
            seq_occ = ",".join([str(i) for i in sorted(seqDup_count.keys())])
            print(f"seq_occ=c({seq_occ})", file=RS)
            seq_uniq = ",".join([str(seqDup_count[i]) for i in sorted(seqDup_count.keys())])
            print(f"seq_uniqRead=c({seq_uniq})", file=RS)
            pos_occ = ",".join([str(i) for i in sorted(posDup_count.keys())])
            print(f"pos_occ=c({pos_occ})", file=RS)
            pos_uniq = ",".join([str(posDup_count[i]) for i in sorted(posDup_count.keys())])
            print(f"pos_uniqRead=c({pos_uniq})", file=RS)
            print(
                f"plot(pos_occ,log10(pos_uniqRead),"
                f"ylab='Number of Reads (log10)',"
                f"xlab='Occurrence of read',"
                f"pch=4,cex=0.8,col='blue',xlim=c(1,{up_bound}),yaxt='n')",
                file=RS,
            )
            print("points(seq_occ,log10(seq_uniqRead),pch=20,cex=0.8,col='red')", file=RS)
            print("ym=floor(max(log10(pos_uniqRead)))", file=RS)
            legend_x = max(up_bound - 200, 1)
            print(
                f"legend({legend_x},ym,legend=c('Sequence-based','Mapping-based'),col=c('blue','red'),pch=c(4,20))",
                file=RS,
            )
            print("axis(side=2,at=0:ym,labels=0:ym)", file=RS)
            print(
                "axis(side=4,"
                "at=c(log10(pos_uniqRead[1]),log10(pos_uniqRead[2]),"
                "log10(pos_uniqRead[3]),log10(pos_uniqRead[4])), "
                "labels=c("
                "round(pos_uniqRead[1]*100/sum(pos_uniqRead*pos_occ)),"
                "round(pos_uniqRead[2]*100/sum(pos_uniqRead*pos_occ)),"
                "round(pos_uniqRead[3]*100/sum(pos_uniqRead*pos_occ)),"
                "round(pos_uniqRead[4]*100/sum(pos_uniqRead*pos_occ))))",
                file=RS,
            )
            print('mtext(4, text = "Reads %", line = 2)', file=RS)
            print("dev.off()", file=RS)

    def clipping_profile(self, outfile: str, q_cut: int, PE: bool, type: str = "S") -> None:
        """calculate profile of soft clipping or insertion"""

        out_file1 = f"{outfile}.clipping_profile.xls"
        out_file2 = f"{outfile}.clipping_profile.r"
        with open(out_file1, "w") as OUT, open(out_file2, "w") as ROUT:
            print("Position\tClipped_nt\tNon_clipped_nt", file=OUT)

            if self.bam_format:
                print("Load BAM file ... ", end=" ", file=sys.stderr)
            else:
                print("Load SAM file ... ", end=" ", file=sys.stderr)

            # Map type character to CIGAR op code
            type_op = 4 if type == "S" else 1  # S=4, I=1
            last_read_len = 0

            # single end sequencing
            if PE is False:
                total_read = 0.0
                soft_clip_profile: dict[int, float] = collections.defaultdict(int)
                for aligned_read in _pysam_iter(self.samfile):
                    if not _passes_qc(aligned_read, q_cut):
                        continue

                    total_read += 1
                    cigar = aligned_read.cigar

                    # Single pass: compute read length and check for target op
                    read_len = 0
                    has_op = False
                    for c, s in cigar:
                        if c in (0, 1, 4, 7, 8):
                            read_len += s
                        if c == type_op:
                            has_op = True
                    last_read_len = read_len
                    if not has_op:
                        continue

                    is_reverse = aligned_read.is_reverse

                    pos = 0
                    for c, s in cigar:
                        if c in (0, 1, 4, 7, 8):  # ops that consume query
                            if c == type_op:
                                for p in range(pos, pos + s):
                                    indx = (read_len - 1 - p) if is_reverse else p
                                    soft_clip_profile[indx] += 1.0
                            pos += s
                print("Done", file=sys.stderr)

                print(f"Total reads used: {int(total_read)}", file=sys.stderr)
                read_pos = list(range(0, last_read_len))
                clip_count = []
                for i in read_pos:
                    print(
                        f"{i}\t{soft_clip_profile[i]}\t{total_read - soft_clip_profile[i]}",
                        file=OUT,
                    )
                    clip_count.append(soft_clip_profile[i])

                print(f'pdf("{outfile}.clipping_profile.pdf")', file=ROUT)
                rp = ",".join([str(i) for i in read_pos])
                print(f"read_pos=c({rp})", file=ROUT)
                cc = ",".join([str(i) for i in clip_count])
                print(f"clip_count=c({cc})", file=ROUT)
                print(f"nonclip_count= {int(total_read)} - clip_count", file=ROUT)
                print(
                    "plot(read_pos, "
                    "nonclip_count*100/(clip_count+nonclip_count),"
                    'col="blue",main="clipping profile",'
                    'xlab="Position of read",'
                    'ylab="Non-clipped %",type="b")',
                    file=ROUT,
                )
                print("dev.off()", file=ROUT)

            if PE is True:
                total_read1 = 0.0
                total_read2 = 0.0
                r1_soft_clip_profile: dict[int, float] = collections.defaultdict(int)
                r2_soft_clip_profile: dict[int, float] = collections.defaultdict(int)
                for aligned_read in _pysam_iter(self.samfile):
                    if not _passes_qc(aligned_read, q_cut):
                        continue
                    if not aligned_read.is_paired:
                        continue
                    if aligned_read.is_read1:
                        total_read1 += 1
                    if aligned_read.is_read2:
                        total_read2 += 1
                    cigar = aligned_read.cigar

                    # Single pass: compute read length and check for target op
                    read_len = 0
                    has_op = False
                    for c, s in cigar:
                        if c in (0, 1, 4, 7, 8):
                            read_len += s
                        if c == type_op:
                            has_op = True
                    last_read_len = read_len
                    if not has_op:
                        continue

                    is_reverse = aligned_read.is_reverse
                    target_profile = r1_soft_clip_profile if aligned_read.is_read1 else r2_soft_clip_profile

                    pos = 0
                    for c, s in cigar:
                        if c in (0, 1, 4, 7, 8):
                            if c == type_op:
                                for p in range(pos, pos + s):
                                    indx = (read_len - 1 - p) if is_reverse else p
                                    target_profile[indx] += 1.0
                            pos += s
                print("Done", file=sys.stderr)

                read_pos = list(range(0, last_read_len))
                r1_clip_count = []
                r2_clip_count = []

                print(f"Total read-1 used: {int(total_read1)}", file=sys.stderr)
                print(f"Total read-2 used: {int(total_read2)}", file=sys.stderr)
                print("Read-1:", file=OUT)
                for i in read_pos:
                    print(
                        f"{i}\t{r1_soft_clip_profile[i]}\t{total_read1 - r1_soft_clip_profile[i]}",
                        file=OUT,
                    )
                    r1_clip_count.append(r1_soft_clip_profile[i])

                print("Read-2:", file=OUT)
                for i in read_pos:
                    print(
                        f"{i}\t{r2_soft_clip_profile[i]}\t{total_read2 - r2_soft_clip_profile[i]}",
                        file=OUT,
                    )
                    r2_clip_count.append(r2_soft_clip_profile[i])

                rp = ",".join([str(i) for i in read_pos])
                print(f'pdf("{outfile}.clipping_profile.R1.pdf")', file=ROUT)
                print(f"read_pos=c({rp})", file=ROUT)
                r1cc = ",".join([str(i) for i in r1_clip_count])
                print(f"r1_clip_count=c({r1cc})", file=ROUT)
                print(f"r1_nonclip_count = {int(total_read1)} - r1_clip_count", file=ROUT)
                print(
                    "plot(read_pos, "
                    "r1_nonclip_count*100/(r1_clip_count + r1_nonclip_count),"
                    'col="blue",main="clipping profile",'
                    'xlab="Position of read (read-1)",'
                    'ylab="Non-clipped %",type="b")',
                    file=ROUT,
                )
                print("dev.off()", file=ROUT)

                print(f'pdf("{outfile}.clipping_profile.R2.pdf")', file=ROUT)
                print(f"read_pos=c({rp})", file=ROUT)
                r2cc = ",".join([str(i) for i in r2_clip_count])
                print(f"r2_clip_count=c({r2cc})", file=ROUT)
                print(f"r2_nonclip_count = {int(total_read2)} - r2_clip_count", file=ROUT)
                print(
                    "plot(read_pos, "
                    "r2_nonclip_count*100/(r2_clip_count + r2_nonclip_count),"
                    'col="blue",main="clipping profile",'
                    'xlab="Position of read (read-2)",'
                    'ylab="Non-clipped %",type="b")',
                    file=ROUT,
                )
                print("dev.off()", file=ROUT)

    def insertion_profile(self, outfile: str, q_cut: int, PE: bool, type: str = "I") -> None:
        """calculate profile of insertion"""

        out_file1 = f"{outfile}.insertion_profile.xls"
        out_file2 = f"{outfile}.insertion_profile.r"
        with open(out_file1, "w") as OUT, open(out_file2, "w") as ROUT:
            print("Position\tInsert_nt\tNon_insert_nt", file=OUT)

            if self.bam_format:
                print("Load BAM file ... ", end=" ", file=sys.stderr)
            else:
                print("Load SAM file ... ", end=" ", file=sys.stderr)

            # Map type character to CIGAR op code
            type_op = 4 if type == "S" else 1  # S=4, I=1
            last_read_len = 0

            # single end sequencing
            if PE is False:
                total_read = 0.0
                soft_clip_profile: dict[int, float] = collections.defaultdict(int)
                for aligned_read in _pysam_iter(self.samfile):
                    if not _passes_qc(aligned_read, q_cut):
                        continue

                    total_read += 1
                    cigar = aligned_read.cigar

                    # Single pass: compute read length and check for target op
                    read_len = 0
                    has_op = False
                    for c, s in cigar:
                        if c in (0, 1, 4, 7, 8):
                            read_len += s
                        if c == type_op:
                            has_op = True
                    last_read_len = read_len
                    if not has_op:
                        continue

                    is_reverse = aligned_read.is_reverse

                    pos = 0
                    for c, s in cigar:
                        if c in (0, 1, 4, 7, 8):
                            if c == type_op:
                                for p in range(pos, pos + s):
                                    indx = (read_len - 1 - p) if is_reverse else p
                                    soft_clip_profile[indx] += 1.0
                            pos += s
                print("Done", file=sys.stderr)

                print(f"Total reads used: {int(total_read)}", file=sys.stderr)
                read_pos = list(range(0, last_read_len))
                clip_count = []
                for i in read_pos:
                    print(
                        f"{i}\t{soft_clip_profile[i]}\t{total_read - soft_clip_profile[i]}",
                        file=OUT,
                    )
                    clip_count.append(soft_clip_profile[i])

                print(f'pdf("{outfile}.insertion_profile.pdf")', file=ROUT)
                rp = ",".join([str(i) for i in read_pos])
                print(f"read_pos=c({rp})", file=ROUT)
                ic = ",".join([str(i) for i in clip_count])
                print(f"insert_count=c({ic})", file=ROUT)
                print(f"noninsert_count= {int(total_read)} - insert_count", file=ROUT)
                print(
                    "plot(read_pos, "
                    "insert_count*100/(insert_count+noninsert_count),"
                    'col="blue",main="Insertion profile",'
                    'xlab="Position of read",'
                    'ylab="Insertion %",type="b")',
                    file=ROUT,
                )
                print("dev.off()", file=ROUT)

            if PE is True:
                total_read1 = 0.0
                total_read2 = 0.0
                r1_soft_clip_profile: dict[int, float] = collections.defaultdict(int)
                r2_soft_clip_profile: dict[int, float] = collections.defaultdict(int)
                for aligned_read in _pysam_iter(self.samfile):
                    if not _passes_qc(aligned_read, q_cut):
                        continue
                    if not aligned_read.is_paired:
                        continue
                    if aligned_read.is_read1:
                        total_read1 += 1
                    if aligned_read.is_read2:
                        total_read2 += 1
                    cigar = aligned_read.cigar

                    # Single pass: compute read length and check for target op
                    read_len = 0
                    has_op = False
                    for c, s in cigar:
                        if c in (0, 1, 4, 7, 8):
                            read_len += s
                        if c == type_op:
                            has_op = True
                    last_read_len = read_len
                    if not has_op:
                        continue

                    is_reverse = aligned_read.is_reverse
                    target_profile = r1_soft_clip_profile if aligned_read.is_read1 else r2_soft_clip_profile

                    pos = 0
                    for c, s in cigar:
                        if c in (0, 1, 4, 7, 8):
                            if c == type_op:
                                for p in range(pos, pos + s):
                                    indx = (read_len - 1 - p) if is_reverse else p
                                    target_profile[indx] += 1.0
                            pos += s
                print("Done", file=sys.stderr)

                read_pos = list(range(0, last_read_len))
                r1_clip_count = []
                r2_clip_count = []

                print(f"Total read-1 used: {int(total_read1)}", file=sys.stderr)
                print(f"Total read-2 used: {int(total_read2)}", file=sys.stderr)
                print("Read-1:", file=OUT)
                for i in read_pos:
                    print(
                        f"{i}\t{r1_soft_clip_profile[i]}\t{total_read1 - r1_soft_clip_profile[i]}",
                        file=OUT,
                    )
                    r1_clip_count.append(r1_soft_clip_profile[i])

                print("Read-2:", file=OUT)
                for i in read_pos:
                    print(
                        f"{i}\t{r2_soft_clip_profile[i]}\t{total_read2 - r2_soft_clip_profile[i]}",
                        file=OUT,
                    )
                    r2_clip_count.append(r2_soft_clip_profile[i])

                rp = ",".join([str(i) for i in read_pos])
                print(f'pdf("{outfile}.insertion_profile.R1.pdf")', file=ROUT)
                print(f"read_pos=c({rp})", file=ROUT)
                r1ic = ",".join([str(i) for i in r1_clip_count])
                print(f"r1_insert_count=c({r1ic})", file=ROUT)
                print(f"r1_noninsert_count = {int(total_read1)} - r1_insert_count", file=ROUT)
                print(
                    "plot(read_pos, "
                    "r1_insert_count*100/"
                    "(r1_insert_count + r1_noninsert_count),"
                    'col="blue",main="Insertion profile",'
                    'xlab="Position of read (read-1)",'
                    'ylab="Insertion %",type="b")',
                    file=ROUT,
                )
                print("dev.off()", file=ROUT)

                print(f'pdf("{outfile}.insertion_profile.R2.pdf")', file=ROUT)
                print(f"read_pos=c({rp})", file=ROUT)
                r2ic = ",".join([str(i) for i in r2_clip_count])
                print(f"r2_insert_count=c({r2ic})", file=ROUT)
                print(f"r2_noninsert_count = {int(total_read2)} - r2_insert_count", file=ROUT)
                print(
                    "plot(read_pos, "
                    "r2_insert_count*100/"
                    "(r2_insert_count + r2_noninsert_count),"
                    'col="blue",main="Insertion profile",'
                    'xlab="Position of read (read-2)",'
                    'ylab="Insertion %",type="b")',
                    file=ROUT,
                )
                print("dev.off()", file=ROUT)

    def mRNA_inner_distance(
        self,
        outfile: str,
        refbed: str,
        low_bound: int = 0,
        up_bound: int = 1000,
        step: int = 10,
        sample_size: int = 1000000,
        q_cut: int = 30,
    ) -> None:
        """estimate the inner distance of mRNA pair end fragment. fragment size = insert_size + 2 x read_length"""

        out_file1 = f"{outfile}.inner_distance.txt"
        out_file2 = f"{outfile}.inner_distance_freq.txt"
        out_file3 = f"{outfile}.inner_distance_plot.r"

        with open(out_file1, "w") as FO, open(out_file2, "w") as FQ, open(out_file3, "w") as RS:
            fchrom = "chr100"  # this is the fake chromosome
            ranges = {}
            ranges[fchrom] = Intersecter()

            window_left_bound = list(range(low_bound, up_bound, step))

            inner_distance_bitsets = BinnedBitSet()
            tmp = BinnedBitSet()
            tmp.set_range(0, 0)
            pair_num = 0
            sizes = []
            counts = []
            count = 0

            print(f"Get exon regions from {refbed} ...", file=sys.stderr)
            bed_obj = BED.ParseBED(refbed)
            ref_exons = []

            for exn in bed_obj.getExon():
                ref_exons.append([exn[0].upper(), exn[1], exn[2]])
            exon_bitsets = binned_bitsets_from_list(ref_exons)

            transcript_ranges = {}
            for i_chr, i_st, i_end, i_strand, i_name in bed_obj.getTranscriptRanges():
                i_chr = i_chr.upper()  # type: ignore[union-attr]
                if i_chr not in transcript_ranges:
                    transcript_ranges[i_chr] = Intersecter()
                else:
                    transcript_ranges[i_chr].add_interval(Interval(i_st, i_end, value=i_name))

            if self.bam_format:
                print("Load BAM file ... ", end=" ", file=sys.stderr)
            else:
                print("Load SAM file ... ", end=" ", file=sys.stderr)

            for aligned_read in _pysam_iter(self.samfile):
                if pair_num >= sample_size:
                    break
                splice_intron_size = 0
                if not _passes_qc(aligned_read, q_cut):
                    continue
                if not aligned_read.is_paired:
                    continue
                if aligned_read.mate_is_unmapped:
                    continue

                read1_len = aligned_read.qlen
                read1_start = aligned_read.pos
                read2_start = aligned_read.mpos  # 0-based, not included
                if read2_start < read1_start:
                    continue  # because BAM file is sorted, mate_read is already processed if its coordinate is smaller
                if read2_start == read1_start and aligned_read.is_read1:
                    continue

                pair_num += 1

                # check if reads were mapped to diff chromosomes
                R_read1_ref = self.samfile.get_reference_name(aligned_read.tid)  # type: ignore[attr-defined]
                R_read2_ref = self.samfile.get_reference_name(aligned_read.rnext)  # type: ignore[attr-defined]
                if R_read1_ref != R_read2_ref:
                    FO.write(f"{aligned_read.qname}\tNA\tsameChrom=No\n")  # reads mapped to different chromosomes
                    continue

                chrom = R_read1_ref.upper()
                intron_blocks = bam_cigar.fetch_intron(read1_start, aligned_read.cigar)
                for intron in intron_blocks:
                    splice_intron_size += intron[1] - intron[0]
                read1_end = read1_start + read1_len + splice_intron_size

                if read2_start >= read1_end:
                    inner_distance = read2_start - read1_end
                else:
                    overlap = 0
                    exon_blocks = aligned_read.get_blocks()
                    for ex in exon_blocks:
                        # Count positions in (read2_start, read1_end] that fall within exon (ex[0]+1, ex[1]]
                        lo = max(ex[0] + 1, read2_start + 1)
                        hi = min(ex[1], read1_end)
                        if hi >= lo:
                            overlap += hi - lo + 1
                    inner_distance = -overlap

                read1_gene_names = set()  # read1_end
                try:
                    for gene in transcript_ranges[chrom].find(
                        read1_end - 1, read1_end
                    ):  # gene: Interval(0, 10, value=a)
                        read1_gene_names.add(gene.value)
                except KeyError:
                    pass

                read2_gene_names = set()  # read2_start
                try:
                    for gene in transcript_ranges[chrom].find(
                        read2_start, read2_start + 1
                    ):  # gene: Interval(0, 10, value=a)
                        read2_gene_names.add(gene.value)
                except KeyError:
                    pass

                if len(read1_gene_names.intersection(read2_gene_names)) == 0:  # no common gene
                    FO.write(
                        f"{aligned_read.qname}\t{inner_distance}\tsameTranscript=No,dist=genomic\n"
                    )  # reads mapped to different gene
                    ranges[fchrom].add_interval(Interval(inner_distance - 1, inner_distance))
                    continue

                if inner_distance > 0:
                    if chrom in exon_bitsets:
                        size = 0
                        inner_distance_bitsets.set_range(read1_end, read2_start - read1_end)
                        inner_distance_bitsets.iand(exon_bitsets[chrom])
                        end = 0
                        while True:
                            start = inner_distance_bitsets.next_set(end)
                            if start == inner_distance_bitsets.size:
                                break
                            end = inner_distance_bitsets.next_clear(start)
                            size += end - start
                        inner_distance_bitsets.iand(tmp)  # clear BinnedBitSet

                        if size == inner_distance:
                            FO.write(f"{aligned_read.qname}\t{size}\tsameTranscript=Yes,sameExon=Yes,dist=mRNA\n")
                            ranges[fchrom].add_interval(Interval(size - 1, size))
                        elif size > 0 and size < inner_distance:
                            FO.write(f"{aligned_read.qname}\t{size}\tsameTranscript=Yes,sameExon=No,dist=mRNA\n")
                            ranges[fchrom].add_interval(Interval(size - 1, size))
                        elif size <= 0:
                            FO.write(
                                f"{aligned_read.qname}\t{inner_distance}"
                                f"\tsameTranscript=Yes,nonExonic=Yes,dist=genomic\n"
                            )
                            ranges[fchrom].add_interval(Interval(inner_distance - 1, inner_distance))
                    else:
                        FO.write(f"{aligned_read.qname}\t{inner_distance}\tunknownChromosome,dist=genomic\n")
                        ranges[fchrom].add_interval(Interval(inner_distance - 1, inner_distance))
                else:
                    FO.write(f"{aligned_read.qname}\t{inner_distance}\treadPairOverlap\n")
                    ranges[fchrom].add_interval(Interval(inner_distance - 1, inner_distance))
            print("Done", file=sys.stderr)

            print(f"Total read pairs  used {pair_num}", file=sys.stderr)
            if pair_num == 0:
                print("Cannot find paired reads", file=sys.stderr)
                sys.exit(1)

            for st in window_left_bound:
                sizes.append(str(st + step / 2))
                count = str(len(ranges[fchrom].find(st, st + step)))  # type: ignore[assignment]
                counts.append(count)
                print(f"{st}\t{st + step}\t{count}", file=FQ)  # type: ignore[operator]
            print(f"out_file = '{outfile}'", file=RS)
            print(f"pdf('{outfile}.inner_distance_plot.pdf')", file=RS)
            sz = ",".join(sizes)  # type: ignore[arg-type]
            ct = ",".join(counts)  # type: ignore[arg-type]
            print(f"fragsize=rep(c({sz}),times=c({ct}))", file=RS)
            print("frag_sd = sd(fragsize)", file=RS)
            print("frag_mean = mean(fragsize)", file=RS)
            print("frag_median = median(fragsize)", file=RS)
            print('write(x=c("Name","Mean","Median","sd"), sep="\t", file=stdout(),ncolumns=4)', file=RS)
            print('write(c(out_file,frag_mean,frag_median,frag_sd),sep="\t", file=stdout(),ncolumns=4)', file=RS)
            n_breaks = len(window_left_bound)
            print(
                f"hist(fragsize,probability=T,breaks={n_breaks},"
                f'xlab="mRNA insert size (bp)",'
                f'main=paste(c("Mean=",frag_mean,";",'
                f'"SD=",frag_sd),collapse=""),'
                f'border="blue")',
                file=RS,
            )
            print(f"lines(density(fragsize,bw={2 * step}),col='red')", file=RS)
            print("dev.off()", file=RS)

    def annotate_junction(self, refgene: str | None, outfile: str, min_intron: int = 50, q_cut: int = 30) -> None:
        """Annotate splicing junctions in BAM or SAM file. Note that a (long) read might have multiple splicing
        events  (splice multiple times), and the same splicing events can be consolidated into a single
        junction"""

        out_file = f"{outfile}.junction.xls"
        out_file2 = f"{outfile}.junction_plot.r"
        if refgene is None:
            print("You must provide reference gene model in bed format.", file=sys.stderr)
            sys.exit(1)
        with open(out_file, "w") as OUT, open(out_file2, "w") as ROUT:
            # reading reference gene model
            refIntronStarts: dict[str, dict[int, int]] = collections.defaultdict(dict)
            refIntronEnds: dict[str, dict[int, int]] = collections.defaultdict(dict)
            total_junc = 0
            novel35_junc = 0
            novel3or5_junc = 0
            known_junc = 0
            filtered_junc = 0
            splicing_events: dict[str, int] = collections.defaultdict(int)

            print("Reading reference bed file: ", refgene, " ... ", end=" ", file=sys.stderr)
            with open(refgene, "r") as _fh:
                for line in _fh:
                    if line.startswith(("#", "track", "browser")):
                        continue
                    # Parse fields from gene tabls
                    fields = line.split()
                    if len(fields) < 12:
                        print("Invalid bed line (skipped):", line, end=" ", file=sys.stderr)
                        continue
                    chrom = fields[0].upper()
                    tx_start = int(fields[1])
                    if int(fields[9]) == 1:
                        continue

                    exon_starts = list(map(int, fields[11].rstrip(",\n").split(",")))
                    exon_starts = [x + tx_start for x in exon_starts]
                    exon_ends = list(map(int, fields[10].rstrip(",\n").split(",")))
                    exon_ends = [x + y for x, y in zip(exon_starts, exon_ends)]
                    intron_start = exon_ends[:-1]
                    intron_end = exon_starts[1:]
                    for i_st, i_end in zip(intron_start, intron_end):
                        refIntronStarts[chrom][i_st] = i_st
                        refIntronEnds[chrom][i_end] = i_end
            print("Done", file=sys.stderr)

            # reading input SAM file
            if self.bam_format:
                print("Load BAM file ... ", end=" ", file=sys.stderr)
            else:
                print("Load SAM file ... ", end=" ", file=sys.stderr)

            for aligned_read in _pysam_iter(self.samfile):
                if not _passes_qc(aligned_read, q_cut):
                    continue

                chrom = self.samfile.getrname(aligned_read.tid).upper()  # type: ignore[attr-defined]
                hit_st = aligned_read.pos
                intron_blocks = bam_cigar.fetch_intron(hit_st, aligned_read.cigar)
                if len(intron_blocks) == 0:
                    continue
                for intrn in intron_blocks:
                    total_junc += 1
                    if intrn[1] - intrn[0] < min_intron:
                        filtered_junc += 1
                        continue
                    splicing_events[f"{chrom}:{intrn[0]}:{intrn[1]}"] += 1
                    if intrn[0] in refIntronStarts[chrom] and intrn[1] in refIntronEnds[chrom]:
                        known_junc += 1  # known both
                    elif intrn[0] not in refIntronStarts[chrom] and intrn[1] not in refIntronEnds[chrom]:
                        novel35_junc += 1
                    else:
                        novel3or5_junc += 1
            print("Done", file=sys.stderr)

            print(f"total = {total_junc}")
            if total_junc == 0:
                print("No splice junction found.", file=sys.stderr)
                sys.exit(1)

            print(f'pdf("{outfile}.splice_events.pdf")', file=ROUT)
            evt_pcts = ",".join([str(i * 100.0 / total_junc) for i in (novel3or5_junc, novel35_junc, known_junc)])
            print(f"events=c({evt_pcts})", file=ROUT)
            pn_pct = round(novel3or5_junc * 100.0 / total_junc)
            cn_pct = round(novel35_junc * 100.0 / total_junc)
            k_pct = round(known_junc * 100.0 / total_junc)
            print(
                f"pie(events,col=c(2,3,4),init.angle=30,"
                f"angle=c(60,120,150),density=c(70,70,70),"
                f'main="splicing events",'
                f'labels=c("partial_novel {pn_pct}%",'
                f'"complete_novel {cn_pct}%","known {k_pct}%"))',
                file=ROUT,
            )
            print("dev.off()", file=ROUT)

            print("\n===================================================================", file=sys.stderr)
            print(f"Total splicing  Events:\t{total_junc}", file=sys.stderr)
            print(f"Known Splicing Events:\t{known_junc}", file=sys.stderr)
            print(f"Partial Novel Splicing Events:\t{novel3or5_junc}", file=sys.stderr)
            print(f"Novel Splicing Events:\t{novel35_junc}", file=sys.stderr)
            print(f"Filtered Splicing Events:\t{filtered_junc}", file=sys.stderr)

            # reset variables
            total_junc = 0
            novel35_junc = 0
            novel3or5_junc = 0
            known_junc = 0

            print("chrom\tintron_st(0-based)\tintron_end(1-based)\tread_count\tannotation", file=OUT)
            for i in splicing_events:
                total_junc += 1
                (chrom, i_st, i_end) = i.split(":")  # type: ignore[assignment]
                print(
                    f"{chrom.replace('CHR', 'chr')}\t{i_st}\t{i_end}\t{splicing_events[i]}\t",
                    end=" ",
                    file=OUT,
                )
                i_st = int(i_st)
                i_end = int(i_end)
                if i_st in refIntronStarts[chrom] and i_end in refIntronEnds[chrom]:
                    print("annotated", file=OUT)
                    known_junc += 1
                elif i_st not in refIntronStarts[chrom] and i_end not in refIntronEnds[chrom]:
                    print("complete_novel", file=OUT)
                    novel35_junc += 1
                else:
                    print("partial_novel", file=OUT)
                    novel3or5_junc += 1

            if total_junc == 0:
                print("No splice read found", file=sys.stderr)
                sys.exit(1)
            print(f"\nTotal splicing  Junctions:\t{total_junc}", file=sys.stderr)
            print(f"Known Splicing Junctions:\t{known_junc}", file=sys.stderr)
            print(f"Partial Novel Splicing Junctions:\t{novel3or5_junc}", file=sys.stderr)
            print(f"Novel Splicing Junctions:\t{novel35_junc}", file=sys.stderr)
            print("\n===================================================================", file=sys.stderr)

            print(f'pdf("{outfile}.splice_junction.pdf")', file=ROUT)
            junc_pcts = ",".join([str(i * 100.0 / total_junc) for i in (novel3or5_junc, novel35_junc, known_junc)])
            print(f"junction=c({junc_pcts})", file=ROUT)
            pn_pct = round(novel3or5_junc * 100.0 / total_junc)
            cn_pct = round(novel35_junc * 100.0 / total_junc)
            k_pct = round(known_junc * 100.0 / total_junc)
            print(
                f"pie(junction,col=c(2,3,4),init.angle=30,"
                f"angle=c(60,120,150),density=c(70,70,70),"
                f'main="splicing junctions",'
                f'labels=c("partial_novel {pn_pct}%",'
                f'"complete_novel {cn_pct}%","known {k_pct}%"))',
                file=ROUT,
            )
            print("dev.off()", file=ROUT)

    def saturation_junction(
        self,
        refgene: str | None,
        outfile: str | None = None,
        sample_start: int = 5,
        sample_step: int = 5,
        sample_end: int = 100,
        min_intron: int = 50,
        recur: int = 1,
        q_cut: int = 30,
    ) -> None:
        """check if an RNA-seq experiment is saturated in terms of detecting known splicing junction"""

        out_file = f"{outfile}.junctionSaturation_plot.r"
        if refgene is None:
            print("You must provide reference gene model in bed format.", file=sys.stderr)
            sys.exit(1)

        with open(out_file, "w") as OUT:
            # reading reference gene
            knownSpliceSites = set()
            chrom_list = set()
            print("reading reference bed file: ", refgene, " ... ", end=" ", file=sys.stderr)
            with open(refgene, "r") as _fh:
                for line in _fh:
                    if line.startswith(("#", "track", "browser")):
                        continue
                    fields = line.split()
                    if len(fields) < 12:
                        print("Invalid bed line (skipped):", line, end=" ", file=sys.stderr)
                        continue
                    chrom = fields[0].upper()
                    chrom_list.add(chrom)
                    tx_start = int(fields[1])
                    if int(fields[9]) == 1:
                        continue

                    exon_starts = list(map(int, fields[11].rstrip(",\n").split(",")))
                    exon_starts = [x + tx_start for x in exon_starts]
                    exon_ends = list(map(int, fields[10].rstrip(",\n").split(",")))
                    exon_ends = [x + y for x, y in zip(exon_starts, exon_ends)]
                    intron_start = exon_ends[:-1]
                    intron_end = exon_starts[1:]
                    for st, end in zip(intron_start, intron_end):
                        knownSpliceSites.add(f"{chrom}:{st}-{end}")
            print(f"Done! Total {len(knownSpliceSites)} known splicing junctions.", file=sys.stderr)

            # read SAM file
            samSpliceSites = []
            intron_start = []
            intron_end = []
            uniqSpliceSites: dict[str, int] = collections.defaultdict(int)

            if self.bam_format:
                print("Load BAM file ... ", end=" ", file=sys.stderr)
            else:
                print("Load SAM file ... ", end=" ", file=sys.stderr)
            for aligned_read in _pysam_iter(self.samfile):
                if not _passes_qc(aligned_read, q_cut):
                    continue
                try:
                    chrom = self.samfile.getrname(aligned_read.tid).upper()  # type: ignore[attr-defined]
                except (ValueError, KeyError, AttributeError):
                    continue
                if chrom not in chrom_list:
                    continue

                hit_st = aligned_read.pos
                intron_blocks = bam_cigar.fetch_intron(hit_st, aligned_read.cigar)
                if len(intron_blocks) == 0:
                    continue
                for intrn in intron_blocks:
                    if intrn[1] - intrn[0] < min_intron:
                        continue
                    samSpliceSites.append(f"{chrom}:{intrn[0]}-{intrn[1]}")
            print("Done", file=sys.stderr)

            print("shuffling alignments ...", end=" ", file=sys.stderr)
            random.shuffle(samSpliceSites)
            print("Done", file=sys.stderr)

            # resampling
            SR_num = len(samSpliceSites)
            sample_size = 0
            known_junctionNum = 0
            unknown_junctionNum = 0
            known_junc = []
            all_junc = []
            unknown_junc = []
            # =========================sampling uniquely mapped reads from population
            tmp = list(range(sample_start, sample_end, sample_step))
            tmp.append(100)
            for pertl in tmp:  # [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95,100]
                index_st = int(SR_num * ((pertl - sample_step) / 100.0))
                index_end = int(SR_num * (pertl / 100.0))
                if index_st < 0:
                    index_st = 0
                sample_size += index_end - index_st

                print(f"sampling {pertl}% ({sample_size}) splicing reads.", end=" ", file=sys.stderr)

                # Incrementally update counts as new splice sites are added
                for i in range(index_st, index_end):
                    sj = samSpliceSites[i]
                    old_count = uniqSpliceSites[sj]
                    uniqSpliceSites[sj] = old_count + 1
                    if old_count == 0:
                        # Brand new junction
                        if sj not in knownSpliceSites:
                            unknown_junctionNum += 1
                    # Junction just crossed the recur threshold → count as known
                    if sj in knownSpliceSites and old_count < recur <= old_count + 1:
                        known_junctionNum += 1

                all_junctionNum = len(uniqSpliceSites)
                all_junc.append(str(all_junctionNum))
                print(f"{all_junctionNum} splicing junctions.", end=" ", file=sys.stderr)

                print(f"{known_junctionNum} known splicing junctions.", end=" ", file=sys.stderr)
                known_junc.append(str(known_junctionNum))

                unknown_junc.append(str(unknown_junctionNum))
                print(f"{unknown_junctionNum} novel splicing junctions.", file=sys.stderr)

            print(f"pdf('{outfile}.junctionSaturation_plot.pdf')", file=OUT)
            print(f"x=c({','.join([str(i) for i in tmp])})", file=OUT)
            print(f"y=c({','.join(known_junc)})", file=OUT)
            print(f"z=c({','.join(all_junc)})", file=OUT)
            print(f"w=c({','.join(unknown_junc)})", file=OUT)
            print(
                f"m=max({int(known_junc[-1]) // 1000},{int(all_junc[-1]) // 1000},{int(unknown_junc[-1]) // 1000})",
                file=OUT,
            )
            print(
                f"n=min({int(known_junc[0]) // 1000},{int(all_junc[0]) // 1000},{int(unknown_junc[0]) // 1000})",
                file=OUT,
            )
            print(
                "plot(x,z/1000,"
                "xlab='percent of total reads',"
                "ylab='Number of splicing junctions (x1000)',"
                "type='o',col='blue',ylim=c(n,m))",
                file=OUT,
            )
            print("points(x,y/1000,type='o',col='red')", file=OUT)
            print("points(x,w/1000,type='o',col='green')", file=OUT)
            print(
                f"legend(5,{int(all_junc[-1]) // 1000}, "
                f'legend=c("All junctions","known junctions",'
                f' "novel junctions"),'
                f'col=c("blue","red","green"),lwd=1,pch=1)',
                file=OUT,
            )
            print("dev.off()", file=OUT)

    def saturation_RPKM(
        self,
        refbed: str | None,
        outfile: str,
        sample_start: int = 5,
        sample_step: int = 5,
        sample_end: int = 100,
        skip_multi: bool = True,
        strand_rule: str | None = None,
        q_cut: int = 30,
    ) -> None:
        """for each gene, check if its RPKM (epxresion level) has already been saturated or not"""

        if refbed is None:
            print("You must specify a bed file representing gene model\n", file=sys.stderr)
            sys.exit(1)
        rpkm_file = f"{outfile}.eRPKM.xls"
        raw_file = f"{outfile}.rawCount.xls"

        with open(rpkm_file, "w") as RPKM_OUT, open(raw_file, "w") as RAW_OUT:
            cUR_num = 0  # number of fragments
            cUR_plus = 0
            cUR_minus = 0
            midpoints_plus: list[tuple[str, int]] = []
            midpoints_minus: list[tuple[str, int]] = []
            midpoints: list[tuple[str, int]] = []
            strandRule = _parse_strand_rule(strand_rule)

            # read SAM or BAM
            if self.bam_format:
                print("Load BAM file ... ", end=" ", file=sys.stderr)
            else:
                print("Load SAM file ... ", end=" ", file=sys.stderr)
            for aligned_read in _pysam_iter(self.samfile):
                if aligned_read.is_qcfail:
                    continue
                if aligned_read.is_duplicate:
                    continue
                if aligned_read.is_secondary:
                    continue
                if aligned_read.is_unmapped:
                    continue
                if skip_multi and aligned_read.mapq < q_cut:
                    continue
                chrom = self.samfile.getrname(aligned_read.tid).upper()  # type: ignore[attr-defined]

                # determine read_id and read_strand
                if aligned_read.is_paired:  # pair end
                    if aligned_read.is_read1:
                        read_id = "1"
                    if aligned_read.is_read2:
                        read_id = "2"
                else:
                    read_id = ""  # single end

                if aligned_read.is_reverse:
                    map_strand = "-"
                else:
                    map_strand = "+"
                strand_key = read_id + map_strand  # used to determine if a read should assign to gene(+) or gene(-)

                exon_blocks = aligned_read.get_blocks()
                cUR_num += len(exon_blocks)

                # strand specific
                if strand_rule is not None:
                    if strandRule[strand_key] == "+":
                        cUR_plus += len(exon_blocks)
                    if strandRule[strand_key] == "-":
                        cUR_minus += len(exon_blocks)
                    for exn in exon_blocks:
                        mid = exn[0] + (exn[1] - exn[0]) // 2
                        if strandRule[strand_key] == "+":
                            midpoints_plus.append((chrom, mid))
                        if strandRule[strand_key] == "-":
                            midpoints_minus.append((chrom, mid))
                # Not strand specific
                else:
                    for exn in exon_blocks:
                        midpoints.append((chrom, exn[0] + (exn[1] - exn[0]) // 2))
            print("Done", file=sys.stderr)

            print("shuffling alignments ...", end=" ", file=sys.stderr)
            random.shuffle(midpoints_plus)
            random.shuffle(midpoints_minus)
            random.shuffle(midpoints)
            print("Done", file=sys.stderr)

            accum: dict[str, list[int]] = collections.defaultdict(list)
            accum_plus: dict[str, list[int]] = collections.defaultdict(list)
            accum_minus: dict[str, list[int]] = collections.defaultdict(list)
            sample_size: float = 0
            RPKM_table = collections.defaultdict(list)
            rawCount_table = collections.defaultdict(list)
            RPKM_head = ["#chr", "start", "end", "name", "score", "strand"]

            # Parse gene model once before sampling loop (#10 perf fix)
            gene_models: list[tuple[str, str, list[tuple[int, int]], int]] = []
            with open(refbed, "r") as _fh:
                for line in _fh:
                    try:
                        if line.startswith(("#", "track", "browser")):
                            continue
                        fields = line.split()
                        chrom = fields[0].upper()
                        tx_start = int(fields[1])
                        tx_end = int(fields[2])
                        geneName = fields[3]
                        strand = fields[5]
                        exon_starts = list(map(int, fields[11].rstrip(",\n").split(",")))
                        exon_starts = [x + tx_start for x in exon_starts]
                        exon_sizes = list(map(int, fields[10].rstrip(",\n").split(",")))
                        exon_ends = [x + y for x, y in zip(exon_starts, exon_sizes)]
                        key = "\t".join((chrom.lower(), str(tx_start), str(tx_end), geneName, "0", strand))
                        mRNA_len = sum(exon_sizes)
                    except (IndexError, ValueError):
                        print(f"[NOTE:input bed must be 12-column] skipped this line: {line}", file=sys.stderr)
                        continue
                    if mRNA_len == 0:
                        print(f"{geneName} has 0 nucleotides. Exit!", file=sys.stderr)
                        sys.exit(1)
                    gene_models.append((key, strand, list(zip(exon_starts, exon_ends)), mRNA_len))

            tmp = list(range(sample_start, sample_end, sample_step))
            tmp.append(100)
            # =========================sampling uniquely mapped reads from population
            for pertl in tmp:  # [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95,100]
                percent_st = (pertl - sample_step) / 100.0
                percent_end = pertl / 100.0
                if percent_st < 0:
                    percent_st = 0
                sample_size = cUR_num * percent_end
                RPKM_head.append(f"{pertl}%")

                if strand_rule is not None:
                    print(
                        f"sampling {pertl}% ({int(cUR_plus * percent_end)}) forward strand fragments ...",
                        file=sys.stderr,
                    )
                    for ch, coord in midpoints_plus[int(cUR_plus * percent_st) : int(cUR_plus * percent_end)]:
                        accum_plus[ch].append(coord)

                    print(
                        f"sampling {pertl}% ({int(cUR_minus * percent_end)}) reverse strand fragments ...",
                        file=sys.stderr,
                    )
                    for ch, coord in midpoints_minus[int(cUR_minus * percent_st) : int(cUR_minus * percent_end)]:
                        accum_minus[ch].append(coord)

                else:
                    print(f"sampling {pertl}% ({int(sample_size)}) fragments ...", file=sys.stderr)
                    for ch, coord in midpoints[int(cUR_num * percent_st) : int(cUR_num * percent_end)]:
                        accum[ch].append(coord)

                # Build sorted numpy arrays for interval counting
                sorted_accum = {ch: np.sort(np.array(coords)) for ch, coords in accum.items()}
                sorted_plus = {ch: np.sort(np.array(coords)) for ch, coords in accum_plus.items()}
                sorted_minus = {ch: np.sort(np.array(coords)) for ch, coords in accum_minus.items()}

                # ========================= calculating RPKM based on sub-population
                print(f"assign reads to transcripts in {refbed} ...", file=sys.stderr)
                for key, strand, exon_intervals, mRNA_len in gene_models:
                    chrom = key.split("\t")[0].upper()
                    mRNA_count = 0
                    for st, end in exon_intervals:
                        if strand_rule is not None:
                            if (strand == "+") and (chrom in sorted_plus):
                                arr = sorted_plus[chrom]
                                mRNA_count += int(np.searchsorted(arr, end, "left") - np.searchsorted(arr, st, "left"))
                            if (strand == "-") and (chrom in sorted_minus):
                                arr = sorted_minus[chrom]
                                mRNA_count += int(np.searchsorted(arr, end, "left") - np.searchsorted(arr, st, "left"))
                        else:
                            if chrom in sorted_accum:
                                arr = sorted_accum[chrom]
                                mRNA_count += int(np.searchsorted(arr, end, "left") - np.searchsorted(arr, st, "left"))
                    if sample_size == 0:
                        print("Too few reads to sample. Exit!", file=sys.stderr)
                        sys.exit(1)
                    mRNA_RPKM = (mRNA_count * 1000000000.0) / (mRNA_len * sample_size)
                    RPKM_table[key].append(str(mRNA_RPKM))
                    rawCount_table[key].append(str(mRNA_count))
                print("", file=sys.stderr)

            print("\t".join(RPKM_head), file=RPKM_OUT)
            print("\t".join(RPKM_head), file=RAW_OUT)
            for key in RPKM_table:
                print(f"{key}\t", end=" ", file=RPKM_OUT)
                print("\t".join(RPKM_table[key]), file=RPKM_OUT)
                print(f"{key}\t", end=" ", file=RAW_OUT)
                print("\t".join(rawCount_table[key]), file=RAW_OUT)

    def mismatchProfile(self, read_length: int, read_num: int, outfile: str, q_cut: int = 30) -> None:
        """
        Calculate mismatch profile. Note that the "MD" tag must exist.
        """

        with open(f"{outfile}.mismatch_profile.xls", "w") as DOUT, open(f"{outfile}.mismatch_profile.r", "w") as ROUT:
            # reading input SAM file
            if self.bam_format:
                print("Process BAM file ... ", end=" ", file=sys.stderr)
            else:
                print("Process SAM file ... ", end=" ", file=sys.stderr)

            MD_pat = _MD_PAT

            count = 0
            # data[read_coord][genotype] = geno_type_number
            data: dict[int, dict[str, int]] = collections.defaultdict(dict)
            for aligned_read in _pysam_iter(self.samfile):
                if count >= read_num:
                    print(f"Total reads used: {count}", file=sys.stderr)
                    break
                if not _passes_qc(aligned_read, q_cut):
                    continue

                # Skip if there is no mismatch, or there is deletion
                tags = aligned_read.tags
                skip = False
                for tag in tags:
                    if tag[0] == "NM" and tag[1] == 0:
                        skip = True  # skip reads with no mismatches
                    if (tag[0] == "MD") and ("^" in tag[1]):
                        skip = True  # skip as there is deletion from the reference
                if skip is True:
                    continue

                # skip partially mapped read
                read_seq = aligned_read.seq
                if len(read_seq) != read_length:
                    continue
                if "N" in read_seq:
                    continue

                matched_portion_size = 0
                for op, value in aligned_read.cigar:
                    if op == 0:
                        matched_portion_size += value
                if matched_portion_size != read_length:
                    continue

                count += 1
                is_reverse = aligned_read.is_reverse
                for tag in tags:
                    if tag[0] == "MD":
                        a = MD_pat.findall(tag[1])  # tag[1] = "5G19T75"; a = [('5', 'G'), ('19', 'T')]
                        read_coord = 0
                        for match_number, ref_base in a:
                            read_coord += int(match_number)
                            read_base = read_seq[read_coord]
                            if read_base == ref_base:
                                continue
                            genotype = f"{ref_base}2{read_base}"
                            idx = (read_length - read_coord - 1) if is_reverse else read_coord
                            if genotype not in data[idx]:
                                data[idx][genotype] = 1
                            else:
                                data[idx][genotype] += 1
                            read_coord += 1
            print(f"Total reads used: {count}", file=sys.stderr)

            if len(data) == 0:
                print("No mismatches found", file=sys.stderr)
                sys.exit(1)
            # write data out
            all_genotypes = ["A2C", "A2G", "A2T", "C2A", "C2G", "C2T", "G2A", "G2C", "G2T", "T2A", "T2C", "T2G"]
            gt_header = "\t".join(all_genotypes)
            print(f"read_pos\tsum\t{gt_header}", file=DOUT)
            for indx in sorted(data):
                tmp = [indx, sum(data[indx].values())]  # read position and sum of mismatches
                for i in all_genotypes:
                    if i in data[indx]:
                        tmp.append(data[indx][i])
                    else:
                        tmp.append(0)
                print("\t".join([str(i) for i in tmp]), file=DOUT)

            # write Rscript
            r_data = collections.defaultdict(list)
            for gt in all_genotypes:
                for indx in sorted(data):
                    if gt in data[indx]:
                        r_data[gt].append(data[indx][gt])
                    else:
                        r_data[gt].append(0)
            for k in sorted(r_data):
                print(f"{k}=c({','.join([str(i) for i in r_data[k]])})", file=ROUT)

            print(
                'color_code = c("green","powderblue",'
                '"lightseagreen","red","violetred4",'
                '"mediumorchid1","blue","royalblue",'
                '"steelblue1","orange","gold","black")',
                file=ROUT,
            )

            print(f"y_up_bound = max(c({','.join([f'log10({i}+1)' for i in all_genotypes])}))", file=ROUT)
            print(f"y_low_bound = min(c({','.join([f'log10({i}+1)' for i in all_genotypes])}))", file=ROUT)

            print(f'pdf("{outfile}.mismatch_profile.pdf")', file=ROUT)
            count = 1
            for gt in all_genotypes:
                if count == 1:
                    print(
                        f'plot(log10({gt}+1),type="l",'
                        f"col=color_code[{count}],"
                        f"ylim=c(y_low_bound,y_up_bound),"
                        f'ylab="log10(# of mismatch)",'
                        f"xlab=\"Read position (5'->3')\")",
                        file=ROUT,
                    )
                else:
                    print(f"lines(log10({gt}+1), col=color_code[{count}])", file=ROUT)
                count += 1
            gt_legend = ",".join([f'"{i}"' for i in all_genotypes])
            print(
                f"legend(13,y_up_bound,legend=c({gt_legend}), fill=color_code, border=color_code, ncol=4)",
                file=ROUT,
            )
            print("dev.off()", file=ROUT)

    def deletionProfile(self, read_length: int, read_num: int, outfile: str, q_cut: int = 30) -> None:
        """
        Calculate deletion profile.
        Deletion: Deletion from the read (relative to the reference), CIGAR operator 'D'
        """

        with open(f"{outfile}.deletion_profile.txt", "w") as DOUT, open(f"{outfile}.deletion_profile.r", "w") as ROUT:
            # reading input SAM file
            if self.bam_format:
                print("Process BAM file ... ", end=" ", file=sys.stderr)
            else:
                print("Process SAM file ... ", end=" ", file=sys.stderr)

            count = 0
            del_postns: dict[int, int] = collections.defaultdict(int)  # key: position of read. value: deletion times
            for aligned_read in _pysam_iter(self.samfile):
                if count >= read_num:
                    print(f"Total reads used: {count}", file=sys.stderr)
                    break
                if not _passes_qc(aligned_read, q_cut):
                    continue

                # skip if read doesn't have deletion
                read_cigar = aligned_read.cigar
                if not any(c == 2 for c, _s in read_cigar):
                    continue

                # skip partially mapped read
                read_seq = aligned_read.seq
                if len(read_seq) != read_length:
                    continue

                matched_portion_size = 0
                for op, value in aligned_read.cigar:
                    if op == 0:
                        matched_portion_size += value
                    if op == 4:
                        matched_portion_size += value
                    if op == 1:
                        matched_portion_size += value
                if matched_portion_size != read_length:
                    continue

                count += 1
                is_reverse = aligned_read.is_reverse
                del_positions = bam_cigar.fetch_deletion_range(read_cigar)
                for p, s in del_positions:
                    if is_reverse:
                        p = read_length - p
                    del_postns[p] += 1
            print(f"Total reads used: {count}", file=sys.stderr)

            del_count = []
            print("read_position\tdeletion_count", file=DOUT)
            for k in range(0, read_length):
                if k in del_postns:
                    print(f"{k}\t{del_postns[k]}", file=DOUT)
                    del_count.append(str(del_postns[k]))
                else:
                    print(f"{k}\t0", file=DOUT)
                    del_count.append("0")

            print(f'pdf("{outfile}.deletion_profile.pdf")', file=ROUT)
            print(f"pos=c({','.join([str(i) for i in range(0, read_length)])})", file=ROUT)
            print(f"value=c({','.join([i for i in del_count])})", file=ROUT)
            print(
                "plot(pos,value,type='b', col='blue',xlab=\"Read position (5'->3')\", ylab='Deletion count')", file=ROUT
            )
            print("dev.off()", file=ROUT)

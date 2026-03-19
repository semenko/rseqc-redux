"""
Analzye 10X genomics single cell BAM files.
"""

from __future__ import annotations

import collections
import csv
import glob
import logging
import os
import re
import subprocess
import sys
from collections.abc import Generator
from typing import Any

import pysam


def _pysam_iter(
    samfile: pysam.AlignmentFile | pysam.IteratorRow,
) -> Generator[pysam.AlignedSegment, None, None]:
    """Iterate over pysam AlignmentFile or fetch iterator, handling ValueError on Python 3.13+."""
    try:
        yield from samfile
    except ValueError:
        return


def _write_edits_csv(mat: dict[int, dict[str, int]], outfile: str) -> None:
    """Write a nucleotide editing matrix (dict of dicts) to CSV.

    Input *mat* maps position (int) -> {edit_type (str): count}.
    Output CSV has rows = edit types (sorted), columns = positions (sorted),
    matching the original pandas ``DataFrame.from_dict(mat)`` layout where
    outer keys become columns and inner keys become rows.
    """
    # Columns = outer keys (positions)
    col_keys = sorted(mat.keys())
    # Rows = union of inner keys (edit types)
    row_keys_set: set[str] = set()
    for d in mat.values():
        row_keys_set.update(d.keys())
    row_keys = sorted(row_keys_set)

    with open(outfile, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Index"] + [str(c) for c in col_keys])
        for rk in row_keys:
            writer.writerow([rk] + [mat[ck].get(rk, 0) for ck in col_keys])


def diff_str(s1: str, s2: str) -> list[list[Any]]:
    """
    Comparing orignal barcode to the corrected barcode
    find the index and the nucleotide that has been corrected.


    Parameters
    ----------
    s1 : str
            the original barcode
    s2 : str
            the corrrected barcode

    """
    results: list[list[Any]] = []
    if len(s1) != len(s2):
        return results
    diff_positions = [i for i in range(len(s1)) if s1[i] != s2[i]]
    for pos in diff_positions:
        results.append([pos, s1[pos], s2[pos]])
    return results


def read_match_type(cigar_str: str) -> str:
    """return the matching type between read and ref"""
    match_type = ""
    if bool(re.search(r"\A\d+M\Z", cigar_str)):
        match_type = "Map_consecutively"
    elif bool(re.search(r"\A\d+M\d+N\d+M\Z", cigar_str)):
        match_type = "Map_with_splicing"
    elif bool(re.search(r"\A\d+S\d+M\Z", cigar_str)):
        match_type = "Map_with_clipping"
    elif bool(re.search(r"\A\d+M\d+S\Z", cigar_str)):
        match_type = "Map_with_clipping"
    elif bool(re.search(r"\A\d+M\d+N\d+M\d+S\Z", cigar_str)):
        match_type = "Map_with_splicing_and_clipping"
    elif bool(re.search(r"\A\d+S\d+M\d+N\d+M\Z", cigar_str)):
        match_type = "Map_with_splicing_and_clipping"
    else:
        match_type = "Others"
    return match_type


def list2str(lst: list[tuple[int, int]]) -> str:
    """
    translate samtools returned cigar_list into cigar_string
    """
    code2Char = {"0": "M", "1": "I", "2": "D", "3": "N", "4": "S", "5": "H", "6": "P", "7": "=", "8": "X"}
    cigar_str = ""
    for i in lst:
        cigar_str += str(i[1]) + code2Char[str(i[0])]
    return cigar_str


def barcode_edits(
    infile: str,
    outfile: str,
    step_size: int = 10000,
    limit: int = 2000000,
    CR_tag: str = "CR",
    CB_tag: str = "CB",
    UR_tag: str = "UR",
    UB_tag: str = "UB",
) -> None:
    """
    Analzye barcode in BAM file.

    Parameters
    ----------
    infile : str
            Input BAM file. Must be sorted and indexed.
    outfile : str
            Prefix of output files.
    step_size: int
            Output progress report when step_size alignments have been processed.
    limit : int
            Only process this number of alignments and stop.
    """
    logging.info('Reading BAM file "%s" ...' % infile)
    samfile = pysam.AlignmentFile(infile, "rb")

    CB_miss = 0  # number of reads without cell barcode
    CB_same = 0  # number of reads whose original cell barcode same as edited barcode
    CB_diff = 0  # number of reads whose cell barcode has been edited
    CB_freq: dict[str, int] = collections.defaultdict(int)  # cell barcode: raw reads
    CB_corrected_bases: dict[int, dict[str, int]] = collections.defaultdict(lambda: collections.defaultdict(int))

    UMI_miss = 0
    UMI_same = 0
    UMI_diff = 0
    UMI_freq: dict[str, int] = collections.defaultdict(int)  # UMI : raw reads
    UMI_corrected_bases: dict[int, dict[str, int]] = collections.defaultdict(lambda: collections.defaultdict(int))

    total_alignments = 0
    try:
        for aligned_read in _pysam_iter(samfile):
            total_alignments += 1
            tag_dict = dict(aligned_read.tags)  # type: ignore[attr-defined]  # {'NM': 1, 'RG': 'L1'}

            original_CB = ""
            corrected_CB = ""
            if CR_tag in tag_dict and CB_tag in tag_dict:
                original_CB = tag_dict[CR_tag].replace("-1", "")
                corrected_CB = tag_dict[CB_tag].replace("-1", "")
                CB_freq[corrected_CB] += 1
                if original_CB != corrected_CB:
                    CB_diff += 1
                    for diff in diff_str(original_CB, corrected_CB):
                        CB_corrected_bases[diff[0]][diff[1] + ":" + diff[2]] += 1
                else:
                    CB_same += 1
            else:
                CB_miss += 1

            original_UMI = ""
            corrected_UMI = ""
            if UR_tag in tag_dict and UB_tag in tag_dict:
                original_UMI = tag_dict[UR_tag].replace("-1", "")
                corrected_UMI = tag_dict[UB_tag].replace("-1", "")
                UMI_freq[corrected_UMI] += 1
                if original_UMI != corrected_UMI:
                    UMI_diff += 1
                    for diff in diff_str(original_UMI, corrected_UMI):
                        UMI_corrected_bases[diff[0]][diff[1] + ":" + diff[2]] += 1
                else:
                    UMI_same += 1
            else:
                UMI_miss += 1

            if total_alignments % step_size == 0:
                print("%d alignments processed.\r" % total_alignments, end=" ", file=sys.stderr)
            if limit is not None:
                if total_alignments >= limit:
                    break
    finally:
        samfile.close()
    logging.info("Total alignments processed: %d" % total_alignments)

    logging.info("Number of alignmenets with <cell barcode> kept AS IS: %d" % CB_same)
    logging.info("Number of alignmenets with <cell barcode> edited: %d" % CB_diff)
    logging.info("Number of alignmenets with <cell barcode> missing: %d" % CB_miss)
    logging.info("Number of alignmenets with UMI kept AS IS: %d" % UMI_same)
    logging.info("Number of alignmenets with UMI edited: %d" % UMI_diff)
    logging.info("Number of alignmenets with UMI missing: %d" % UMI_miss)

    # writing cell barcode
    logging.info('Writing cell barcode frequencies to "%s"' % (outfile + ".CB_freq.tsv"))
    with open(outfile + ".CB_freq.tsv", "w") as CB_OUT:
        for bc, count in sorted(CB_freq.items(), key=lambda item: item[1], reverse=True):
            CB_OUT.write(bc + "\t" + str(count) + "\n")

    # writing UMI
    logging.info('Writing UMI frequencies to "%s"' % (outfile + ".UMI_freq.tsv"))
    with open(outfile + ".UMI_freq.tsv", "w") as UMI_OUT:
        for bc, count in sorted(UMI_freq.items(), key=lambda item: item[1], reverse=True):
            UMI_OUT.write(bc + "\t" + str(count) + "\n")

    CB_mat_file = outfile + ".CB_edits_count.csv"
    logging.info('Writing the nucleotide editing matrix (count) of cell barcode to "%s"' % CB_mat_file)
    _write_edits_csv(CB_corrected_bases, CB_mat_file)

    UMI_mat_file = outfile + ".UMI_edits_count.csv"
    logging.info('Writing the nucleotide editing matrix of molecular barcode (UMI) to "%s"' % UMI_mat_file)
    _write_edits_csv(UMI_corrected_bases, UMI_mat_file)


def mapping_stat(
    infile: str,
    step_size: int = 50000,
    CB_tag: str = "CB",
    UMI_tag: str = "UB",
    RE_tag: str = "RE",
    TX_tag: str = "TX",
    AN_tag: str = "AN",
    xf_tag: str = "xf",
    chrM_id: str = "chrM",
    n_thread: int = 1,
) -> None:
    """
    Reads mapping statistics

    Parameters
    ----------
    infile : str
            Input BAM file. Must be sorted and indexed.
    outfile : str
            Prefix of output files. If outfile is None, only do counting and do not generate BAM files.
    step_size: int
            Output progress report when step_size alignments have been processed.
    limit : int
            Only process this number of alignments and stop.
    """
    logging.info('Reading BAM file "%s" ...' % infile)
    try:
        # older versions of pysam
        samfile = pysam.AlignmentFile(infile, mode="rb", require_index=True, thread=n_thread)  # type: ignore[call-arg]
    except (TypeError, ValueError):
        # latest verion of pysam (v0.19.1)
        samfile = pysam.AlignmentFile(infile, mode="rb", require_index=True, threads=n_thread)
    if not samfile.check_index():
        logging.error("Cannot find the index file")
        samfile.close()
        sys.exit(1)
    chrom_info = list(zip(samfile.references, samfile.lengths))  # [('chr1', 195471971), ('chr10', 130694993),...]

    total_alignments = 0
    confi_alignments = 0

    total_reads_n = 0
    confi_reads_n = 0

    ## for confidently mapped reads
    # PCR duplicate or not
    confi_reads_dup = 0
    confi_reads_nondup = 0

    # Reverse or forward
    confi_reads_rev = 0
    confi_reads_fwd = 0

    # with error-corrected CB or UMI
    confi_CB = 0
    confi_UB = 0

    exon_reads = 0
    intron_reads = 0
    intergenic_reads = 0
    other_reads1 = 0

    # sense or antisense
    sense_reads = 0
    anti_reads = 0
    other_reads2 = 0

    chrM_reads = 0
    # read match type
    read_type: dict[str, int] = collections.defaultdict(int)

    try:
        for chr_id, chr_len in chrom_info:
            logging.info('Processing "%s" ...' % chr_id)
            chrom_count = 0
            chrom_total_reads = set()  # total reads in BAM file
            chrom_confi_reads = set()  # reads marked as confidently mapped to transcriptome by xf:i:1 tag

            with open(chr_id + ".all_reads_id.txt", "w") as ALL, open(chr_id + ".confident_reads_id.txt", "w") as CONF:
                for aligned_read in _pysam_iter(samfile.fetch(chr_id)):
                    total_alignments += 1
                    chrom_count += 1
                    read_id = aligned_read.query_name
                    tag_dict = dict(aligned_read.tags)  # type: ignore[attr-defined]  # {'NM': 1, 'RG': 'L1'}
                    cigar_str = list2str(aligned_read.cigar)  # type: ignore[attr-defined]
                    chrom_total_reads.add(read_id)

                    # confident alignments
                    if xf_tag in tag_dict and tag_dict[xf_tag] & 0x1 != 0:
                        if chr_id == chrM_id:
                            chrM_reads += 1
                        # with or without CB/UMI barcode
                        if CB_tag in tag_dict:
                            confi_CB += 1
                        if UMI_tag in tag_dict:
                            confi_UB += 1

                        # duplicate or not
                        if aligned_read.is_duplicate:
                            confi_reads_dup += 1
                        else:
                            confi_reads_nondup += 1

                        # forward or reverse
                        if aligned_read.is_reverse:
                            confi_reads_rev += 1
                        else:
                            confi_reads_fwd += 1

                        # Single character indicating the region type of this alignment
                        # (E = exonic, N = intronic, I = intergenic).
                        if RE_tag in tag_dict:
                            if tag_dict[RE_tag] == "E":
                                exon_reads += 1
                            elif tag_dict[RE_tag] == "I":
                                intron_reads += 1
                            elif tag_dict[RE_tag] == "N":
                                intergenic_reads += 1
                            else:
                                other_reads1 += 1
                        else:
                            other_reads1 += 1

                        # sense or antisense
                        if TX_tag in tag_dict:
                            sense_reads += 1
                        elif AN_tag in tag_dict:
                            anti_reads += 1
                        else:
                            other_reads2 += 1

                        # map type
                        cigar_str = list2str(aligned_read.cigar)  # type: ignore[attr-defined]
                        tmp = read_match_type(cigar_str)
                        read_type[tmp] += 1

                        confi_alignments += 1
                        chrom_confi_reads.add(read_id)

                    if chrom_count % step_size == 0:
                        print("%d alignments processed.\r" % chrom_count, end=" ", file=sys.stderr)
                logging.info('Processed %d alignments from "%s"' % (chrom_count, chr_id))
                for i in chrom_total_reads:
                    print(i, file=ALL)
                for i in chrom_confi_reads:
                    print(i, file=CONF)
    finally:
        samfile.close()

    logging.info("Processing total %d alignments mapped to all chromosomes." % total_alignments)

    logging.info("Count total mapped reads ...")
    output1 = subprocess.check_output("awk '!a[$0]++' *.all_reads_id.txt |tee All_reads_uniqID.txt | wc -l", shell=True)

    logging.info("Count confidently mapped reads ...")
    output2 = subprocess.check_output(
        "awk '!a[$0]++' *.confident_reads_id.txt |tee confident_reads_uniqID.txt | wc -l", shell=True
    )

    logging.info("Removing intermediate files ...")
    for f in glob.glob("*.all_reads_id.txt"):
        os.unlink(f)
    for f in glob.glob("*.confident_reads_id.txt"):
        os.unlink(f)

    total_reads_n = int(output1.decode("utf-8").strip())
    confi_reads_n = int(output2.decode("utf-8").strip())
    non_confi_reads = total_reads_n - confi_reads_n

    print("")
    print("\nTotal_alignments: %d" % total_alignments)
    print("└--Confident_alignments: %d" % confi_alignments)
    print("")
    print("Total_mapped_reads:\t%d" % total_reads_n)
    print("|--Non_confidently_mapped_reads:\t%d\t(%.2f%%)" % (non_confi_reads, non_confi_reads * 100.0 / total_reads_n))
    print("└--Confidently_mapped_reads:\t%d\t(%.2f%%)" % (confi_reads_n, confi_reads_n * 100.0 / total_reads_n))

    print("   |--Reads_with_PCR_duplicates:\t%d\t(%.2f%%)" % (confi_reads_dup, confi_reads_dup * 100.0 / confi_reads_n))
    print(
        "   └--Reads_no_PCR_duplicates:\t%d\t(%.2f%%)"
        % (confi_reads_nondup, confi_reads_nondup * 100.0 / confi_reads_n)
    )
    print("")

    print(
        "   |--Reads_map_to_forward(Waston)_strand:\t%d\t(%.2f%%)"
        % (confi_reads_fwd, confi_reads_fwd * 100.0 / confi_reads_n)
    )
    print(
        "   └--Reads_map_to_Reverse(Crick)_strand:\t%d\t(%.2f%%)"
        % (confi_reads_rev, confi_reads_rev * 100.0 / confi_reads_n)
    )
    print("")

    print("   |--Reads_map_to_sense_strand:\t%d\t(%.2f%%)" % (sense_reads, sense_reads * 100.0 / confi_reads_n))
    print("   └--Reads_map_to_antisense_strand:\t%d\t(%.2f%%)" % (anti_reads, anti_reads * 100.0 / confi_reads_n))
    print("   └--Other:\t%d\t(%.2f%%)" % (other_reads2, other_reads2 * 100.0 / confi_reads_n))
    print("")

    print("   |--Reads_map_to_exons:\t%d\t(%.2f%%)" % (exon_reads, exon_reads * 100.0 / confi_reads_n))
    print("   └--Reads_map_to_introns:\t%d\t(%.2f%%)" % (intron_reads, intron_reads * 100.0 / confi_reads_n))
    print("   └--Reads_map_to_intergenic:\t%d\t(%.2f%%)" % (intergenic_reads, intergenic_reads * 100.0 / confi_reads_n))
    print("   └--Other:\t%d\t(%.2f%%)" % (other_reads1, other_reads1 * 100.0 / confi_reads_n))
    print("")

    print("   |--Reads_with_error-corrected_barcode:\t%d\t(%.2f%%)" % (confi_CB, confi_CB * 100.0 / confi_reads_n))
    print(
        "   └--Reads_no_error-corrected_barcode:\t%d\t(%.2f%%)"
        % ((confi_reads_n - confi_CB), (confi_reads_n - confi_CB) * 100.0 / confi_reads_n)
    )
    print("")
    print("   |--Reads_with_error-corrected_UMI:\t%d\t(%.2f%%)" % (confi_UB, confi_UB * 100.0 / confi_reads_n))
    print(
        "   └--Reads_no_error-corrected_UMI:\t%d\t(%.2f%%)"
        % ((confi_reads_n - confi_UB), (confi_reads_n - confi_UB) * 100.0 / confi_reads_n)
    )
    print("")
    print("   |--Reads_map_to_mitochonrial_genome:\t%d\t(%.2f%%)" % (chrM_reads, chrM_reads * 100.0 / confi_reads_n))
    print(
        "   └--Reads_map_to_nuclear_genome:\t%d\t(%.2f%%)"
        % ((confi_reads_n - chrM_reads), (confi_reads_n - chrM_reads) * 100.0 / confi_reads_n)
    )

    print("")

    for i in sorted(read_type):
        if i == "Others":
            continue
        print("   |--%s:\t%d\t(%.2f%%)" % (i, read_type[i], read_type[i] * 100.0 / confi_reads_n))
    print("   └--%s:\t%d\t(%.2f%%)" % ("Others", read_type["Others"], read_type["Others"] * 100.0 / confi_reads_n))
    print("")

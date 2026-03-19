"""
Analzye 10X genomics single cell BAM files.
"""

from __future__ import annotations

import collections
import csv
import logging
import re
import sys
from typing import Any

import pysam

from rseqc.cli_common import _pysam_iter  # noqa: F401

# Pre-compiled CIGAR classification patterns for read_match_type()
_RE_CONSECUTIVE = re.compile(r"\A\d+M\Z")
_RE_SPLICING = re.compile(r"\A\d+M\d+N\d+M\Z")
_RE_CLIP_LEFT = re.compile(r"\A\d+S\d+M\Z")
_RE_CLIP_RIGHT = re.compile(r"\A\d+M\d+S\Z")
_RE_SPLICE_CLIP_RIGHT = re.compile(r"\A\d+M\d+N\d+M\d+S\Z")
_RE_SPLICE_CLIP_LEFT = re.compile(r"\A\d+S\d+M\d+N\d+M\Z")

# CIGAR operation code → character (indexed by int op code)
_CIGAR_CHAR = ("M", "I", "D", "N", "S", "H", "P", "=", "X")


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
    if _RE_CONSECUTIVE.search(cigar_str):
        return "Map_consecutively"
    if _RE_SPLICING.search(cigar_str):
        return "Map_with_splicing"
    if _RE_CLIP_LEFT.search(cigar_str) or _RE_CLIP_RIGHT.search(cigar_str):
        return "Map_with_clipping"
    if _RE_SPLICE_CLIP_RIGHT.search(cigar_str) or _RE_SPLICE_CLIP_LEFT.search(cigar_str):
        return "Map_with_splicing_and_clipping"
    return "Others"


def list2str(lst: list[tuple[int, int]]) -> str:
    """Translate pysam cigar_list into a CIGAR string."""
    return "".join(f"{size}{_CIGAR_CHAR[op]}" for op, size in lst)


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
    logging.info(f'Reading BAM file "{infile}" ...')
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
                print(f"{total_alignments} alignments processed.\r", end=" ", file=sys.stderr)
            if limit is not None:
                if total_alignments >= limit:
                    break
    finally:
        samfile.close()
    logging.info(f"Total alignments processed: {total_alignments}")

    logging.info(f"Number of alignmenets with <cell barcode> kept AS IS: {CB_same}")
    logging.info(f"Number of alignmenets with <cell barcode> edited: {CB_diff}")
    logging.info(f"Number of alignmenets with <cell barcode> missing: {CB_miss}")
    logging.info(f"Number of alignmenets with UMI kept AS IS: {UMI_same}")
    logging.info(f"Number of alignmenets with UMI edited: {UMI_diff}")
    logging.info(f"Number of alignmenets with UMI missing: {UMI_miss}")

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
    logging.info(f'Writing the nucleotide editing matrix (count) of cell barcode to "{CB_mat_file}"')
    _write_edits_csv(CB_corrected_bases, CB_mat_file)

    UMI_mat_file = outfile + ".UMI_edits_count.csv"
    logging.info(f'Writing the nucleotide editing matrix of molecular barcode (UMI) to "{UMI_mat_file}"')
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
    logging.info(f'Reading BAM file "{infile}" ...')
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

    all_read_ids: set[str] = set()
    confi_read_ids: set[str] = set()

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
            logging.info(f'Processing "{chr_id}" ...')
            chrom_count = 0

            for aligned_read in _pysam_iter(samfile.fetch(chr_id)):
                total_alignments += 1
                chrom_count += 1
                read_id = aligned_read.query_name
                tag_dict = dict(aligned_read.tags)  # type: ignore[attr-defined]  # {'NM': 1, 'RG': 'L1'}
                all_read_ids.add(read_id)

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
                    confi_read_ids.add(read_id)

                if chrom_count % step_size == 0:
                    print(f"{chrom_count} alignments processed.\r", end=" ", file=sys.stderr)
            logging.info(f'Processed {chrom_count} alignments from "{chr_id}"')
    finally:
        samfile.close()

    logging.info(f"Processing total {total_alignments} alignments mapped to all chromosomes.")

    total_reads_n = len(all_read_ids)
    confi_reads_n = len(confi_read_ids)
    non_confi_reads = total_reads_n - confi_reads_n

    print("")
    print(f"\nTotal_alignments: {total_alignments}")
    print(f"└--Confident_alignments: {confi_alignments}")
    print("")
    print(f"Total_mapped_reads:\t{total_reads_n}")
    non_confi_pct = non_confi_reads * 100.0 / total_reads_n
    print(f"|--Non_confidently_mapped_reads:\t{non_confi_reads}\t({non_confi_pct:.2f}%)")
    confi_pct = confi_reads_n * 100.0 / total_reads_n
    print(f"└--Confidently_mapped_reads:\t{confi_reads_n}\t({confi_pct:.2f}%)")

    dup_pct = confi_reads_dup * 100.0 / confi_reads_n
    print(f"   |--Reads_with_PCR_duplicates:\t{confi_reads_dup}\t({dup_pct:.2f}%)")
    nondup_pct = confi_reads_nondup * 100.0 / confi_reads_n
    print(f"   └--Reads_no_PCR_duplicates:\t{confi_reads_nondup}\t({nondup_pct:.2f}%)")
    print("")

    fwd_pct = confi_reads_fwd * 100.0 / confi_reads_n
    print(f"   |--Reads_map_to_forward(Waston)_strand:\t{confi_reads_fwd}\t({fwd_pct:.2f}%)")
    rev_pct = confi_reads_rev * 100.0 / confi_reads_n
    print(f"   └--Reads_map_to_Reverse(Crick)_strand:\t{confi_reads_rev}\t({rev_pct:.2f}%)")
    print("")

    sense_pct = sense_reads * 100.0 / confi_reads_n
    print(f"   |--Reads_map_to_sense_strand:\t{sense_reads}\t({sense_pct:.2f}%)")
    anti_pct = anti_reads * 100.0 / confi_reads_n
    print(f"   └--Reads_map_to_antisense_strand:\t{anti_reads}\t({anti_pct:.2f}%)")
    other2_pct = other_reads2 * 100.0 / confi_reads_n
    print(f"   └--Other:\t{other_reads2}\t({other2_pct:.2f}%)")
    print("")

    exon_pct = exon_reads * 100.0 / confi_reads_n
    print(f"   |--Reads_map_to_exons:\t{exon_reads}\t({exon_pct:.2f}%)")
    intron_pct = intron_reads * 100.0 / confi_reads_n
    print(f"   └--Reads_map_to_introns:\t{intron_reads}\t({intron_pct:.2f}%)")
    intergenic_pct = intergenic_reads * 100.0 / confi_reads_n
    print(f"   └--Reads_map_to_intergenic:\t{intergenic_reads}\t({intergenic_pct:.2f}%)")
    other1_pct = other_reads1 * 100.0 / confi_reads_n
    print(f"   └--Other:\t{other_reads1}\t({other1_pct:.2f}%)")
    print("")

    cb_pct = confi_CB * 100.0 / confi_reads_n
    print(f"   |--Reads_with_error-corrected_barcode:\t{confi_CB}\t({cb_pct:.2f}%)")
    no_cb = confi_reads_n - confi_CB
    no_cb_pct = no_cb * 100.0 / confi_reads_n
    print(f"   └--Reads_no_error-corrected_barcode:\t{no_cb}\t({no_cb_pct:.2f}%)")
    print("")
    ub_pct = confi_UB * 100.0 / confi_reads_n
    print(f"   |--Reads_with_error-corrected_UMI:\t{confi_UB}\t({ub_pct:.2f}%)")
    no_ub = confi_reads_n - confi_UB
    no_ub_pct = no_ub * 100.0 / confi_reads_n
    print(f"   └--Reads_no_error-corrected_UMI:\t{no_ub}\t({no_ub_pct:.2f}%)")
    print("")
    chrm_pct = chrM_reads * 100.0 / confi_reads_n
    print(f"   |--Reads_map_to_mitochonrial_genome:\t{chrM_reads}\t({chrm_pct:.2f}%)")
    nuclear = confi_reads_n - chrM_reads
    nuclear_pct = nuclear * 100.0 / confi_reads_n
    print(f"   └--Reads_map_to_nuclear_genome:\t{nuclear}\t({nuclear_pct:.2f}%)")

    print("")

    for i in sorted(read_type):
        if i == "Others":
            continue
        rt_pct = read_type[i] * 100.0 / confi_reads_n
        print(f"   |--{i}:\t{read_type[i]}\t({rt_pct:.2f}%)")
    others_pct = read_type["Others"] * 100.0 / confi_reads_n
    print(f"   └--Others:\t{read_type['Others']}\t({others_pct:.2f}%)")
    print("")

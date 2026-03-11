#!/usr/bin/env python
"""deal with Kmer. DNA sequence should only A, C, G, T."""

import itertools
import math
import re
from collections import Counter
from collections.abc import Generator


def word_generator(seq: str, word_size: int, step_size: int, frame: int = 0) -> Generator[str, None, None]:
    """generate DNA word from sequence using word_size and step_size. Frame is 0, 1 or2"""
    for i in range(frame, len(seq), step_size):
        word = seq[i : i + word_size]
        if len(word) == word_size:
            yield word


def seq_generator(fastafile: str) -> Generator[list[str], None, None]:
    """DNA sequence only contains A,C,G,T,N. sequence with other characters will be removed"""
    tmpseq = ""
    name = ""
    DNA_pat = re.compile(r"^[ACGTN]+$")
    for line in open(fastafile, "r"):
        line = line.strip().upper()
        if line.startswith(("#", " ", "\n")):
            continue
        if line.startswith((">", "@")):
            if tmpseq:
                yield [name, tmpseq]
                tmpseq = ""
            name = line.split()[0][1:]
        elif DNA_pat.match(line):
            tmpseq += line
    yield [name, tmpseq]


def all_possible_kmer(l: int) -> Generator[str, None, None]:
    """return all possible combinations of A,C,G,T,N. only support A,C,G,T,N. l is length of kmer"""
    for i in itertools.product(["A", "C", "G", "T", "N"], repeat=l):
        yield "".join(i)


def kmer_freq_file(
    fastafile: str, word_size: int, step_size: int = 1, frame: int = 0, min_count: int = 0
) -> dict[str, int]:
    """Calculate kmer frequency from fasta file"""
    seq_num = 0
    ret_dict = {}
    for n, s in seq_generator(fastafile):
        seq_num += 1
        if seq_num == 1:
            count_table = Counter(word_generator(s, word_size=word_size, step_size=step_size, frame=frame))
        else:
            count_table.update(word_generator(s, word_size=word_size, step_size=step_size, frame=frame))

    for kmer in all_possible_kmer(word_size):
        if kmer not in count_table:
            count_table[kmer] = 0
        if count_table[kmer] >= min_count:
            if "N" in kmer:
                continue
            ret_dict[kmer] = count_table[kmer]
    return ret_dict


def kmer_freq_seq(seq: str, word_size: int, step_size: int = 1, frame: int = 0, min_count: int = 0) -> None:
    """Calculate kmer frequency from DNA sequence. coding. genome is hexamer table calculated
    from coding region and whole genome (as background control)
    """
    count_table = Counter(word_generator(seq, word_size=word_size, step_size=step_size, frame=frame))
    for kmer in all_possible_kmer(word_size):
        if kmer not in count_table:
            count_table[kmer] = 0
        if count_table[kmer] >= min_count:
            print(kmer + "\t" + str(count_table[kmer]))


def kmer_ratio(seq: str, word_size: int, step_size: int, coding: dict[str, int], noncoding: dict[str, int]) -> float:
    if len(seq) < word_size:
        return 0

    sum_of_log_ratio_0 = 0.0
    frame0_count = 0.0
    for k in word_generator(seq=seq, word_size=word_size, step_size=step_size, frame=0):
        if (k not in coding) or (k not in noncoding):
            continue
        if coding[k] > 0 and noncoding[k] > 0:
            sum_of_log_ratio_0 += math.log(coding[k] / noncoding[k])
        elif coding[k] > 0 and noncoding[k] == 0:
            sum_of_log_ratio_0 += 1
        elif coding[k] == 0 and noncoding[k] == 0:
            continue
        elif coding[k] == 0 and noncoding[k] > 0:
            sum_of_log_ratio_0 -= 1
        else:
            continue
        frame0_count += 1
    return sum_of_log_ratio_0 / frame0_count

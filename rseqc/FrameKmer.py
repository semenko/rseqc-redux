#!/usr/bin/env python
"""deal with Kmer. DNA sequence should only A, C, G, T."""

import itertools
import re
from collections import Counter
from collections.abc import Generator

_DNA_PAT = re.compile(r"^[ACGTN]+$")


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
    with open(fastafile, "r") as _fh:
        for line in _fh:
            line = line.strip().upper()
            if line.startswith(("#", " ", "\n")):
                continue
            if line.startswith((">", "@")):
                if tmpseq:
                    yield [name, tmpseq]
                    tmpseq = ""
                name = line.split()[0][1:]
            elif _DNA_PAT.match(line):
                tmpseq += line
    yield [name, tmpseq]


def all_possible_kmer(length: int) -> Generator[str, None, None]:
    """return all possible combinations of A,C,G,T,N. only support A,C,G,T,N. length is length of kmer"""
    for i in itertools.product(["A", "C", "G", "T", "N"], repeat=length):
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

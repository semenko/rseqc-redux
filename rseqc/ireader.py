"""
read compressed (.gz .bz) files
"""

import bz2
import gzip
import sys
import urllib
from collections.abc import Generator
from subprocess import PIPE, Popen
from typing import IO, Any


def nopen(f: str | IO, mode: str = "rb") -> Any:
    if not isinstance(f, str):
        return f
    if f.startswith("|"):
        p = Popen(f[1:], stdout=PIPE, stdin=PIPE, shell=True)
        if mode[0] == "r":
            return p.stdout
        return p
    return (
        {"r": sys.stdin, "w": sys.stdout}[mode[0]]
        if f == "-"
        else gzip.open(f, mode)
        if f.endswith((".gz", ".Z", ".z"))
        else bz2.BZ2File(f, mode)  # type: ignore[call-overload]
        if f.endswith((".bz", ".bz2", ".bzip2"))
        else urllib.urlopen(f)  # type: ignore[attr-defined]  # known bug: should be urllib.request.urlopen
        if f.startswith(("http://", "https://", "ftp://"))
        else open(f, mode)
    )


def reader(fname: str) -> Generator[str]:
    for l in nopen(fname):
        yield l.decode("utf8").strip().replace("\r", "")

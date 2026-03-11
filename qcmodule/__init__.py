"""Backward-compatibility shim: ``from qcmodule import SAM`` keeps working."""

import importlib as _importlib

_SUBMODULES = {
    "SAM", "BED", "bam_cigar", "scbam", "fasta", "fastq",
    "heatmap", "mystat", "getBamFiles", "twoList", "FrameKmer",
}


def __getattr__(name: str):
    if name in _SUBMODULES:
        return _importlib.import_module(f"rseqc.{name}")
    raise AttributeError(f"module 'qcmodule' has no attribute {name!r}")

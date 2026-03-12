# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Removed

- Removed `ParseSAM` and `QCSAM` dead code classes from `rseqc/SAM.py` (~2,970 lines). All scripts use `ParseBAM` via pysam.
- Removed `import string`, `import os`, and `from bx_extras.fpconst import isNaN` from `rseqc/SAM.py` (only used by dead classes).
- Removed `qcmodule/` backward-compatibility shim — all scripts now import directly from `rseqc`.
- Removed Python 3 version checks from all 26 scripts (dead code on Python >=3.10).
- Removed dead `open(outfile, "w")` call in `annoGene.py` `annotateBed()`.
- Removed unnecessary `list()` wrappers on dict view iteration in `fasta.py` (10 sites) and `orf.py` (2 sites).

### Fixed

- Fixed `cdsEdn_float` typo in `BED.py` `unionBED()` — was `NameError` when `boundary="cds"` (Bug #1).
- Fixed `name` before assignment in `fasta.py` `Fasta.__init__()` — was `UnboundLocalError` on empty/malformed FASTA (Bug #2).
- Fixed `if int(fields[9] == 1)` operator precedence in `BED.py` (6 sites) — single-exon gene guard now works correctly (Bug #3).
- Fixed `sc`/`tc` parameter shadowing `start_coden`/`stop_coden` in `orf.py` — was `UnboundLocalError` when calling without explicit codons (Bug #4).
- Fixed `chromm` typo in `annoGene.py` `getUTRExonFromLine()` and `getCDSExonFromLine()` — was `NameError` (Bug #5).
- Fixed `txstart`/`txEnd` undefined and `tmp.append()` call in `annoGene.py` `getExonFromFile2()` — was `NameError` (Bug #6).
- Fixed `input_file` undefined in `bam_stat.py` — should be `options.input_file` (Bug #7).
- Fixed `urllib.urlopen` in `ireader.py` — replaced with `urllib.request.urlopen` (Bug #8).
- Fixed `subtractBed3` no-op guard in `BED.py` — dead `if chrom not in bitsets1` check removed (Bug #9).
- Fixed `Hill_number` q=1 in `mystat.py` — now correctly parses comma-separated string before passing to `shannon_entropy()` (Bug #11).
- Fixed `sys.exit()` → `sys.exit(1)` in 6 library error paths (`SAM.py`, `BED.py`) so errors exit non-zero.
- Fixed Python 3.13 compatibility: catch `ValueError` alongside `StopIteration` when iterating pysam alignments (pysam bug with PEP 745 changes). 29 sites across 7 files.

### Changed

- Migrated all 33 CLI scripts from `optparse` to `argparse`.
- Renamed ambiguous variable `l` to descriptive names across 12 files (E741).
- `--help` output formatting now uses argparse style.
- `class ParseBAM(object):` → `class ParseBAM:` (modern Python 3 style).
- `(v1 + 1).__div__(v2 + 1)` → `(v1 + 1) / (v2 + 1)` in `twoList.py` (Python 2 artifact).
- Improved type annotations in `changePoint.py` (removed 2 `type: ignore` suppressions).

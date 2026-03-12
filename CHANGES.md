# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Removed

- Removed `ParseSAM` and `QCSAM` dead code classes from `rseqc/SAM.py` (~2,970 lines). All scripts use `ParseBAM` via pysam.
- Removed `import string`, `import os`, and `from bx_extras.fpconst import isNaN` from `rseqc/SAM.py` (only used by dead classes).
- Removed `qcmodule/` backward-compatibility shim — all scripts now import directly from `rseqc`.
- Removed Python 3 version checks from all 26 scripts (dead code on Python >=3.10).
- Removed dead `open(outfile, "w")` call in `annoGene.py` `annotateBed()`.
- Removed unnecessary `list()` wrappers on dict view iteration in `fasta.py` (10 sites), `orf.py` (2 sites), `SAM.py` (5 sites), `BED.py` (2 sites), `geneBody_coverage.py` (1 site).
- Removed dead `bootstrap()` function from `changePoint.py` (contained a list-to-float comparison bug, never called).
- Removed dead `outfile` parameter from `annoGene.py` `annotateBed()` (never used).
- Removed ~70 dead bare field expressions (`fields[N]`, `int(fields[N])`) across `BED.py`, `SAM.py`, `annoGene.py`, `wiggle.py`, and 4 scripts.

### Fixed

- Fixed `heatmap.py` string concatenation bug in `logging.error()` — `+` replaced with `%` format operator, was producing garbled error messages.
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
- Fixed `exit(0)`/`sys.exit(0)` → `sys.exit(1)` in ~90 error/validation paths across all scripts and library code.
- Fixed bare `exit()` → `sys.exit()` in `BED.py`, `SAM.py`, `geneBody_coverage.py`, `geneBody_coverage2.py`, `tin.py`.
- Fixed Python 3.13 compatibility: `_pysam_iter()` helper wraps all pysam BAM iteration to handle `ValueError` bug (PEP 745 / coverage instrumentation).

### Changed

- Replaced `while 1: next(samfile)` anti-pattern with `for` loops across `SAM.py` (20 sites), `scbam.py` (4 sites), and 5 scripts.
- Migrated all 33 CLI scripts from `optparse` to `argparse`.
- Renamed ambiguous variable `l` to descriptive names across 12 files (E741).
- `--help` output formatting now uses argparse style.
- `class ParseBAM(object):` → `class ParseBAM:` (modern Python 3 style).
- `(v1 + 1).__div__(v2 + 1)` → `(v1 + 1) / (v2 + 1)` in `twoList.py` (Python 2 artifact).
- Improved type annotations in `changePoint.py` (removed 2 `type: ignore` suppressions).
- E501 (line length > 120) fully resolved and enabled in ruff config — 0 violations remaining.
- Converted 95 bare `open()` calls to `with` statements across `rseqc/` and `scripts/` to prevent resource leaks.
- Replaced ~80 `list(map((lambda ...), ...))` with list comprehensions across all BED-parsing code.
- Added type hints to all 20 `rseqc/` library modules (SAM.py, BED.py, annoGene.py, fasta.py, fastq.py, scbam.py, wiggle.py, getBamFiles.py, heatmap.py + previously done modules).

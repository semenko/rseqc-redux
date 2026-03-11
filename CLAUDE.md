# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

rseqc-redux is a modernization of RSeQC 5.0.1 (RNA-seq Quality Control), originally by Liguo Wang. Licensed under GPLv3+. The original code has no tests, uses deprecated APIs (distutils, nose, psyco), and targets Python >=3.5.

## Goals (in priority order)

1. **Establish comprehensive test coverage** — write tests for all library modules and CLI scripts before changing behavior. Use pytest. Test with small, synthetic BAM/BED fixtures where possible.
2. **Identify and document bugs** — use tests to surface incorrect logic, off-by-one errors, and edge cases in the original code. File issues for each bug found.
3. **Modernize the codebase** — only after test coverage is solid:
   - Replace `setup.py`/distutils with `pyproject.toml` (using hatch or setuptools with modern config)
   - Use `uv` as the package manager for development workflows
   - Target Python >=3.10, drop legacy version checks
   - Remove dead code (`psyco_full.py`, `distribute_setup.py`)
   - Convert CLI scripts to proper entry points
   - Add type hints incrementally
4. **Improve performance** — only after extensive tests exist to verify correctness.

## Current State (as of 2026-03-11)

**Infrastructure:** Done — pyproject.toml, CI (3.10–3.13), PyPI publishing.

**Tests:** 274 passing, 1 xfail. Coverage: utility modules 78–97%, BED.py 13%, SAM.py 2%.

**Lint/Type:** CI green — ruff (0 errors), mypy (0 errors), ruff format clean.

**What's been modernized:**
- Star imports replaced with explicit imports in all `rseqc/` and `scripts/` files
- Bare `except:` → `except Exception:` across entire codebase (97 instances)
- Unused variables removed (147 instances via ruff F841)
- Type hints added to `cigar.py`, `ireader.py`, `bam_cigar.py`, `mystat.py`, `quantile.py`, `dotProduct.py`, `changePoint.py`, `twoList.py`, `FrameKmer.py`, `orf.py`
- Legacy boilerplate removed: Python 3 version checks (32 scripts), `__author__`/`__version__` metadata (44 files), UTF-8 encoding declarations, shebangs from library modules, "converted from python2.7" docstrings, dead code and commented-out prints
- CLI smoke tests for all 33 scripts (`--help` flag)

**What still needs work:**
- SAM.py/BED.py methods need integration tests with BAM fixtures (pysam-built)
- E501 (line length) and E741 (ambiguous vars) still suppressed

## Known Bugs (documented, NOT fixed — need regression tests first)

1. **`cdsEdn_float` typo** — `rseqc/BED.py:1370,1441` — `NameError` when `boundary == "cds"` in `unionBED()`. Should be `cdsEnd_float`. Marked with `# noqa: F821`.
2. **`name` before assignment** — `rseqc/fasta.py:48` — `UnboundLocalError` on empty/malformed FASTA input. Marked with `# noqa: F821`.
3. **`if int(fields[9] == 1)`** — `rseqc/BED.py:720` — operator precedence bug. `fields[9] == 1` evaluates to bool, `int(bool)` is 0 or 1, so the `continue` for single-exon genes never fires. Intron extraction works by accident because single-exon genes produce empty intron lists anyway.
4. **`sc` parameter shadows `start_coden`** — `rseqc/orf.py:63` — `UnboundLocalError` when calling `longest_orf()` without explicit `sc`/`tc` args. The `for sc in start_coden` loop variable shadows the function parameter, and Python detects `start_coden` as a local. Documented with xfail test.
5. **`chromm` typo** — `rseqc/annoGene.py:189,231` — `NameError` in `getUTRExonFromLine()` and `getCDSExonFromLine()`. Should be `chrom`. These functions are dead code (never called from scripts).
6. **`txstart`/`txEnd` undefined** — `rseqc/annoGene.py:151` — `NameError` in `getExonFromFile2()`. Variable name case mismatch + undefined `txEnd`. Dead code.
7. **`input_file` undefined** — `scripts/bam_stat.py:67` — Should be `options.input_file`. Marked with `# noqa: F821`.
8. **`urllib.urlopen`** — `rseqc/ireader.py:31` — Python 2 API, should be `urllib.request.urlopen`. URL-based file reading is broken.
9. **`subtractBed3` no-op guard** — `rseqc/BED.py:3003` — `if chrom not in bitsets1` is always False since we iterate over `bitsets1`. Likely should be `bitsets2`.
10. **`string.atoi`/`string.join`/`string.maketrans`** — ~60 call sites in `rseqc/SAM.py` ParseSAM class. These Python 2 `string` module functions don't exist in Python 3. The entire ParseSAM class is dead code (all scripts use ParseBAM via pysam).
11. **`Hill_number` q=1 bug** — `rseqc/mystat.py` — When qvalue=1, ZeroDivisionError falls through to `shannon_entropy(arg)` where `arg` is a comma-separated string, but `shannon_entropy()` expects a list. Iterating over string characters causes `ValueError`.

## Commands

```bash
# Install dev dependencies (using uv)
uv sync

# Run all tests
uv run pytest

# Run a single test file
uv run pytest tests/test_bed.py

# Run a single test
uv run pytest tests/test_bed.py::test_parse_bed_line -v

# Lint
uv run ruff check .

# Format
uv run ruff format .

# Type check
uv run mypy rseqc/
```

## Architecture

### Library (`rseqc/`)
Core modules imported by the CLI scripts:
- **SAM.py** (4600 lines) — the largest module; BAM/SAM parsing, QC metrics computation, gene model overlap. Most scripts depend on this.
- **BED.py** (2600 lines) — BED format parsing, gene model representation, exon/intron/UTR operations.
- **scbam.py** — single-cell BAM utilities (cell barcode demux, UMI handling).
- Smaller utilities: `annoGene.py`, `bam_cigar.py`, `cigar.py`, `fasta.py`, `fastq.py`, `ireader.py` (transparent gz/bz2 reader), `wiggle.py`, `orf.py`, `mystat.py`, `quantile.py`.

### CLI Scripts (`scripts/`)
33 standalone scripts, each with `optparse`-based CLI and a `main()`. Key ones:
- `infer_experiment.py` — infer RNA-seq strandedness
- `tin.py` — transcript integrity number
- `bam_stat.py` — basic BAM statistics
- `geneBody_coverage.py` — gene body coverage profile
- `junction_annotation.py` / `junction_saturation.py` — splice junction analysis
- `read_distribution.py` — reads over genomic features
- `FPKM_count.py` / `FPKM-UQ.py` — expression quantification
- `sc_*.py` — single-cell QC scripts

### Dependencies
Runtime: `pysam`, `bx-python`, `numpy`, `pyBigWig`, `cython`

### Key Patterns
- Almost all scripts follow: parse args with `optparse` → open BAM with pysam → iterate reads → call into `SAM.py` or `BED.py` for logic → write results + optional matplotlib plots.
- `BED.py` and `SAM.py` are tightly coupled — SAM.py imports BED models extensively.

## CI/CD

GitHub Actions workflows should:
- Run pytest on push/PR against Python 3.10, 3.11, 3.12, 3.13
- Run ruff for linting and formatting checks
- Publish to PyPI on tagged releases

## Development Conventions

- Use `uv` for all dependency management and virtual environments
- Write tests in `tests/` mirroring the source structure (e.g., `tests/test_sam.py` for `rseqc/SAM.py`)
- Test fixtures (small BAM/BED files) go in `tests/fixtures/`
- Do not change behavior without a test proving the old behavior first
- When fixing a bug, add a regression test before applying the fix

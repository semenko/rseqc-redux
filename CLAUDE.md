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

**Tests:** 301 passing. Coverage: utility modules 78–97%, BED.py 13%, SAM.py has integration tests with BAM fixture.

**Lint/Type:** CI green — ruff (0 errors, E741/E712 enabled), mypy (0 errors), ruff format clean. Only E501 (line length) still suppressed.

**What's been modernized:**
- Star imports replaced with explicit imports in all `rseqc/` and `scripts/` files
- Bare `except:` → `except Exception:` across entire codebase (97 instances)
- Unused variables removed (147 instances via ruff F841)
- Type hints added to `cigar.py`, `ireader.py`, `bam_cigar.py`, `mystat.py`, `quantile.py`, `dotProduct.py`, `changePoint.py`, `twoList.py`, `FrameKmer.py`, `orf.py`
- Legacy boilerplate removed: Python 3 version checks (all scripts), `__author__`/`__version__` metadata (44 files), UTF-8 encoding declarations, shebangs from library modules, "converted from python2.7" docstrings, dead code and commented-out prints
- CLI smoke tests for all 33 scripts (`--help` flag)
- All 33 CLI scripts migrated from `optparse` to `argparse`
- Ambiguous variable `l` renamed across 12 files (E741)
- `ParseSAM` and `QCSAM` dead code classes removed from SAM.py (~2,970 lines)
- All 11 known bugs fixed with regression tests (see CHANGES.md)
- Session-scoped pysam-built BAM fixture for integration tests
- CLI integration tests (bam_stat, infer_experiment, read_distribution, read_GC, read_quality with real BAM+BED)
- `qcmodule` backward-compat shim removed; all scripts import directly from `rseqc`
- Syntax modernizations: `class Foo(object)` → `class Foo:`, unnecessary `list()` wrappers removed, Python 2 `__div__` replaced with `/` operator
- `sys.exit()` → `sys.exit(1)` in library error paths (SAM.py, BED.py)
- Dead `open()` call removed in annoGene.py
- `type: ignore` comments cleaned up in changePoint.py
- Python 3.13 compatibility: pysam iterator `ValueError` workaround (29 sites across 7 files)

**What still needs work:**
- E501 (line length) — 930 violations, needs incremental reformatting
- Python 3.14 blocked on pysam and pyBigWig releasing 3.14 wheels
- More SAM.py/BED.py method-level integration tests

## Commands

```bash
# Install dev dependencies (using uv)
uv sync

# Run all tests
uv run pytest

# Run a single test file
uv run pytest tests/test_BED.py

# Run a single test
uv run pytest tests/test_BED.py::test_parse_bed_line -v

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
- **SAM.py** (~2,870 lines) — BAM/SAM parsing via pysam, QC metrics computation, gene model overlap. Most scripts depend on this. Contains only `ParseBAM` class (dead `ParseSAM`/`QCSAM` removed).
- **BED.py** (2600 lines) — BED format parsing, gene model representation, exon/intron/UTR operations.
- **scbam.py** — single-cell BAM utilities (cell barcode demux, UMI handling).
- Smaller utilities: `annoGene.py`, `bam_cigar.py`, `cigar.py`, `fasta.py`, `fastq.py`, `ireader.py` (transparent gz/bz2 reader), `wiggle.py`, `orf.py`, `mystat.py`, `quantile.py`.

### CLI Scripts (`scripts/`)
33 standalone scripts, each with `argparse`-based CLI and a `main()`. Key ones:
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
- Almost all scripts follow: parse args with `argparse` → open BAM with pysam → iterate reads → call into `SAM.py` or `BED.py` for logic → write results + optional matplotlib plots.
- `BED.py` and `SAM.py` are tightly coupled — SAM.py imports BED models extensively.

## CI/CD

GitHub Actions workflows should:
- Run pytest on push/PR against Python 3.10, 3.11, 3.12, 3.13
- Run ruff for linting and formatting checks
- Publish to PyPI on tagged releases

## Development Conventions

- Use `uv` for all dependency management and virtual environments
- Write tests in `tests/` mirroring the source structure (e.g., `tests/test_SAM.py` for `rseqc/SAM.py`)
- Test fixtures (small BAM/BED files) go in `tests/fixtures/`
- Do not change behavior without a test proving the old behavior first
- When fixing a bug, add a regression test before applying the fix

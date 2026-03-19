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

## Current State (as of 2026-03-12)

**Infrastructure:** Done — pyproject.toml, CI (3.10–3.13), PyPI publishing.

**Tests:** 476 passing. SAM.py 41%, BED.py 92%, scbam.py 66%, utility modules 83–100% (bam_cigar 100%, twoList 100%). Script helper functions tested (RPKM_saturation, geneBody_coverage, read_distribution, RNA_fragment_size, tin, FPKM_count, read_hexamer, junction_annotation).

**Lint/Type:** CI green — ruff (0 errors, E741/E712/E501 all enabled), mypy (0 errors), ruff format clean. Type hints on all 33 scripts' `main()` and all helper functions.

**What's been modernized:**
- Star imports replaced with explicit imports in all `rseqc/` and `scripts/` files
- Bare `except:` → `except Exception:` → narrowed to specific types (~85 sites) across entire codebase
- Unused variables removed (147 instances via ruff F841)
- Type hints added to all 12 `rseqc/` library modules and all 33 scripts
- Legacy boilerplate removed: Python 3 version checks (all scripts), `__author__`/`__version__` metadata (44 files), UTF-8 encoding declarations, shebangs from library modules, "converted from python2.7" docstrings, dead code and commented-out prints
- CLI smoke tests for all 33 scripts (`--help` flag)
- All 33 CLI scripts migrated from `optparse` to `argparse`
- Ambiguous variable `l` renamed across 12 files (E741)
- `ParseSAM` and `QCSAM` dead code classes removed from SAM.py (~2,970 lines)
- All 11 known bugs fixed with regression tests (see CHANGES.md)
- Session-scoped pysam-built BAM fixture for integration tests
- CLI integration tests for 32 of 33 scripts (all except FPKM_UQ which requires external htseq-count)
- `qcmodule` backward-compat shim removed; all scripts import directly from `rseqc`
- Syntax modernizations: `class Foo(object)` → `class Foo:`, unnecessary `list()` wrappers removed, Python 2 `__div__` replaced with `/` operator
- `exit(0)`/`sys.exit(0)` → `sys.exit(1)` in ~90 error paths across all scripts and library code
- Dead bare field expressions removed (~70 sites across BED.py, SAM.py, scripts)
- `list(map((lambda ...), ...))` → list comprehensions (~80 sites across all BED-parsing code)
- Python 3.13 compatibility: `_pysam_iter()` helper wraps all pysam BAM iteration (handles ValueError bug)
- `while 1: next(samfile)` anti-pattern replaced with `for` loops across all files
- E501 (line length) fully resolved and enabled in ruff — 0 violations
- 95 bare `open()` calls converted to `with` statements (resource leak prevention)
- `heatmap.py` string concatenation bug fixed in `logging.error()`
- `subprocess.call(..., shell=True)` replaced with `subprocess.run([...], check=False)` (23 sites: Rscript, gzip, wigToBigWig, htseq-count)
- `subprocess.run("rm -rf *.pattern", shell=True)` replaced with `glob.glob()` + `os.unlink()` in `scbam.py`
- `BED.ParseBED` supports context managers (`__enter__`/`__exit__`/`close()`) to fix file handle leaks
- Dead ParseBAM methods removed from SAM.py (5 methods: `calculate_rpkm`, `coverageGeneBody`, `junction_freq`, `shuffle_RPKM`, `fetchAlignments` + `print_bits_as_bed`)
- BED.py pruned from ~2,960 lines to ~267 lines: `CompareBED` class deleted, `intersectBed3()` removed, 14 unused `ParseBED` methods removed
- `pandas` removed as direct dependency (logomaker still pulls it transitively; `fastq.py` has one lazy import for logomaker integration)
- Performance: `bamTowig()` inner loop replaced with numpy array slice ops; `readsNVC()` string-key dict → 2D numpy array; `clipping_profile()`/`insertion_profile()` optimized to single-pass CIGAR loop; `stat()` splice check uses inline `any()`, integer `tid` comparison; `mystat.py` all 6 statistical functions cache `sum(lst)` (O(n²) → O(n)); `bam_cigar.py` `list2str`/`list2longstr` use `"".join()` + module-level tuple
- Shared `cli_common.py` module created with `printlog`, `build_bitsets`, `load_chromsize`, `run_rscript` — consolidated from 14+ scripts
- 14 Rscript try/except blocks consolidated into `run_rscript()` calls
- `pysam.Samfile` → `pysam.AlignmentFile` (16 sites across rseqc/ and scripts/)
- Dead functions removed: `searchit` (split_bam), `normalize` (RPKM_saturation), `kmer_freq_seq` (FrameKmer)
- Dead modules removed: `dotProduct.py`, `annoGene.py`, `fasta.py`, `orf.py`, `quantile.py`, `wiggle.py`, `cigar.py`, `changePoint.py`
- Coverage config fixed: CI now measures `--cov=rseqc --cov=scripts`
- `from __future__ import annotations` added to `tin.py`, `RNA_fragment_size.py`
- Performance: `tin.py` and `geneBody_coverage.py` pileup loops replaced with pysam `count_coverage()` (C-level, ~50-100x faster); `readsNVC()` per-character Python loop replaced with numpy vectorized lookup; `calWigSum()` redundant double-fetch removed; `deletionProfile()` list comprehension → `any()` generator; `shannon_entropy()` vectorized with numpy
- `%`-formatting (~92 sites) and string concatenation (~77 sites) converted to f-strings across all `rseqc/` and `scripts/`; ruff rules `UP031`/`UP032` enabled to prevent regression
- `bam_cigar.fetch_exon()` soft clip bug fixed — S ops no longer advance reference coordinate (was double-counting offset for reads with leading soft clips)
- Pre-compiled regex: `scbam.read_match_type()` (6 patterns), `FrameKmer.seq_generator()`, `SAM.mismatchProfile()` — all moved to module-level compiled constants
- `_pysam_iter()` consolidated into `cli_common.py` (was duplicated in SAM.py and scbam.py); backward-compat re-export kept in SAM.py
- `overlay_bigwig.py`: 8 trivial wrapper functions replaced with lambdas in `_ACTIONS` dict + `_apply()` wrapper
- `geneBody_coverage._printlog()` consolidated into `cli_common.printlog(logfile=)` parameter
- `overlap_length2()` renamed to `overlap_length()` in RNA_fragment_size.py
- `bam_cigar.fetch_exon()` removed — all call sites replaced with pysam's `aligned_read.get_blocks()` (handles =/X CIGAR ops more correctly)
- `scbam.mapping_stat()` temp-file + awk antipattern replaced with in-memory `set()` for unique read counting; removed `subprocess`/`glob`/`os` imports
- `mystat.percentile_list()` manual interpolation loop replaced with `np.percentile()`
- `heatmap.make_heatmap()` removed `install.packages("pheatmap")`, replaced manual subprocess with `run_rscript()`
- BED.py header skipping standardized to tuple-form `startswith(("#", "track", "browser"))` in all 6 methods
- `fastq.fasta_iter()` and `fastq.fastq_iter()` replaced with `pysam.FastxFile` (handles gzip, multiline FASTA, robust FASTQ parsing); removed hand-rolled `_open_file()` helper; bz2 support dropped (unused in practice)
- `FrameKmer.seq_generator()` replaced with `pysam.FastxFile` (handles multiline FASTA, gzip, comment lines natively)
- `scbam.list2str()` and `_CIGAR_CHAR` removed — replaced with pysam's native `aligned_read.cigarstring` property
- `tests/fixtures/mini.fq` fixed: seq/qual length mismatch (was silently accepted by hand-rolled parser, caught by pysam)

**What still needs work:**
- Python 3.14 blocked on pysam and pyBigWig releasing 3.14 wheels
- More SAM.py method-level tests (many methods mix computation with file I/O)
- Standardize logging vs `print(file=sys.stderr)` (430+ print-to-stderr calls)

**Potential future directions (lower priority):**
- `SAM.bamTowig()` writes `.wig` then shells out to `wigToBigWig` — could write `.bw` directly via `pyBigWig` (already a dependency), eliminating the external tool dependency and intermediate file. Significant refactor, not a quick cleanup.
- `SAM.bam2fq()` and `SAM.readsNVC()` use manual `str.maketrans()`/`[::-1]` for reverse complementing. The `readsNVC()` base-counting inner loop is already numpy-vectorized, but the reverse complement step remains manual Python.

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
Core modules imported by the CLI scripts (12 modules):
- **SAM.py** (~2,121 lines) — BAM/SAM parsing via pysam, QC metrics computation, gene model overlap. Most scripts depend on this. Contains `ParseBAM` class.
- **BED.py** (~267 lines) — BED format parsing: `ParseBED` class (7 public methods: `getUTR`, `getExon`, `getTranscriptRanges`, `getCDSExon`, `getIntron`, `getIntergenic`, plus context manager support), and module-level functions `unionBed3`, `subtractBed3`, `tillingBed`.
- **scbam.py** — single-cell BAM utilities (cell barcode demux, UMI handling).
- **cli_common.py** — shared CLI utilities: `_pysam_iter`, `printlog`, `build_bitsets`, `load_chromsize`, `run_rscript`.
- Smaller utilities: `bam_cigar.py`, `fastq.py`, `FrameKmer.py`, `getBamFiles.py`, `heatmap.py`, `ireader.py` (transparent gz/bz2 reader), `mystat.py`, `twoList.py`.

### CLI Scripts (`scripts/`)
33 standalone scripts, each with `argparse`-based CLI and a `main()`. Key ones:
- `infer_experiment.py` — infer RNA-seq strandedness
- `tin.py` — transcript integrity number
- `bam_stat.py` — basic BAM statistics
- `geneBody_coverage.py` — gene body coverage profile
- `junction_annotation.py` / `junction_saturation.py` — splice junction analysis
- `read_distribution.py` — reads over genomic features
- `FPKM_count.py` / `FPKM_UQ.py` — expression quantification
- `sc_*.py` — single-cell QC scripts

### Important: Users only call `scripts/`
Users invoke the CLI scripts as entry points. No external code imports `rseqc/` as a library. This means internal APIs can be freely refactored without breaking downstream users.

### Dependencies
Runtime: `pysam`, `bx-python`, `numpy`, `pyBigWig`, `logomaker`, `matplotlib`

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
- **CHANGES.md**: Add a bullet under `[Unreleased]` for every user-visible change (bug fix, new feature, performance improvement, removed code, changed behavior). When bumping the version for a release, rename `[Unreleased]` to the new version/date and add a fresh empty `[Unreleased]` section at the top. Never leave released changes tagged as `[Unreleased]`.

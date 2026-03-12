# Changelog

All notable changes to this project will be documented in this file.

For the full changelog, see [CHANGES.md](https://github.com/semenko/rseqc-redux/blob/main/CHANGES.md) in the repository root.

## [Unreleased]

### Highlights

- 476 automated tests covering all library modules and 32 of 33 CLI scripts
- Type hints on all 33 CLI scripts and all library modules
- 11 bug fixes with regression tests
- 7 dead library modules removed (annoGene, fasta, orf, quantile, wiggle, cigar, changePoint)
- Dead code classes removed: `ParseSAM`, `QCSAM`, `CompareBED`, 6 unused `ParseBAM` methods
- Performance: numpy vectorization in `bamTowig()`, `readsNVC()`, single-pass CIGAR loops
- All 33 CLI scripts migrated from `optparse` to `argparse`
- `pandas` removed as direct dependency
- Shared `cli_common.py` module consolidates utilities from 14+ scripts
- Python 3.13 compatibility fixes

## [5.0.3] - 2026-03-08

Initial modernized release.

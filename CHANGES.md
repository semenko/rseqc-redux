# Changelog

All notable changes to this project will be documented in this file.

## [5.0.3] - 2026-03-12

### Tests

- Expanded unit test coverage from 409 to 550 tests (+141 tests).
- **bam_cigar.py**: 75% → 100% — 35 edge-case CIGAR tests covering all operation types in every fetch function.
- **quantile.py**: 84% → 100% — boundary conditions (`j < 0`, `j >= n`), `issorted=True`, single-element input.
- **twoList.py**: 81% → 100% — `Division`, `Max`, `Min` functions.
- **fasta.py**: 35% → 91% — `addSeq` duplicate, `countBase`, `cal_entropy`, `revComp` (all seqs), `getUniqSeqs`, `findPattern` (with/without rev comp, specific seqID), `fetchSeq` (from file with BED3/BED6, missing chrom).
- **FrameKmer.py**: 40% → 98% — `seq_generator`, `kmer_freq_file` (basic, multi-seq, min_count, frame), `kmer_ratio` edge cases (coding-only, noncoding-only, both-zero, missing kmer).
- **orf.py**: 35% → 83% — `longest_orf_bed` (+/- strand, no ORF), custom codons, `_reverse_comp` edge cases.
- **getBamFiles.py**: 75% → 93% — single BAM, comment lines, comma-separated input, `printit` flag.
- **ireader.py**: 79% → 95% — pipe commands, stdin/stdout, bz2 files.
- **scripts/ helpers**: 47 new tests for directly-importable helper functions across 8 scripts: `square_error` (RPKM_saturation), `valid_name`/`pearson_moment_coefficient` (geneBody_coverage), `cal_size`/`foundone` (read_distribution), `overlap_length2` (RNA_fragment_size), `uniqify`/`shannon_entropy`/`tin_score` (tin), `build_range` (FPKM_count), `file_exist` (read_hexamer), `generate_bed12`/`generate_interact` (junction_annotation).

### Added

- Type hints on all 33 CLI scripts: `main() -> None` and all helper function signatures.
- `from __future__ import annotations` added to `tin.py`, `RNA_fragment_size.py`, and `FPKM_UQ.py`.

### Changed

- Replaced `list(map(int, ...))` with list comprehensions in 6 scripts (tin, FPKM_count, RNA_fragment_size, geneBody_coverage, geneBody_coverage2).
- Removed unnecessary `list()` wrappers on dict view iteration in `normalize_bigwig.py` (2 sites) and `overlay_bigwig.py` (1 site).
- Replaced `list(map(itemgetter(1), g))` with list comprehension in `normalize_bigwig.py`; removed unused `operator.itemgetter` import.
- Removed commented-out dead code in `split_paired_bam.py` and `RNA_fragment_size.py`.
- `sum([comprehension])` → `sum(generator)` in `tin.py`.
- Fixed 3 bare `sys.exit()` → `sys.exit(1)` in `FPKM_UQ.py`.

### Performance

- **bamTowig()**: Replaced per-position `dict` accumulation with numpy array slice operations (`Fwig[start:end] += 1.0`), eliminating the inner Python loop over every base of every exon of every read. ~50-100x faster for large BAM files.
- **readsNVC()**: Replaced `defaultdict` with string keys (`str(i) + "A"`) with a 2D numpy array (`counts[position][base_index]`), eliminating string allocation per base per read.
- **clipping_profile() / insertion_profile()**: Eliminated `list2longstr()` string expansion (which built a per-base character string from CIGAR) — now works directly from CIGAR tuples with integer op codes. Both SE and PE code paths updated. Further optimized: single-pass CIGAR loop computes read length + op presence simultaneously (was 2-3 separate passes per read). PE paths now write directly to target profile dict instead of building intermediate `clip_positions` list.
- **stat()**: Replaced `fetch_intron()` list allocation (just to check `len() == 0`) with `any(c == 3 for c, _s in cigar)` inline check; replaced `getrname()` string comparison with integer `tid != rnext`. Removed redundant `if mapq >= q_cut` guard (always true after prior `continue`).
- **calWigSum()**: Removed dead `read_id + map_strand` expression and now-unused variable assignments.
- **mystat.py**: Cached `sum(lst)` before loops in all 6 statistical functions (`shannon_entropy`, `shannon_entropy_es`, `shannon_entropy_ht`, `simpson_index`, `simpson_index_es`, `Hill_number`) — was O(n²), now O(n). Simplified `RSS()` to single `sum(generator)` expression.
- **cigar.py**: `sum([list comprehension])` → `sum(generator)` to avoid intermediate list allocation (5 sites).
- **cigar.py / bam_cigar.py**: `list2str()` and `list2longstr()` use `"".join()` with module-level tuple lookup instead of string concatenation with per-call dict construction.
- **wiggle.py**: Replaced regex + `str()` NaN check with `math.isnan()` in all fetch_avg/fetch_sum methods; replaced `len(list(range(st, end)))` with `(end - st)`.
- Removed remaining unnecessary `list()` wrappers on dict view iteration in `SAM.py` (4 sites: `bamTowig`, `calWigSum`, `readDupRate`).

### Added

- Snapshot tests for `readsNVC()`, `clipping_profile()`, `bamTowig()`, and `stat()` splice counting — verify exact output values, not just file existence.
- Tests for `wiggle.ParseWig2` sum/avg methods (`fetch_avg_scores_by_range`, `fetch_avg_scores_by_positions`, `fetch_sum_scores_by_range`, `fetch_sum_scores_by_positions`).
- Test for `mystat.shannon_entropy_es` with multi-element input.
- CLI integration tests for 17 additional scripts: `bam2fq`, `bam2wig`, `divide_bam`, `split_bam`, `split_paired_bam`, `FPKM_count`, `RNA_fragment_size`, `RPKM_saturation`, `read_hexamer`, `sc_seqLogo`, `sc_seqQual`, `sc_bamStat`, `sc_editMatrix`, `geneBody_coverage2`, `normalize_bigwig`, `overlay_bigwig`. All 32 of 33 scripts now have integration tests (only `FPKM_UQ` excluded — requires external `htseq-count`).
- Test fixtures: `mini.chrom.sizes`, `mini.fa`, `mini.fq`, session-scoped `mini_bigwig` BigWig fixture.

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
- Narrowed ~85 broad `except Exception:` to specific types (`OSError`, `KeyError`, `IndexError`, `ValueError`, `ZeroDivisionError`, etc.) across all `rseqc/` and `scripts/` files.
- Replaced 23 `subprocess.call(..., shell=True)` with `subprocess.run([...], check=False)` for Rscript (17), gzip (3), wigToBigWig (3), and htseq-count (1) calls.
- Replaced `subprocess.run("rm -rf *.pattern", shell=True)` with `glob.glob()` + `os.unlink()` in `scbam.py`.
- Added `__enter__`/`__exit__`/`close()` context manager support to `BED.ParseBED` and `BED.CompareBED` to fix file handle resource leaks.
- Added 6 CLI integration tests: `read_duplication`, `clipping_profile`, `insertion_profile`, `read_NVC`, `tin`, `geneBody_coverage`.

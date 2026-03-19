# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Fixed

- **bam_cigar.py**: Fix `fetch_exon()` soft clipping bug — soft clip CIGAR ops (S/4) incorrectly advanced the reference coordinate, shifting exon boundaries rightward for reads with leading soft clips. The sibling function `fetch_intron()` already handled this correctly (the equivalent line was commented out). This is a **breaking change** in output for any script that calls `fetch_exon()` — see details below.

  **What was wrong:** In the BAM specification, soft clips consume query sequence but **not** reference positions. pysam's `reference_start` already points past any leading soft clips, so the S op in `fetch_exon()` was double-counting the offset — shifting all exon block coordinates rightward by the clip length. The block *length* was preserved, only the *position* was wrong.

  **Which reads are affected:** Soft clipping occurs when the aligner (e.g. STAR, HISAT2) cannot place part of a read against the reference — typically due to adapter read-through, low-quality bases at read ends, or short overhangs at exon–exon junctions that don't map uniquely. In typical RNA-seq data, 5–15% of reads carry soft clips; roughly half are leading clips. Only **leading** soft clips (e.g. `5S45M`) triggered the bug; trailing clips (`45M5S`) were harmless because the S op came after all exon blocks were already recorded.

  **Practical impact per caller:**
  | Script / method | What it uses exon blocks for | Effect of the bug | Severity |
  |---|---|---|---|
  | `bam2wig` (`bamTowig`) | Per-base coverage array | Coverage deposited a few bp to the right of the true position; total signal magnitude unchanged | Very low — invisible at typical genome-browser zoom |
  | `RPKM_saturation` (`calWigSum`) | Sum of exon block lengths | Block lengths were correct (only position was shifted) | **None** |
  | `read_duplication` (`readDupRate`) | Position-based duplicate key includes exon boundary string | Reads with different clip lengths but identical alignment could get distinct keys, slightly underestimating position-based duplication | Very low — PCR duplicates almost always share the same CIGAR |
  | `inner_distance` (`mRNA_inner_distance`) | Overlap calculation between read1/read2 exon blocks | Off by at most the clip length (1–10 bp) in the overlap count | Negligible vs. typical insert sizes of 200–500 bp |
  | `RPKM_saturation` (`saturation_RPKM`) | Midpoint of exon blocks → gene interval-tree lookup | Midpoint shifted by clip length; could misassign a read at a gene boundary | Very low — affects reads within ~10 bp of a gene edge |
  | `read_distribution` | Midpoint of exon blocks → feature classification (CDS, UTR, intron, intergenic) | Same midpoint shift; could reclassify a read at a feature boundary | Very low — probability proportional to clip_size / feature_size |
  | `FPKM_count` | Midpoint of exon blocks → gene-level counting | Same midpoint shift as above | Very low |

  **Bottom line:** The bug produced incorrect coordinates, but the practical impact on any single analysis was very small to undetectable. Shifts of 1–10 bp for ~5% of reads do not meaningfully change gene-level counts, coverage profiles, or quality metrics. No biological conclusions from prior RNA-seq experiments should need revisiting.

- **tests/fixtures/mini.fq**: Fix sequence/quality length mismatch in reads 1 and 2 — was silently accepted by the hand-rolled line-counting parser, now caught by pysam's htslib-backed FASTQ validation.

### Changed

- **SAM.py**: `bamTowig()` now writes BigWig files natively via `pyBigWig` instead of shelling out to the external `wigToBigWig` binary — removes the only non-Python binary runtime dependency. The `chrom_file` parameter has been removed (internal API only; no downstream callers). `subprocess` import removed from SAM.py.
- **scbam.py**: Pre-compile 6 regex patterns in `read_match_type()` at module level instead of calling `re.search()` with literal patterns per read.
- **FrameKmer.py**: Move `DNA_pat` regex to module-level `_DNA_PAT` (was re-compiled every call to `seq_generator()`)
- **SAM.py**: Move MD-tag regex to module-level `_MD_PAT` (was re-compiled every call to `mismatchProfile()`)
- Consolidate `_pysam_iter()` into `cli_common.py` — was duplicated identically in `SAM.py` and `scbam.py`; 7 scripts updated to import from `cli_common` (backward-compat re-export kept in `SAM.py`)
- **overlay_bigwig.py**: Replace 8 trivial wrapper functions (Add, Subtract, Product, Division, Average, geometricMean, Max, Min) and `_check_list` with inline lambdas in `_ACTIONS` dict + single `_apply()` validation wrapper (~40 lines removed)
- **geneBody_coverage.py**: Delete local `_printlog()`, use `cli_common.printlog()` with new `logfile=` parameter; simplify `pearson_moment_coefficient()` to list comprehension and fix `int(len/2)` → `len//2`
- **RNA_fragment_size.py**: Rename `overlap_length2()` → `overlap_length()` (the `2` suffix was a historical artifact from a deleted predecessor)
- **SAM.py**, **read_distribution.py**: Replace all `bam_cigar.fetch_exon()` calls (6 sites) with pysam's native `aligned_read.get_blocks()`, which also handles `=`/`X` CIGAR ops correctly.
- **scbam.py**: Replace `mapping_stat()` temp-file + awk antipattern (`subprocess.check_output("awk '!a[$0]++' *.reads_id.txt | wc -l", shell=True)`) with in-memory `set()` for unique read counting. Removes file I/O, subprocess calls, and intermediate file cleanup.
- **mystat.py**: Replace manual `percentile_list()` interpolation loop with `np.percentile(..., method='linear')`.
- **heatmap.py**: Remove `install.packages("pheatmap")` from generated R script (fails silently in non-interactive R sessions); replace manual subprocess call with `run_rscript()` from `cli_common`.
- **BED.py**: Standardize header skipping in `getUTR()`, `getTranscriptRanges()`, `getIntron()` from three separate `startswith` checks to tuple form `startswith(("#", "track", "browser"))` (matching `getExon`, `getCDSExon`, `getIntergenic`).
- **fastq.py**: Rewrite `fasta_iter()` and `fastq_iter()` using `pysam.FastxFile` instead of hand-rolled line-counting parser. Gains native gzip support via htslib, correct multiline FASTA concatenation, and strict FASTQ validation. **Breaking:** bz2-compressed input is no longer supported (was handled by the now-removed `_open_file()` helper; bz2 FASTA/FASTQ is practically unused).
- **FrameKmer.py**: Rewrite `seq_generator()` using `pysam.FastxFile` instead of manual line-by-line parsing. Gains multiline FASTA support, native gzip handling, and proper comment/header handling. DNA pattern filter now applied to full concatenated sequences rather than individual lines.
- **scbam.py**: Replace `list2str(aligned_read.cigar)` with pysam's native `aligned_read.cigarstring` property in `mapping_stat()`.

### Removed

- **bam_cigar.py**: Delete `fetch_exon()` function — replaced by pysam's `aligned_read.get_blocks()` at all call sites.
- **fastq.py**: Delete `_open_file()` helper — replaced by `pysam.FastxFile` (handles plain and gzip natively). bz2 decompression support dropped.
- **scbam.py**: Delete `list2str()` function and `_CIGAR_CHAR` tuple — replaced by pysam's native `aligned_read.cigarstring`.
- **scbam.py**: Remove `subprocess`, `glob`, `os` imports (no longer needed after `mapping_stat()` rewrite).

## [6.1.0] - 2026-03-19

### Fixed

- **geneBody_coverage.py**: Guard `Rcode_write()` against empty dataset — `names[0]` crashed with `IndexError` when no BAM files produced coverage signal.

### Changed

- Delete 3 utility modules, reducing `rseqc/` from 12 to 9 modules:
  - **getBamFiles.py**: `isbamfile()` and `get_bam_files()` moved into `cli_common.py`
  - **ireader.py**: Replaced with private `_open_file()` generator in `fastq.py` (handles plain/gz/bz2; drops unused URL, pipe, and stdin support)
  - **twoList.py**: Array operation functions inlined into `overlay_bigwig.py` with `_ACTIONS` dispatch dict
- **tin.py**: Inline `uniqify()` as `list(dict.fromkeys(...))` (canonical Python 3.7+ order-preserving dedup)
- **overlay_bigwig.py**: `Max`/`Min` now use `numpy.maximum`/`numpy.minimum` instead of `map(max, zip(...))`; `--action` argument uses `choices=` for better error messages
- Convert `%`-formatting (~92 sites) and string concatenation (~77 sites) to f-strings across all `rseqc/` and `scripts/` files; enable ruff rules `UP031`/`UP032` to prevent regression
- **FPKM_count.py**: Replace inline strand rule parsing with shared `_parse_strand_rule()`, deduplicate triple print block, use truthiness checks, simplify `elif`→`else`
- **read_hexamer.py**: Remove custom `file_exist()` function, use `os.path.exists()` directly
- **infer_experiment.py**: Replace manual file existence loop with `validate_files_exist()`
- **RNA_fragment_size.py**: Remove redundant `str()` calls on string literals in header output
- Extract shared `iter_bed12()` BED12 parser into `cli_common.py`, replacing 6 duplicated parsing loops across 5 scripts
- Extract `validate_bam_index()` into `cli_common.py`, consolidating 3 duplicate BAM index checks
- Replace manual file existence checks with `validate_files_exist()` in `geneBody_coverage2.py` and `FPKM_UQ.py`
- **read_distribution.py**: Remove dead `foundone()` function, 10 unused `build_bitsets()` calls, and replace 21-element return tuple with `GeneModelResult` NamedTuple
- **read_distribution.py**: Replace 24 repetitive `subtractBed3()` calls with a loop over intergenic regions
- Convert pysam try/finally to context managers in 5 scripts (`split_paired_bam.py`, `RNA_fragment_size.py`, `geneBody_coverage.py`, `divide_bam.py`, `tin.py`)

### Removed

- **getBamFiles.py**, **ireader.py**, **twoList.py**: Deleted (functionality moved to `cli_common.py`, `fastq.py`, and `overlay_bigwig.py` respectively)
- 7 tests for dead `ireader` features (URL, pipe, stdin, stdout, passthrough)
- 5 `uniqify` tests (function inlined)

## [6.0.1] - 2026-03-19

### Fixed

- **SAM.py**: Fix `int(fields[9] == 1)` operator-precedence bug in `annotate_junction()` and `saturation_junction()` — single-exon genes were never excluded from junction analysis, corrupting known/novel classification (same Bug #3 already fixed in BED.py).
- **scbam.py**: Fix `KeyError` crash in `mapping_stat()` when read lacks `RE` tag — `elif tag_dict[RE_tag]` was outside the `if RE_tag in tag_dict` block.
- **SAM.py**: Fix missing `\n` in `mRNA_inner_distance()` output for unknownChromosome case.
- **sc_bamStat.py**: Fix BAM index validation no-op — `if not (file + ".bai")` always evaluated to False; now uses `os.path.exists()`.
- **BED.py**: Fix `getExon()` and `getCDSExon()` crash (IndexError) on BED files with `#`/`track`/`browser` header lines.
- **SAM.py**: Fix `mismatchProfile()` and `deletionProfile()` for-else bug — "Total reads used" message was written to data file instead of stderr when loop completed without break.
- **SAM.py**: Remove dead `inner_distance = 0` assignment before `continue` in `mRNA_inner_distance()`.
- **BED.py**: Remove identical if/else branches in `getIntron()` (strand check had no effect).
- **FPKM_UQ.py**: Fix `--info` argument `default=5` (int for filename) → `default=None`.
- **split_paired_bam.py**: Replace deprecated `pysam.AlignedRead()` with `pysam.AlignedSegment(header)`.
- **junction_saturation.py**, **RPKM_saturation.py**: Fix typo "samller" → "smaller" in error messages.
- **twoList.py**: Fix `Min()` docstring that incorrectly said "max".
- **SAM.py**: Fix typo "Totoal" → "Total" in 6 user-facing stderr messages (`clipping_profile`, `insertion_profile`).
- **SAM.py**: Fix `if`/`if` → `if`/`elif` for `is_read1`/`is_read2` in `configure_experiment()` — both flags set on a malformed read could silently overwrite `read_id`.

### Changed

- CI: auto-create GitHub releases on tag push.
- **divide_bam.py**, **split_paired_bam.py**: Wrap BAM file handles in try/finally to prevent resource leaks.
- **scbam.py**: Close BAM file handles in `barcode_edits()` and `mapping_stat()` via try/finally blocks.
- **scbam.py**: Replace `defaultdict(dict)` + try/except KeyError with `defaultdict(lambda: defaultdict(int))` in `barcode_edits()`.
- **geneBody_coverage.py**: Close pysam.AlignmentFile in `genebody_coverage()` via try/finally.
- **RNA_fragment_size.py**: Close pysam.AlignmentFile in `main()` via try/finally.
- **normalize_bigwig.py**: Close pyBigWig file handle via try/finally.
- **overlay_bigwig.py**: Close two pyBigWig file handles via try/finally.
- **tin.py**: Remove duplicate `build_bitsets()` — now imports from `cli_common`.
- **geneBody_coverage.py**: Rename local `printlog()` → `_printlog()` to distinguish from `cli_common.printlog` (local version also writes to log.txt).
- **sc_bamStat.py**: Replace single-item `for file in [args.bam_file]` loop with direct validation.
- **SAM.py**: Simplify `int(int(x) / 1000)` → `int(x) // 1000` in `saturation_junction()` R output (7 sites).
- **SAM.py**: Remove dead commented-out `current_pos = self.samfile.tell()` in `configure_experiment()`.

### Performance

- **saturation_junction()**: Replaced O(n) rescan of all unique junctions at each percentile step with incremental counters. Known/unknown counts updated as each splice site observation is added.
- **configure_experiment()**: Fast-path the common single-hit case (read overlaps one gene) to avoid `set()` + `':'.join()` allocation per read.
- **mRNA_inner_distance()**: Replace O(n) list of every exon position with O(k) range-overlap arithmetic for overlapping read pairs.
- **RNA_fragment_size.py**: Replace `len(list(range(...)))` O(n) with `max(0, min(...) - max(...))` O(1) in `overlap_length()`.

## [6.0.0] - 2026-03-12

Major modernization release — first release under the `rseqc-redux` name.

### Tests

- Expanded unit test coverage from 409 to 550 tests (+141 tests), then consolidated during module cleanup (476 tests currently).
- **bam_cigar.py**: 75% → 100% — 35 edge-case CIGAR tests covering all operation types in every fetch function.
- **quantile.py**: 84% → 100% — boundary conditions (`j < 0`, `j >= n`), `issorted=True`, single-element input. _(module later removed)_
- **twoList.py**: 81% → 100% — `Division`, `Max`, `Min` functions.
- **fasta.py**: 35% → 91% — `addSeq` duplicate, `countBase`, `cal_entropy`, `revComp` (all seqs), `getUniqSeqs`, `findPattern` (with/without rev comp, specific seqID), `fetchSeq` (from file with BED3/BED6, missing chrom). _(module later removed)_
- **FrameKmer.py**: 40% → 98% — `seq_generator`, `kmer_freq_file` (basic, multi-seq, min_count, frame), `kmer_ratio` edge cases (coding-only, noncoding-only, both-zero, missing kmer).
- **orf.py**: 35% → 83% — `longest_orf_bed` (+/- strand, no ORF), custom codons, `_reverse_comp` edge cases. _(module later removed)_
- **getBamFiles.py**: 75% → 93% — single BAM, comment lines, comma-separated input, `printit` flag.
- **ireader.py**: 79% → 95% — pipe commands, stdin/stdout, bz2 files.
- **scripts/ helpers**: 47 new tests for directly-importable helper functions across 8 scripts: `square_error` (RPKM_saturation), `valid_name`/`pearson_moment_coefficient` (geneBody_coverage), `cal_size`/`foundone` (read_distribution), `overlap_length2` (RNA_fragment_size), `uniqify`/`shannon_entropy`/`tin_score` (tin), `build_range` (FPKM_count), `file_exist` (read_hexamer), `generate_bed12`/`generate_interact` (junction_annotation).

### Added

- Type hints on all 33 CLI scripts: `main() -> None` and all helper function signatures.
- `from __future__ import annotations` added to `tin.py`, `RNA_fragment_size.py`, and `FPKM_UQ.py`.
- Snapshot tests for `readsNVC()`, `clipping_profile()`, `bamTowig()`, and `stat()` splice counting — verify exact output values, not just file existence.
- CLI integration tests for 17 additional scripts: `bam2fq`, `bam2wig`, `divide_bam`, `split_bam`, `split_paired_bam`, `FPKM_count`, `RNA_fragment_size`, `RPKM_saturation`, `read_hexamer`, `sc_seqLogo`, `sc_seqQual`, `sc_bamStat`, `sc_editMatrix`, `geneBody_coverage2`, `normalize_bigwig`, `overlay_bigwig`. All 32 of 33 scripts now have integration tests (only `FPKM_UQ` excluded — requires external `htseq-count`).
- Test fixtures: `mini.chrom.sizes`, `mini.fa`, `mini.fq`, session-scoped `mini_bigwig` BigWig fixture.

### Removed

- Removed `ParseSAM` and `QCSAM` dead code classes from `rseqc/SAM.py` (~2,970 lines). All scripts use `ParseBAM` via pysam.
- Removed `qcmodule/` backward-compatibility shim — all scripts now import directly from `rseqc`.
- Removed Python 3 version checks from all 26 scripts (dead code on Python >=3.10).
- Removed 7 dead library modules: `annoGene.py`, `fasta.py`, `orf.py`, `quantile.py`, `wiggle.py`, `cigar.py`, `changePoint.py` (functionality either unused or inlined).
- Removed `CompareBED` class and 14 unused `ParseBED` methods from `BED.py`. Removed `intersectBed3()`.
- Removed dead `ParseBAM` methods: `calculate_rpkm`, `coverageGeneBody`, `junction_freq`, `shuffle_RPKM`, `fetchAlignments`, `print_bits_as_bed`.
- Removed dead functions: `searchit` (split_bam), `normalize` (RPKM_saturation), `kmer_freq_seq` (FrameKmer).
- Removed dead modules: `dotProduct.py` (benchmark-only).
- Removed ~70 dead bare field expressions (`fields[N]`, `int(fields[N])`) across BED.py, SAM.py, and scripts.
- Removed `pandas` as a direct dependency (only used transitively via logomaker).

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
- Fixed Python 3.13 compatibility: `_pysam_iter()` helper wraps all pysam BAM iteration to handle `ValueError` bug (PEP 745 / coverage instrumentation).

### Changed

- Replaced `while 1: next(samfile)` anti-pattern with `for` loops across `SAM.py` (20 sites), `scbam.py` (4 sites), and 5 scripts.
- Migrated all 33 CLI scripts from `optparse` to `argparse`.
- Renamed ambiguous variable `l` to descriptive names across 12 files (E741).
- `class ParseBAM(object):` → `class ParseBAM:` (modern Python 3 style).
- E501 (line length > 120) fully resolved and enabled in ruff config — 0 violations remaining.
- Converted 95 bare `open()` calls to `with` statements across `rseqc/` and `scripts/` to prevent resource leaks.
- Replaced ~80 `list(map((lambda ...), ...))` with list comprehensions across all BED-parsing code.
- Added type hints to all `rseqc/` library modules.
- Narrowed ~85 broad `except Exception:` to specific types (`OSError`, `KeyError`, `IndexError`, `ValueError`, `ZeroDivisionError`, etc.) across all `rseqc/` and `scripts/` files.
- Replaced 23 `subprocess.call(..., shell=True)` with `subprocess.run([...], check=False)` for Rscript (17), gzip (3), wigToBigWig (3), and htseq-count (1) calls.
- Replaced `subprocess.run("rm -rf *.pattern", shell=True)` with `glob.glob()` + `os.unlink()` in `scbam.py`.
- Added `__enter__`/`__exit__`/`close()` context manager support to `BED.ParseBED` to fix file handle resource leaks.
- Shared `cli_common.py` module created with `printlog`, `build_bitsets`, `load_chromsize`, `run_rscript` — consolidated from 14+ scripts.

### Performance

- **bamTowig()**: Replaced per-position `dict` accumulation with numpy array slice operations (`Fwig[start:end] += 1.0`), eliminating the inner Python loop over every base of every exon of every read. ~50-100x faster for large BAM files.
- **readsNVC()**: Replaced `defaultdict` with string keys (`str(i) + "A"`) with a 2D numpy array (`counts[position][base_index]`), eliminating string allocation per base per read.
- **clipping_profile() / insertion_profile()**: Eliminated `list2longstr()` string expansion — now works directly from CIGAR tuples with integer op codes. Single-pass CIGAR loop (was 2-3 passes per read). PE paths write directly to target dict (no intermediate list).
- **stat()**: Replaced `fetch_intron()` list allocation with `any(c == 3 for c, _s in cigar)` inline check; replaced `getrname()` string comparison with integer `tid != rnext`.
- **mystat.py**: Cached `sum(lst)` before loops in all 6 statistical functions — was O(n²), now O(n).
- **bam_cigar.py**: `list2str()` and `list2longstr()` use `"".join()` with module-level tuple lookup instead of string concatenation.

## [5.0.3] - 2026-03-12

Initial modernized release — performance optimizations, bug fixes, dead code removal, argparse migration, and test infrastructure. Most changes consolidated into the v6.0.0 release notes above.

## [5.0.2] - 2026-03-11

Project infrastructure: `pyproject.toml`, GitHub Actions CI (Python 3.10–3.13), PyPI publishing workflow, ruff linting/formatting.

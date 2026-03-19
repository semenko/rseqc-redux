# rseqc-redux

[![PyPI](https://img.shields.io/pypi/v/rseqc-redux.svg)][pypi status]
[![Status](https://img.shields.io/pypi/status/rseqc-redux.svg)][pypi status]
[![Python Version](https://img.shields.io/pypi/pyversions/rseqc-redux)][pypi status]
[![License](https://img.shields.io/pypi/l/rseqc-redux)][license]

[![Tests](https://github.com/semenko/rseqc-redux/actions/workflows/ci.yml/badge.svg)][tests]
[![Coverage](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/semenko/rseqc-redux/python-coverage-comment-action-data/endpoint.json)][coverage]
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

[pypi status]: https://pypi.org/project/rseqc-redux/
[tests]: https://github.com/semenko/rseqc-redux/actions?workflow=CI
[coverage]: https://htmlpreview.github.io/?https://github.com/semenko/rseqc-redux/blob/python-coverage-comment-action-data/htmlcov/index.html
[license]: https://github.com/semenko/rseqc-redux/blob/main/LICENSE

**rseqc-redux** is a modernized fork of [RSeQC](http://rseqc.sourceforge.net/) (RNA-seq Quality Control), originally by Liguo Wang. It updates the RSeQC 5.0.1 codebase with modern Python packaging, comprehensive tests, and CI — while preserving the original functionality.

## Requirements

- Python 3.10 or higher
- Indexed BAM files (`.bam` with corresponding `.bam.bai` index files)
- A BED-format gene model / reference annotation

### Dependencies

Core dependencies are automatically installed:
- `pysam` — BAM/SAM file handling
- `bx-python` — Interval indexing and overlap
- `numpy` — Numerical operations
- `pyBigWig` — BigWig file I/O
- `matplotlib` / `logomaker` — Plotting

## Installation

### From PyPI (Recommended)

```bash
pip install rseqc-redux
```

### From Source

```bash
git clone https://github.com/semenko/rseqc-redux.git
cd rseqc-redux
pip install .
```

### Development Installation

```bash
git clone https://github.com/semenko/rseqc-redux.git
cd rseqc-redux
uv sync
```

## Quick Start

```bash
# Basic BAM statistics
bam_stat -i sample.bam

# Infer RNA-seq strandedness
infer_experiment -r gene_model.bed -i sample.bam

# Transcript integrity number
tin -i sample.bam -r gene_model.bed

# Gene body coverage
geneBody_coverage -r gene_model.bed -i sample.bam -o output

# Read distribution over genomic features
read_distribution -r gene_model.bed -i sample.bam

# Junction annotation
junction_annotation -r gene_model.bed -i sample.bam -o output
```

## Available Tools

| Tool | Description |
|------|-------------|
| `bam_stat` | Summarize mapping statistics of a BAM file |
| `bam2fq` | Convert BAM alignments to FASTQ format |
| `bam2wig` | Convert BAM to wiggle/BigWig |
| `divide_bam` | Equally divide BAM file into n parts |
| `split_bam` | Split BAM by chromosome |
| `split_paired_bam` | Split paired-end BAM into two single-end BAMs |
| `infer_experiment` | Infer RNA-seq strandedness |
| `inner_distance` | Inner distance between read pairs |
| `RNA_fragment_size` | Fragment size statistics per gene |
| `tin` | Calculate Transcript Integrity Number |
| `geneBody_coverage` | Gene body coverage profile |
| `geneBody_coverage2` | Gene body coverage from BigWig input |
| `read_distribution` | Reads over genomic features (CDS, UTR, intron, etc.) |
| `read_duplication` | Read duplication rate |
| `read_GC` | GC content of reads |
| `read_NVC` | Nucleotide composition (ACGT) along reads |
| `read_quality` | Per-position quality scores |
| `read_hexamer` | Hexamer frequency analysis |
| `junction_annotation` | Annotate splice junctions |
| `junction_saturation` | Splice junction saturation analysis |
| `FPKM_count` | FPKM expression quantification |
| `FPKM-UQ` | Upper-quartile normalized FPKM |
| `RPKM_saturation` | RPKM saturation analysis |
| `mismatch_profile` | Mismatch profile along reads |
| `insertion_profile` | Insertion profile along reads |
| `deletion_profile` | Deletion profile along reads |
| `clipping_profile` | Clipping profile along reads |
| `normalize_bigwig` | Normalize BigWig signal to fixed wigsum |
| `overlay_bigwig` | Pairwise operations on two BigWig files |
| `sc_bamStat` | Single-cell RNA-seq mapping statistics |
| `sc_editMatrix` | Barcode/UMI error correction heatmaps |
| `sc_seqLogo` | DNA sequence logo from FASTQ/FASTA |
| `sc_seqQual` | Sequencing quality heatmap from FASTQ |

## Performance

Benchmarked against RSeQC 5.0.4 on the official RSeQC test data (chr22 subset, ~430K paired-end reads). See [`benchmark-data/BENCHMARK_REPORT.md`](benchmark-data/BENCHMARK_REPORT.md) for the full report with profiling details.

| Script | Original | Redux | Speedup |
|--------|----------|-------|---------|
| `tin` | 23.7s | 3.4s | **7.1x** |
| `geneBody_coverage` | 19.8s | 3.5s | **5.7x** |
| `bam_stat` | 3.5s | 0.8s | **4.2x** |
| `RPKM_saturation` | 3.9s | 1.4s | **2.8x** |
| `read_NVC` | 2.7s | 1.5s | **1.8x** |
| `clipping_profile` | 1.0s | 0.7s | **1.5x** |
| All other scripts | | | ~1.0x |

Key optimizations: pysam `count_coverage()` replacing manual pileup iteration (tin, geneBody_coverage), numpy vectorization (read_NVC, RPKM_saturation), single-pass CIGAR loops (clipping/insertion profiles), and dead code removal (bam_stat).

## Upgrading from RSeQC 5.x or earlier rseqc-redux releases

### Soft clipping bug fix (post-6.1.0)

`fetch_exon()` in `bam_cigar.py` had a long-standing bug inherited from the original RSeQC: soft clip CIGAR operations incorrectly advanced the reference coordinate, shifting exon boundaries rightward for reads with leading soft clips. The sibling function `fetch_intron()` already handled this correctly. This affects output from any script that resolves read-to-genome exon coordinates: `bam2wig`, `read_distribution`, `read_duplication`, `RPKM_saturation`, `FPKM_count`, and `inner_distance`.

**Should I re-run my analyses?** Almost certainly not. The shifts are small (1–10 bp, matching the soft clip length) and only affect the ~5% of reads with leading soft clips. Gene-level counts, coverage profiles, and quality metrics change by amounts well below meaningful thresholds. See [CHANGES.md](CHANGES.md) for a per-script impact table.

## Contributing

Contributions are welcome!

### Development Commands

```bash
# Install development dependencies
uv sync

# Run tests
uv run pytest

# Run full test matrix
uv run nox

# Lint and format
uv run ruff check .
uv run ruff format .

# Type check
uv run mypy rseqc/
```

## License

Distributed under the terms of the [GPLv3+ license][license], rseqc-redux is free and open source software.

## Issues

If you encounter any problems, please [file an issue] with a description of the problem, steps to reproduce, and relevant error messages.

[file an issue]: https://github.com/semenko/rseqc-redux/issues

## Credits

Original RSeQC by Liguo Wang — [rseqc.sourceforge.net](http://rseqc.sourceforge.net/)

Modernization by [Nick Semenkovich (@semenko)](https://nick.semenkovich.com/).

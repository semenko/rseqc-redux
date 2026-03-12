# rseqc-redux

**Modernized RNA-seq Quality Control** — a fork of [RSeQC](http://rseqc.sourceforge.net/) 5.0.1.

[![PyPI](https://img.shields.io/pypi/v/rseqc-redux)](https://pypi.org/project/rseqc-redux/)
[![Python](https://img.shields.io/pypi/pyversions/rseqc-redux)](https://pypi.org/project/rseqc-redux/)
[![CI](https://github.com/semenko/rseqc-redux/actions/workflows/ci.yml/badge.svg)](https://github.com/semenko/rseqc-redux/actions/workflows/ci.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Quick start

```bash
pip install rseqc-redux
```

```bash
bam_stat.py -i sample.bam
```

## What is rseqc-redux?

rseqc-redux provides **33 command-line tools** for comprehensive quality control of RNA-seq data. It evaluates sequence quality, nucleotide composition, GC content, read distribution over genomic features, junction saturation, gene body coverage, and much more.

This is a modernized fork of [RSeQC](http://rseqc.sourceforge.net/) by Liguo Wang, with:

- Python 3.10+ support (tested on 3.10–3.13)
- 470+ automated tests
- Bug fixes and correctness improvements
- Modern packaging with `pyproject.toml`

## Tool categories

| Category | Tools |
|----------|-------|
| **BAM Statistics** | [bam_stat](tools/bam_stat.md), [bam2fq](tools/bam2fq.md), [bam2wig](tools/bam2wig.md), [divide_bam](tools/divide_bam.md), [split_bam](tools/split_bam.md), [split_paired_bam](tools/split_paired_bam.md) |
| **Read Quality** | [read_quality](tools/read_quality.md), [read_GC](tools/read_GC.md), [read_NVC](tools/read_NVC.md), [read_hexamer](tools/read_hexamer.md), [read_duplication](tools/read_duplication.md), [read_distribution](tools/read_distribution.md) |
| **Alignment Profiles** | [clipping_profile](tools/clipping_profile.md), [deletion_profile](tools/deletion_profile.md), [insertion_profile](tools/insertion_profile.md), [mismatch_profile](tools/mismatch_profile.md) |
| **Gene Body & Expression** | [geneBody_coverage](tools/geneBody_coverage.md), [geneBody_coverage2](tools/geneBody_coverage2.md), [FPKM_count](tools/FPKM_count.md), [FPKM-UQ](tools/FPKM_UQ.md), [RPKM_saturation](tools/RPKM_saturation.md), [tin](tools/tin.md) |
| **Junction Analysis** | [junction_annotation](tools/junction_annotation.md), [junction_saturation](tools/junction_saturation.md) |
| **Fragment & Strandedness** | [infer_experiment](tools/infer_experiment.md), [inner_distance](tools/inner_distance.md), [RNA_fragment_size](tools/RNA_fragment_size.md) |
| **BigWig Utilities** | [normalize_bigwig](tools/normalize_bigwig.md), [overlay_bigwig](tools/overlay_bigwig.md) |
| **Single Cell** | [sc_bamStat](tools/sc_bamStat.md), [sc_editMatrix](tools/sc_editMatrix.md), [sc_seqLogo](tools/sc_seqLogo.md), [sc_seqQual](tools/sc_seqQual.md) |

## Credits

Originally developed by [Liguo Wang](mailto:wangliguo78@gmail.com) as [RSeQC](http://rseqc.sourceforge.net/). Modernized and maintained by [Nick Semenkovich](https://github.com/semenko).

# rseqc-redux

[![CI](https://github.com/semenko/rseqc-redux/actions/workflows/ci.yml/badge.svg)](https://github.com/semenko/rseqc-redux/actions/workflows/ci.yml)

A modernized fork of [RSeQC](http://rseqc.sourceforge.net/) (RNA-seq Quality Control), originally by Liguo Wang.

rseqc-redux updates the RSeQC 5.0.1 codebase with modern Python packaging, comprehensive tests, and CI — while preserving the original functionality.

## Installation

```bash
pip install rseqc-redux
# or
uv add rseqc-redux
```

## Usage

```bash
# Basic BAM statistics
bam_stat -i sample.bam

# Infer RNA-seq strandedness
infer_experiment -r gene_model.bed -i sample.bam

# Transcript integrity number
tin -i sample.bam -r gene_model.bed
```

## Development

```bash
# Clone and install
git clone https://github.com/semenko/rseqc-redux.git
cd rseqc-redux
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

GPLv2 — see [LICENSE](LICENSE) for details.

## Credits

Original RSeQC by Liguo Wang — [rseqc.sourceforge.net](http://rseqc.sourceforge.net/)

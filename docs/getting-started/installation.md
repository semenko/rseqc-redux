# Installation

## Prerequisites

- Python 3.10 or later
- A C compiler (required by pysam and bx-python)

## Install from PyPI

```bash
pip install rseqc-redux
```

Or with [uv](https://docs.astral.sh/uv/):

```bash
uv pip install rseqc-redux
```

## Install from source

```bash
git clone https://github.com/semenko/rseqc-redux.git
cd rseqc-redux
pip install .
```

## Verify installation

```bash
bam_stat.py --version
```

## System dependencies

Some dependencies require system libraries:

- **pysam** requires `htslib` headers. On most systems, pip will build from source automatically. If you encounter build errors, install `htslib` via your package manager:
    - macOS: `brew install htslib`
    - Ubuntu/Debian: `apt install libhts-dev`
- **pyBigWig** requires `libcurl` and `zlib` headers.

## Optional: gene model files

Many tools require a reference gene model in BED12 format. See [Input Formats](input-formats.md) for details on obtaining these files.

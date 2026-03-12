# read_quality.py

Calculate Phred quality score distribution for each position on the read.

!!! note
    Each read should have the same (fixed) length.

## Usage

```bash
read_quality.py -i input.bam -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in BAM or SAM format | Required |
| `-o`, `--out-prefix` | Prefix of output files | Required |
| `-r`, `--reduce` | Ignore nucleotides with a particular phred score occurring fewer than this many times (reduces R vector size) | 1 |
| `-q`, `--mapq` | Minimum mapping quality for "uniquely mapped" | 30 |

## Output

- `prefix.qual.r` — R script for generating boxplot
- `prefix.qual.boxplot.pdf` — quality boxplot (if R is available)

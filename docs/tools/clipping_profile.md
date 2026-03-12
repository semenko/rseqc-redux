# clipping_profile.py

Estimate the soft-clipping profile of RNA-seq reads. CIGAR strings must contain 'S' operations (the aligner must support clipped mapping).

## Usage

```bash
clipping_profile.py -i input.bam -o output_prefix -s SE
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in BAM or SAM format | Required |
| `-o`, `--out-prefix` | Prefix of output files | Required |
| `-s`, `--sequencing` | Sequencing layout: `SE` (single-end) or `PE` (paired-end) | Required |
| `-q`, `--mapq` | Minimum mapping quality for "uniquely mapped" | 30 |

## Output

- `prefix.clipping_profile.xls` — clipping percentages per read position
- `prefix.clipping_profile.r` — R script for clipping profile plot
- `prefix.clipping_profile.pdf` — clipping profile plot (if R is available)

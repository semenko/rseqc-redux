# mismatch_profile.py

Calculate the distribution of mismatches across read positions.

!!! note
    The `MD` tag must be present in the BAM file.

## Usage

```bash
mismatch_profile.py -i input.bam -l 150 -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input` | Input BAM file | Required |
| `-l`, `--read-align-length` | Alignment length of reads (usually the original read length) | Required |
| `-o`, `--out-prefix` | Prefix of output files | Required |
| `-n`, `--read-num` | Number of aligned reads with mismatches to use | 1000000 |
| `-q`, `--mapq` | Minimum mapping quality | 30 |

## Output

- `prefix.mismatch_profile.txt` — mismatch percentages per read position
- `prefix.mismatch_profile.r` — R script for mismatch profile plot
- `prefix.mismatch_profile.pdf` — mismatch profile plot (if R is available)

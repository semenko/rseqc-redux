# deletion_profile.py

Calculate the distribution of deleted nucleotides across read positions.

## Usage

```bash
deletion_profile.py -i input.bam -l 150 -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input` | Input BAM file | Required |
| `-l`, `--read-align-length` | Alignment length of reads (usually the original read length) | Required |
| `-o`, `--out-prefix` | Prefix of output files | Required |
| `-n`, `--read-num` | Number of aligned reads with deletions to use | 1000000 |
| `-q`, `--mapq` | Minimum mapping quality | 30 |

## Output

- `prefix.deletion_profile.txt` — deletion percentages per read position
- `prefix.deletion_profile.r` — R script for deletion profile plot
- `prefix.deletion_profile.pdf` — deletion profile plot (if R is available)

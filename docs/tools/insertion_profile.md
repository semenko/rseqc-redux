# insertion_profile.py

Calculate the distribution of inserted nucleotides across read positions. CIGAR strings must contain 'I' operations.

## Usage

```bash
insertion_profile.py -i input.bam -o output_prefix -s SE
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in BAM or SAM format | Required |
| `-o`, `--out-prefix` | Prefix of output files | Required |
| `-s`, `--sequencing` | Sequencing layout: `SE` (single-end) or `PE` (paired-end) | Required |
| `-q`, `--mapq` | Minimum mapping quality for "uniquely mapped" | 30 |

## Output

- `prefix.insertion_profile.xls` — insertion percentages per read position
- `prefix.insertion_profile.r` — R script for insertion profile plot
- `prefix.insertion_profile.pdf` — insertion profile plot (if R is available)

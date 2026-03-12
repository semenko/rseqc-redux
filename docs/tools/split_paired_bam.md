# split_paired_bam.py

Split a paired-end BAM file into two single-end BAM files.

## Usage

```bash
split_paired_bam.py -i input.bam -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in BAM or SAM format (sorted, indexed) | Required |
| `-o`, `--out-prefix` | Prefix of output BAM files | Required |

## Output

- `prefix.R1.bam` — first reads (read 1)
- `prefix.R2.bam` — second reads (read 2)

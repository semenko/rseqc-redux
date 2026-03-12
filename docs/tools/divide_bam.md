# divide_bam.py

Equally divide a BAM file into n parts. Each part contains roughly m/n alignments randomly sampled from the total.

## Usage

```bash
divide_bam.py -i input.bam -n 4 -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in BAM format (sorted, indexed) | Required |
| `-n`, `--subset-num` | Number of output BAM files | Required |
| `-o`, `--out-prefix` | Prefix of output BAM files (produces `prefix_1.bam`, etc.) | Required |
| `-s`, `--skip-unmap` | Skip unmapped reads | Off |

## Output

Produces n BAM files named `prefix_1.bam`, `prefix_2.bam`, ..., `prefix_n.bam`.

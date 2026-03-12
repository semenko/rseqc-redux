# read_NVC.py

Check nucleotide frequency at each position of the read (5' to 3'). Generates a Nucleotide Versus Cycle (NVC) plot.

## Usage

```bash
read_NVC.py -i input.bam -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in BAM or SAM format | Required |
| `-o`, `--out-prefix` | Prefix of output files | Required |
| `-x`, `--nx` | Include N and X nucleotides in the NVC plot | Off |
| `-q`, `--mapq` | Minimum mapping quality for "uniquely mapped" | 30 |

## Output

- `prefix.NVC.xls` — tab-separated nucleotide frequencies per position
- `prefix.NVC_plot.r` — R script for NVC plot
- `prefix.NVC_plot.pdf` — NVC plot (if R is available)

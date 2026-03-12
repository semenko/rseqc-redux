# read_GC.py

Calculate the GC content distribution of mapped reads.

## Usage

```bash
read_GC.py -i input.bam -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in BAM or SAM format | Required |
| `-o`, `--out-prefix` | Prefix of output files | Required |
| `-q`, `--mapq` | Minimum mapping quality for "uniquely mapped" | 30 |

## Output

- `prefix.GC.xls` — tab-separated: GC%, read count
- `prefix.GC_plot.r` — R script for GC distribution plot
- `prefix.GC_plot.pdf` — GC distribution plot (if R is available)

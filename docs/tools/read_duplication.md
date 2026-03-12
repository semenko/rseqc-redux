# read_duplication.py

Calculate reads duplication rate using two methods:

- **Sequence-based**: reads with identical sequence are duplicates
- **Mapping-based**: reads mapped to the exact same genomic location are duplicates

## Usage

```bash
read_duplication.py -i input.bam -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in BAM or SAM format | Required |
| `-o`, `--out-prefix` | Prefix of output files | Required |
| `-u`, `--up-limit` | Upper limit of read occurrence for plotting | 500 |
| `-q`, `--mapq` | Minimum mapping quality for "uniquely mapped" | 30 |

## Output

- `prefix.pos.DupRate.xls` — position-based duplication rate
- `prefix.seq.DupRate.xls` — sequence-based duplication rate
- `prefix.DupRate_plot.r` — R script for duplication plot
- `prefix.DupRate_plot.pdf` — duplication rate plot (if R is available)

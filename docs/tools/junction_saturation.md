# junction_saturation.py

Assess whether sequencing depth is sufficient for junction discovery. As sequencing depth approaches saturation, fewer new junctions are discovered.

## Usage

```bash
junction_saturation.py -i input.bam -r gene_model.bed -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in BAM or SAM format | Required |
| `-o`, `--out-prefix` | Prefix of output files | Required |
| `-r`, `--refgene` | Reference gene model in BED format | Required |
| `-l`, `--percentile-floor` | Sampling starts from this percentile | 5 |
| `-u`, `--percentile-ceiling` | Sampling ends at this percentile | 100 |
| `-s`, `--percentile-step` | Sampling step size | 5 |
| `-m`, `--min-intron` | Minimum intron size (bp) | 50 |
| `-v`, `--min-coverage` | Minimum supporting reads to call a junction | 1 |
| `-q`, `--mapq` | Minimum mapping quality for "uniquely mapped" | 30 |

## Output

- `prefix.junctionSaturation_plot.r` — R script for saturation plot
- `prefix.junctionSaturation_plot.pdf` — saturation curve (known, novel, all junctions)

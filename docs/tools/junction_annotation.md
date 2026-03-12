# junction_annotation.py

Annotate splicing junctions against a gene model at two levels:

- **Read level**: individual spliced reads
- **Junction level**: consolidated junctions from multiple reads spanning the same intron

## Usage

```bash
junction_annotation.py -i input.bam -r gene_model.bed -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in BAM or SAM format | Required |
| `-r`, `--refgene` | Reference gene model in BED format | Required |
| `-o`, `--out-prefix` | Prefix of output files | Required |
| `-m`, `--min-intron` | Minimum intron length (bp) | 50 |
| `-q`, `--mapq` | Minimum mapping quality for "uniquely mapped" | 30 |

## Output

- `prefix.junction.xls` — annotated junctions (known, novel, partial novel)
- `prefix.junction_plot.r` — R script for junction plots
- `prefix.splice_events.pdf` — splice event pie chart
- `prefix.splice_junction.pdf` — splice junction pie chart

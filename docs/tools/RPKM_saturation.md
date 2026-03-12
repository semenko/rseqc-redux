# RPKM_saturation.py

Assess whether RPKM values are saturated by re-sampling reads at increasing percentages. Strand-specific protocols are supported.

## Usage

```bash
RPKM_saturation.py -i input.bam -r gene_model.bed -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in BAM or SAM format | Required |
| `-o`, `--out-prefix` | Prefix of output files | Required |
| `-r`, `--refgene` | Reference gene model in BED format | Required |
| `-d`, `--strand` | Strand rule: `1++,1--,2+-,2-+` or `1+-,1-+,2++,2--` or `none` | None |
| `-l`, `--percentile-floor` | Sampling starts from this percentile | 5 |
| `-u`, `--percentile-ceiling` | Sampling ends at this percentile | 100 |
| `-s`, `--percentile-step` | Sampling step size | 5 |
| `-c`, `--rpkm-cutoff` | Ignore transcripts with RPKM below this value | 0.01 |
| `-q`, `--mapq` | Minimum mapping quality for "uniquely mapped" | 30 |

## Output

- `prefix.saturation.r` — R script for saturation plot
- `prefix.saturation.pdf` — saturation plot (if R is available)
- `prefix.rawCount.xls` — raw count at each sampling level
- `prefix.eRPKM.xls` — estimated RPKM values

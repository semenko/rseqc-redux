# tin.py

Calculate the Transcript Integrity Number (TIN) for each transcript in a BED file. TIN is conceptually similar to RIN (RNA Integrity Number) but provides transcript-level measurement of RNA quality and is more sensitive to RNA degradation.

## Usage

```bash
tin.py -i input.bam -r gene_model.bed
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input` | Input BAM file(s). Accepts a single BAM, comma-separated BAMs, or a directory of BAM files. All must be sorted and indexed | Required |
| `-r`, `--refgene` | Reference gene model in BED format (standard 12-column) | Required |
| `-c`, `--minCov` | Minimum number of reads mapped to a transcript | 10 |
| `-n`, `--sample-size` | Number of equally-spaced positions sampled from each mRNA | 100 |
| `-s`, `--subtract-background` | Subtract background noise estimated from intronic reads | Off |

## Output

- Tab-separated file to stdout with columns: geneID, chrom, tx_start, tx_end, TIN
- Summary statistics printed to stderr (median TIN)

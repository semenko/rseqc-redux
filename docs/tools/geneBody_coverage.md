# geneBody_coverage.py

Calculate RNA-seq read coverage over the gene body. This is useful for assessing RNA integrity and 3'/5' bias.

!!! note
    - Only sorted and indexed BAM files are supported (not SAM)
    - Genes/transcripts with mRNA length < 100 bp are skipped

## Usage

```bash
geneBody_coverage.py -i input.bam -r gene_model.bed -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input` | Input BAM file(s). Accepts a single BAM, comma-separated BAMs, or a text file listing BAM paths | Required |
| `-r`, `--refgene` | Reference gene model in BED format | Required |
| `-o`, `--out-prefix` | Prefix of output files | Required |
| `-l`, `--minimum_length` | Minimum mRNA length (bp) | 100 |
| `-f`, `--format` | Output figure format: `pdf`, `png`, or `jpeg` | pdf |

## Output

- `prefix.geneBodyCoverage.curves.pdf` — gene body coverage plot
- `prefix.geneBodyCoverage.txt` — coverage percentiles (0–100)
- `prefix.geneBodyCoverage.r` — R script for coverage plot

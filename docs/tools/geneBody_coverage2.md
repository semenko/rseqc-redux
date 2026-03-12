# geneBody_coverage2.py

Calculate RNA-seq read coverage over the gene body using a BigWig file as input (instead of BAM).

## Usage

```bash
geneBody_coverage2.py -i input.bw -r gene_model.bed -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Coverage signal file in BigWig format | Required |
| `-r`, `--refgene` | Reference gene model in BED format | Required |
| `-o`, `--out-prefix` | Prefix of output files | Required |
| `-t`, `--graph-type` | Output format: `pdf`, `jpeg`, `bmp`, `tiff`, or `png` | pdf |

## Output

- `prefix.geneBodyCoverage.curves.pdf` — gene body coverage plot
- `prefix.geneBodyCoverage.txt` — coverage percentiles

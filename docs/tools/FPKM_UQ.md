# FPKM-UQ.py

Calculate count, FPKM, and FPKM-UQ (upper quartile normalized FPKM) values per gene, following the [GDC mRNA pipeline definition](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#upper-quartile-fpkm).

## Usage

```bash
FPKM-UQ.py --bam input.bam --gtf genes.gtf -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `--bam` | Alignment file in BAM format (sorted, indexed) | Required |
| `--gtf` | Gene model in GTF format | Required |
| `--info` | Gene model information file | 5 |
| `-o`, `--output` | Prefix of output file | Required |
| `--log2` | Convert FPKM and FPKM-UQ to log2(x+1) scale | Off |

## Output

Tab-separated file with count, FPKM, and FPKM-UQ values per gene.

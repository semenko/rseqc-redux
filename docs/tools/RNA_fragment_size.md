# RNA_fragment_size.py

Calculate fragment size statistics for each gene/transcript. Reports the number of fragments used, plus mean, median, and standard deviation of fragment size.

## Usage

```bash
RNA_fragment_size.py -i input.bam -r gene_model.bed
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input` | Input BAM file | Required |
| `-r`, `--refgene` | Reference gene model in BED format (standard 12-column) | Required |
| `-q`, `--mapq` | Minimum mapping quality for "uniquely mapped" | 30 |
| `-n`, `--frag-num` | Minimum number of fragments per gene | 3 |

## Output

Tab-separated output to stdout with columns: gene, chrom, start, end, fragment_count, mean_size, median_size, stdev_size.

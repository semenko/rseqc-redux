# FPKM_count.py

Calculate raw read count, FPM (fragments per million), and FPKM (fragments per million mapped reads per kilobase exon) for each gene in a BED file.

## Usage

```bash
FPKM_count.py -i input.bam -r gene_model.bed -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in BAM format (SAM not supported) | Required |
| `-o`, `--out-prefix` | Prefix of output files | Required |
| `-r`, `--refgene` | Reference gene model in BED format | Required |
| `-d`, `--strand` | Strand rule: `1++,1--,2+-,2-+` or `1+-,1-+,2++,2--` or `none` | None |
| `-u`, `--skip-multi-hits` | Skip multi-hit reads | Off |
| `-e`, `--only-exonic` | Only count reads falling within exons | Off |
| `-q`, `--mapq` | Minimum mapping quality for "uniquely mapped" | 30 |
| `-s`, `--single-read` | Weight for read pairs with only one end mapped (0–1) | 1 |

## Output

Tab-separated file with columns: chrom, start, end, name, score, strand, raw_count, FPM, FPKM.

# split_bam.py

Split a BAM file according to an input gene list (BED format).

## Usage

```bash
split_bam.py -i input.bam -r genes.bed -o output_prefix
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in BAM or SAM format (sorted, indexed) | Required |
| `-r`, `--genelist` | Gene list in BED format | Required |
| `-o`, `--out-prefix` | Prefix of output BAM files | Required |

## Output

- `prefix.in.bam` — reads mapped to exon regions defined by the gene list
- `prefix.ex.bam` — reads not mapped to those regions
- `prefix.junk.bam` — unmapped or QC-failed reads

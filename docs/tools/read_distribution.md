# read_distribution.py

Calculate how mapped reads are distributed over genomic features: exons, introns, UTRs, intergenic regions, etc.

The following reads are skipped: QC-failed, PCR duplicates, unmapped, non-primary (secondary).

## Usage

```bash
read_distribution.py -i input.bam -r gene_model.bed
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i`, `--input-file` | Alignment file in BAM or SAM format | Required |
| `-r`, `--refgene` | Reference gene model in BED format | Required |

## Output

Prints a table to stdout with read counts and tags for each genomic feature:

```
Group               Total_Tags         Tag_count          Tags/Kb
CDS_Exons           ...                ...                ...
5'UTR_Exons         ...                ...                ...
3'UTR_Exons         ...                ...                ...
Introns             ...                ...                ...
TSS_up_10kb         ...                ...                ...
TES_down_10kb       ...                ...                ...
```

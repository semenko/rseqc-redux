# Input Formats

## BAM files

Most tools require a sorted, indexed BAM file as input.

```bash
# Sort and index with samtools
samtools sort input.bam -o sorted.bam
samtools index sorted.bam
```

!!! note
    SAM format is accepted by some tools, but BAM is recommended for performance.

## BED12 gene models

Many tools require a reference gene model in **BED12 format** (12-column BED). Each line represents a transcript with exon structure.

The 12 columns are:

| Column | Description |
|--------|-------------|
| chrom | Chromosome |
| chromStart | Transcript start |
| chromEnd | Transcript end |
| name | Transcript name |
| score | Score (0-1000) |
| strand | + or - |
| thickStart | CDS start |
| thickEnd | CDS end |
| itemRgb | Color (usually 0) |
| blockCount | Number of exons |
| blockSizes | Comma-separated exon sizes |
| blockStarts | Comma-separated exon starts (relative to chromStart) |

### Where to download gene models

- **UCSC Table Browser**: [genome.ucsc.edu/cgi-bin/hgTables](https://genome.ucsc.edu/cgi-bin/hgTables) — select "BED" as output format, check "Whole Gene"
- **GENCODE**: Download GTF and convert with `gtfToGenePred` and `genePredToBed` from UCSC tools

!!! tip
    Use a comprehensive gene model (e.g., GENCODE basic) that pools all transcripts. This improves junction and exon coverage analysis.

## Chromosome size files

Required by `bam2wig.py`. A tab-separated file with chromosome name and size:

```
chr1	248956422
chr2	242193529
chr3	198295559
...
```

Generate from a BAM file:

```bash
samtools view -H input.bam | grep @SQ | awk '{print $2"\t"$3}' | sed 's/SN://;s/LN://' > chrom.sizes
```

Or download from UCSC:

```bash
fetchChromSizes hg38 > hg38.chrom.sizes
```

## BigWig files

Required by `geneBody_coverage2.py`, `normalize_bigwig.py`, and `overlay_bigwig.py`. BigWig is a compressed binary format for continuous signal data. Generate from a wiggle file:

```bash
wigToBigWig input.wig chrom.sizes output.bw
```
